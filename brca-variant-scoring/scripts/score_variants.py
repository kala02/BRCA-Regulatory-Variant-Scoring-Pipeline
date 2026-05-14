# scripts/score_variants.py

import pandas as pd
import sys
import os
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# ── 1. Read arguments passed in from Snakemake ──────────────────────────────
# snakemake.input[0]  → path to your VCF
# snakemake.output[0] → path where this script should save its result
# snakemake.config    → the config dict from Snakefile (holds API key)

vcf_path    = snakemake.input[0]
output_path = snakemake.output[0]
api_key     = snakemake.config['api_key']

# ── 2. Connect to AlphaGenome ────────────────────────────────────────────────
print("Connecting to AlphaGenome...")
dna_model = dna_client.create(api_key)

# ── 3. Load your VCF ────────────────────────────────────────────────────────
vcf = pd.read_csv(vcf_path, sep='\t')
print(f"Loaded {len(vcf)} variants from {vcf_path}")

# ── 4. Choose scorers — RNA-seq, ATAC, and histone marks ────────────────────
# These are the three most interpretable output types for regulatory analysis
selected_scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['ATAC'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE'],
]

# ── 5. Score every variant in a loop ────────────────────────────────────────
all_results = []

for i, row in vcf.iterrows():
    print(f"  Scoring {row.variant_id} ({i+1}/{len(vcf)})...")

    # Build the Variant object from the VCF row
    variant = genome.Variant(
        chromosome=str(row.CHROM),
        position=int(row.POS),
        reference_bases=str(row.REF),
        alternate_bases=str(row.ALT),
        name=str(row.variant_id),
    )

    # Resize to 1MB window centered on the variant
    interval = variant.reference_interval.resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )

    # Call AlphaGenome
    scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=selected_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS,
    )

    # Convert to tidy DataFrame and tag with metadata from VCF
    df = variant_scorers.tidy_scores([scores])
    df['gene']         = row.GENE
    df['significance'] = row.SIGNIFICANCE
    all_results.append(df)

# ── 6. Combine all results and save ─────────────────────────────────────────
final_df = pd.concat(all_results, ignore_index=True)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
final_df.to_csv(output_path, index=False)

print(f"\n✅ Saved {len(final_df):,} rows to {output_path}")
print(f"   Columns: {list(final_df.columns)}")