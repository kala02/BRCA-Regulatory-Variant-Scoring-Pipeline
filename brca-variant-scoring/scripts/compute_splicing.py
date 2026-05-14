# scripts/compute_splicing.py

import pandas as pd
import os
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

vcf_path    = snakemake.input[0]
output_path = snakemake.output[0]
api_key     = snakemake.config['api_key']

print("Connecting to AlphaGenome...")
dna_model = dna_client.create(api_key)

vcf = pd.read_csv(vcf_path, sep='\t')

splicing_scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITE_USAGE'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_JUNCTIONS'],
]

all_results = []

for i, row in vcf.iterrows():
    print(f"  Splicing scoring {row.variant_id} ({i+1}/{len(vcf)})...")

    variant = genome.Variant(
        chromosome=str(row.CHROM),
        position=int(row.POS),
        reference_bases=str(row.REF),
        alternate_bases=str(row.ALT),
        name=str(row.variant_id),
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

    scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=splicing_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS,
    )

    df = variant_scorers.tidy_scores([scores])
    df['gene']         = row.GENE
    df['significance'] = row.SIGNIFICANCE
    all_results.append(df)

# ── Combine all results ──────────────────────────────────────────────────────
raw_df = pd.concat(all_results, ignore_index=True)

# FIX: convert any non-string columns that could contain objects to strings
# The 'scored_interval' or other columns may contain Variant/Interval objects
for col in raw_df.columns:
    if raw_df[col].dtype == object:
        raw_df[col] = raw_df[col].astype(str)

# ── Compute merged splicing score ────────────────────────────────────────────
merged = (
    raw_df.groupby(['variant_id', 'gene', 'significance', 'output_type'])['raw_score']
    .max()
    .reset_index()
    .pivot(index=['variant_id', 'gene', 'significance'],
           columns='output_type',
           values='raw_score')
    .fillna(0.0)
    .reset_index()
)

# Compute merged score
merged['alphagenome_splicing'] = (
    merged.get('SPLICE_SITES', pd.Series([0.0]*len(merged))).values
    + merged.get('SPLICE_SITE_USAGE', pd.Series([0.0]*len(merged))).values
    + merged.get('SPLICE_JUNCTIONS', pd.Series([0.0]*len(merged))).values / 5.0
)

# Flag risk level
merged['splicing_risk'] = merged['alphagenome_splicing'].apply(
    lambda x: 'HIGH' if x > 1.0 else ('MODERATE' if x > 0.5 else 'LOW')
)

os.makedirs(os.path.dirname(output_path), exist_ok=True)
merged.to_csv(output_path, index=False)

print(f"\n✅ Saved splicing scores to {output_path}")
print(merged[['variant_id', 'gene', 'alphagenome_splicing', 'splicing_risk']])