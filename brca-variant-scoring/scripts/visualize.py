# scripts/visualize.py

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

from alphagenome.data import genome
from alphagenome.models import dna_client

# ── 1. Read inputs ────────────────────────────────────────────────────────────
ranked_path = snakemake.input[0]
vcf_path    = snakemake.input[1]
plots_dir   = snakemake.output[0]
api_key     = snakemake.config['api_key']

os.makedirs(plots_dir, exist_ok=True)

ranked = pd.read_csv(ranked_path)
vcf    = pd.read_csv(vcf_path, sep='\t')

print("Connecting to AlphaGenome...")
dna_model = dna_client.create(api_key)

# ── 2. Discover available ontology terms for RNA_SEQ ─────────────────────────
# Instead of hardcoding a tissue, query the model's own metadata
# and pick the first available RNA_SEQ ontology term
print("Fetching available track metadata...")
metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

rna_meta = metadata[metadata['output_type'] == 'RNA_SEQ']
atac_meta = metadata[metadata['output_type'] == 'ATAC']

# Pick ontology terms that have data for both RNA_SEQ and ATAC
rna_terms  = rna_meta['ontology_curie'].dropna().unique().tolist()
atac_terms = atac_meta['ontology_curie'].dropna().unique().tolist()

# Terms available in both
shared_terms = list(set(rna_terms) & set(atac_terms))

# Prefer breast-related terms if available, otherwise take first 2
breast_terms = [t for t in shared_terms if any(
    kw in metadata[metadata['ontology_curie'] == t]['biosample_name'].str.lower().str.cat()
    for kw in ['breast', 'mammary', 'mcf', 'k562', 'hela']
)]

# Use breast terms if found, else fall back to first 2 shared terms
# If no shared terms found, pass empty list = use all available tracks
use_terms = [
    'UBERON:0008367',
    'CL:0002327',
    'EFO:0001203',
]
print(f"  Using breast tissue ontology terms: {use_terms}")
# ── 3. Summary bar chart ─────────────────────────────────────────────────────
print("\nGenerating summary bar chart...")

fig, ax = plt.subplots(figsize=(14, 6))

colors = ['#2563eb' if g == 'BRCA1' else '#dc2626' for g in ranked['gene']]
ax.barh(
    y=range(len(ranked)),
    width=ranked['top_quantile_score'],
    color=colors, alpha=0.85, edgecolor='white', linewidth=0.5,
)

for i, row in ranked.iterrows():
    rank_idx = row['rank'] - 1
    if row['splicing_risk'] == 'HIGH':
        ax.text(row['top_quantile_score'] + 0.002, rank_idx,
                '⚡ SPLICE', va='center', fontsize=8,
                color='#7c3aed', fontweight='bold')
    elif row['splicing_risk'] == 'MODERATE':
        ax.text(row['top_quantile_score'] + 0.002, rank_idx,
                '~ splice', va='center', fontsize=8, color='#d97706')

def shorten_sig(s):
    s = str(s)
    if 'Pathogenic' in s and 'Uncertain' not in s and 'Likely' not in s:
        return 'Path.'
    elif 'Uncertain' in s:
        return 'VUS'
    elif 'Conflicting' in s:
        return 'Conflict.'
    elif 'Likely pathogenic' in s:
        return 'Likely Path.'
    return s[:10]

labels = [
    f"{row.variant_id}  [{shorten_sig(row.significance)}]"
    for _, row in ranked.iterrows()
]
ax.set_yticks(range(len(ranked)))
ax.set_yticklabels(labels, fontsize=8)
ax.invert_yaxis()

ax.axvline(x=0.9, color='orange', linestyle='--', linewidth=1, alpha=0.7)
ax.text(0.901, len(ranked) - 0.5, 'HIGH threshold (0.9)',
        fontsize=7, color='orange', va='bottom')

ax.set_xlabel('Top Quantile Score (0 = no effect, 1 = maximum effect)', fontsize=10)
ax.set_title(
    'BRCA1/2 Variant Regulatory Impact — AlphaGenome Scores\n'
    'Ranked by quantile score across RNA-seq, ATAC, and histone tracks',
    fontsize=12, fontweight='bold'
)
ax.legend(handles=[
    mpatches.Patch(color='#2563eb', alpha=0.85, label='BRCA1 (chr17)'),
    mpatches.Patch(color='#dc2626', alpha=0.85, label='BRCA2 (chr13)'),
], loc='lower right', fontsize=9)
ax.set_xlim(0, 1.08)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'summary_bar_chart.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  ✅ Saved: summary_bar_chart.png")

# ── 4. REF vs ALT plots for top 3 variants ───────────────────────────────────
print("\nGenerating REF vs ALT track plots for top 3 variants...")

for _, row in ranked.head(3).iterrows():
    variant_id = row['variant_id']
    print(f"\n  Plotting {variant_id}...")

    # Parse variant_id: chr13:32394749:G>A
    parts        = variant_id.split(':')
    chrom        = parts[0]
    pos          = int(parts[1])
    ref_b, alt_b = parts[2].split('>')

    variant  = genome.Variant(
        chromosome=chrom, position=pos,
        reference_bases=ref_b, alternate_bases=alt_b,
        name=variant_id,
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

    # ── Try to get predictions, skip gracefully if empty ────────────────────
    try:
        output = dna_model.predict_variant(
            interval=interval,
            variant=variant,
            organism=dna_client.Organism.HOMO_SAPIENS,
            requested_outputs=[
                dna_client.OutputType.RNA_SEQ,
                dna_client.OutputType.ATAC,
            ],
            ontology_terms=use_terms,
        )

        # Verify we actually got tracks back
        n_rna  = output.reference.rna_seq.values.shape[1]
        n_atac = output.reference.atac.values.shape[1]
        print(f"    RNA_SEQ tracks: {n_rna}  |  ATAC tracks: {n_atac}")

        if n_rna == 0 and n_atac == 0:
            print(f"    ⚠️  No tracks returned — skipping plot for {variant_id}")
            continue

    except Exception as e:
        print(f"    ⚠️  predict_variant failed for {variant_id}: {e}")
        continue

    # ── Window: 10kb around variant ──────────────────────────────────────────
    half_win   = 5000
    win_start  = pos - half_win
    win_end    = pos + half_win
    idx_start  = max(0, win_start - interval.start)
    idx_end    = min(output.reference.rna_seq.values.shape[0],
                     win_end   - interval.start)
    positions  = np.arange(win_start, win_start + (idx_end - idx_start))

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    track_configs = [
        ('RNA-seq', 'rna_seq', '#94a3b8', '#dc2626', 'RNA-seq signal', n_rna),
        ('ATAC',    'atac',    '#94a3b8', '#2563eb', 'ATAC signal',    n_atac),
    ]

    for ax_idx, (label, attr, c_ref, c_alt, ylabel, n_tracks) in enumerate(track_configs):
        ax = axes[ax_idx]

        if n_tracks == 0:
            ax.text(0.5, 0.5, f'No {label} tracks available for selected tissues',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=10, color='grey')
            ax.set_ylabel(ylabel, fontsize=9)
            continue

        # Average across all returned tracks (may be >1 tissue)
        ref_vals = getattr(output.reference, attr).values[idx_start:idx_end, :].mean(axis=1)
        alt_vals = getattr(output.alternate, attr).values[idx_start:idx_end, :].mean(axis=1)

        ax.fill_between(positions, ref_vals, alpha=0.5, color=c_ref, label='REF')
        ax.fill_between(positions, alt_vals, alpha=0.5, color=c_alt, label='ALT')

        # Shade the difference
        ax.fill_between(positions, ref_vals, alt_vals,
                        where=(alt_vals > ref_vals),
                        alpha=0.35, color=c_alt)
        ax.fill_between(positions, ref_vals, alt_vals,
                        where=(alt_vals < ref_vals),
                        alpha=0.35, color='#7c3aed')

        ax.axvline(x=pos, color='black', linestyle='--',
                   linewidth=1.2, alpha=0.8, label='Variant position')
        ax.set_ylabel(ylabel, fontsize=9)
        ax.legend(fontsize=8, loc='upper right')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_title(
        f'{variant_id}  |  {row.gene}  |  {shorten_sig(row.significance)}\n'
        f'REF vs ALT — averaged across {len(use_terms)} tissue(s)\n'
        f'Quantile: {row.top_quantile_score:.4f}  |  '
        f'Splicing: {row.splicing_risk}  |  '
        f'Impact: {row.overall_impact}',
        fontsize=10, fontweight='bold'
    )
    axes[1].set_xlabel(f'Genomic position ({chrom})', fontsize=9)

    plt.tight_layout()
    safe_id  = variant_id.replace(':', '_').replace('>', 'to')
    out_path = os.path.join(plots_dir, f'REF_vs_ALT_{safe_id}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ Saved: {out_path}")

print(f"\n✅ All plots saved to {plots_dir}")