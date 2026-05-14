# scripts/rank_and_report.py

import pandas as pd
import os

raw_path      = snakemake.input[0]   # results/raw_scores.csv
splicing_path = snakemake.input[1]   # results/splicing_scores.csv
output_path   = snakemake.output[0]  # results/ranked_variants.csv

# ── Load both score tables ───────────────────────────────────────────────────
raw_df      = pd.read_csv(raw_path)
splicing_df = pd.read_csv(splicing_path)

# ── Summarise raw scores: best hit per variant ───────────────────────────────
# For each variant, find the track where it scored highest (quantile_score)
# quantile_score = how extreme this variant is compared to all variants
best_hits = (
    raw_df.sort_values('quantile_score', ascending=False)
    .groupby('variant_id')
    .first()
    .reset_index()
    [['variant_id', 'gene', 'significance', 'output_type',
      'biosample_name', 'quantile_score', 'raw_score']]
    .rename(columns={
        'output_type':    'top_output_type',
        'biosample_name': 'top_tissue',
        'quantile_score': 'top_quantile_score',
        'raw_score':      'top_raw_score',
    })
)

# ── Join splicing risk ───────────────────────────────────────────────────────
splicing_summary = splicing_df[
    ['variant_id', 'alphagenome_splicing', 'splicing_risk']
]

report = best_hits.merge(splicing_summary, on='variant_id', how='left')

# ── Rank by top quantile score (highest = most impactful) ───────────────────
report = report.sort_values('top_quantile_score', ascending=False)
report.insert(0, 'rank', range(1, len(report) + 1))

# ── Add an overall impact flag ───────────────────────────────────────────────
def flag_impact(row):
    if row['top_quantile_score'] > 0.9 or row['splicing_risk'] == 'HIGH':
        return 'HIGH'
    elif row['top_quantile_score'] > 0.7 or row['splicing_risk'] == 'MODERATE':
        return 'MODERATE'
    else:
        return 'LOW'

report['overall_impact'] = report.apply(flag_impact, axis=1)

os.makedirs(os.path.dirname(output_path), exist_ok=True)
report.to_csv(output_path, index=False)

print(f"\n✅ Final ranked report saved to {output_path}")
print(f"\n{'='*60}")
print("TOP 5 VARIANTS BY REGULATORY IMPACT:")
print('='*60)
print(report.head(5)[
    ['rank', 'variant_id', 'gene', 'significance',
     'top_output_type', 'top_quantile_score', 'splicing_risk', 'overall_impact']
].to_string(index=False))