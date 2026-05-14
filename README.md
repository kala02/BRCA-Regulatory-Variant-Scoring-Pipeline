# BRCA Regulatory Variant Scoring Pipeline

A Snakemake pipeline that uses [Google DeepMind's AlphaGenome](https://deepmind.google/science/alphagenome/) to predict the regulatory and splicing impact of non-coding variants near **BRCA1** and **BRCA2** — the two most clinically significant hereditary breast and ovarian cancer genes.

---

## Background

BRCA1 and BRCA2 are tumour suppressor genes that maintain genomic stability through DNA damage repair. Pathogenic variants in these genes confer a lifetime breast cancer risk of 50–70% and ovarian cancer risk of 15–45%. However, a large fraction of variants identified in clinical testing are classified as **Variants of Uncertain Significance (VUS)** — variants for which there is insufficient evidence to determine pathogenicity.

Reclassification of VUS requires functional evidence. Traditional computational tools (e.g. SIFT, PolyPhen) focus on coding variants and cannot assess the regulatory consequences of non-coding changes. AlphaGenome addresses this gap by predicting genome-wide molecular phenotypes — gene expression, chromatin accessibility, histone modifications, and splicing — directly from DNA sequence at single base-pair resolution.

This pipeline applies AlphaGenome to systematically score BRCA1/2 variants for regulatory and splicing disruption, providing a ranked priority list to guide functional validation efforts.

---

## Key Finding

**`chr17:43076605 G>A` (BRCA1) — currently classified as Uncertain Significance in ClinVar — was flagged HIGH by two independent evidence streams:**

| Evidence | Score | Threshold | Classification |
|---|---|---|---|
| Regulatory impact (quantile) | 0.99998 | > 0.9 | HIGH |
| Merged splicing score | 1.572 | > 1.0 | HIGH |
| SPLICE_JUNCTIONS component | 2.975 | — | Highest in dataset |
| SPLICE_SITES component | 0.520 | — | Elevated |
| SPLICE_SITE_USAGE component | 0.458 | — | Elevated |

All three independent splicing models fired simultaneously for this variant, consistent with predicted exon skipping or cryptic splice site activation. This represents strong computational evidence for reclassification as **Likely Pathogenic**.

---

## Pipeline Overview

```
data/brca_variants.vcf
        │
        ├──▶ Rule 1: score_variants     ──▶ results/raw_scores.csv
        │    (AlphaGenome batch scoring)     (430,372 rows — variant × track × score)
        │
        ├──▶ Rule 2: compute_splicing   ──▶ results/splicing_scores.csv
        │    (3 splicing scorers +           (20 rows — merged score per variant)
        │     merged splicing formula)
        │
        └──▶ Rule 3: rank_and_report    ──▶ results/ranked_variants.csv
             (join + rank + flag)            (20 rows — final ranked report)
                    │
                    └──▶ Rule 4: visualize  ──▶ results/plots/
                         (REF vs ALT plots       summary_bar_chart.png
                          + summary chart)        REF_vs_ALT_*.png (top 3)
```

Each rule is a separate Python script with a single responsibility. If ranking logic changes, only Rule 3 re-runs — the expensive API calls in Rules 1 and 2 are skipped automatically by Snakemake.

---

## Results

### Summary bar chart — all 20 variants ranked


All 20 variants exceeded a quantile score of 0.98, reflecting the clinical pre-selection of the variant set. The `⚡ SPLICE` marker identifies `chr17:43076605` as the only variant with both top regulatory and splicing impact simultaneously.

---

### REF vs ALT plots — top 3 variants

**Rank 1 — `chr13:32394749 G>A` | BRCA2 | Pathogenic**

Known pathogenic variant. Scored at the 99.998th quantile — AlphaGenome independently confirms clinical classification, validating the pipeline's accuracy. Strongest effect observed in non-breast tissue (SJCRH30), consistent with BRCA2's ubiquitous role in DNA repair. Splicing score LOW — mechanism is primarily regulatory.

---

**Rank 2 — `chr17:43076605 G>A` | BRCA1 | Uncertain Significance ⭐**


The headline finding. Visual breast-tissue tracks show subtle REF vs ALT overlap — the effect is not primarily in breast chromatin accessibility. The splicing scores tell the real story: `SPLICE_JUNCTIONS = 2.97`, the highest junction disruption score in the dataset, combined with elevated SPLICE_SITES (0.52) and SPLICE_SITE_USAGE (0.46). The convergence of all three independent splicing models constitutes strong computational evidence for aberrant splicing — likely exon skipping — that would produce a truncated, non-functional BRCA1 protein.

> **Clinical implication:** This VUS should be prioritised for functional validation via minigene splicing assay or patient-derived cell lines. Confirmation would support reclassification as Likely Pathogenic under ACMG/AMP guidelines (PS3 evidence criterion).

---

**Rank 3 — `chr17:43092314 C>T` | BRCA1 | Conflicting classifications**


ClinVar has conflicting submissions for this variant. AlphaGenome scores it at the 99.93rd quantile — supporting the pathogenic submissions. The RNA-seq plot clearly shows the BRCA1 gene body as an active expression block (~5kb wide), with the variant positioned at its left boundary near a splice junction. The strongest predicted effect is in HUES64 (embryonic stem cells), consistent with BRCA1's critical role in early development and proliferating cell types.

---

### Full ranked report — all 20 variants

| Rank | Variant | Gene | ClinVar | Quantile | Splicing | Impact |
|---|---|---|---|---|---|---|
| 1 | chr13:32394749:G>A | BRCA2 | Pathogenic | 0.99998 | LOW | HIGH |
| 2 | chr17:43076605:G>A | BRCA1 | **VUS** | 0.99998 | **HIGH** | HIGH |
| 3 | chr17:43092314:C>T | BRCA1 | Conflicting | 0.9993 | LOW | HIGH |
| 4 | chr17:43104949:G>A | BRCA1 | Pathogenic | 0.9992 | MODERATE | HIGH |
| 5 | chr17:43051110:C>G | BRCA1 | Conflicting | 0.9992 | LOW | HIGH |
| 6–20 | … | … | … | > 0.987 | LOW | HIGH |


---

## Methods

### Variant selection

Non-coding single nucleotide variants near BRCA1 and BRCA2 were sourced from ClinVar (GRCh38 assembly). Variants were filtered to Pathogenic and Uncertain Significance classifications with validated VCF-format alleles. 10 variants per gene were sampled randomly (seed = 42) for a balanced 20-variant test set.

### Regulatory scoring

Each variant was scored using AlphaGenome's `score_variant()` function with three output types: RNA_SEQ, ATAC, and CHIP_HISTONE. The model predicts reference and alternate allele sequences across a 1Mb window centered on each variant and computes the difference across 5,563 human tissue tracks. Variants were ranked by `quantile_score` — the genome-wide percentile rank of the predicted effect size.

### Splicing scoring

Three dedicated splicing scorers were applied independently: SPLICE_SITES (changes in donor/acceptor site probabilities), SPLICE_SITE_USAGE (changes in relative splice site utilisation), and SPLICE_JUNCTIONS (changes in exon-exon junction expression). A merged splicing score was computed as:

```
alphagenome_splicing = SPLICE_SITES_max
                     + SPLICE_SITE_USAGE_max
                     + SPLICE_JUNCTIONS_max / 5.0
```

Variants with a merged score > 1.0 were classified HIGH risk, > 0.5 MODERATE, and < 0.5 LOW, following AlphaGenome's recommended thresholds.

### Visualisation

REF vs ALT track plots were generated for the top 3 variants using three breast-relevant ontology terms: `UBERON:0008367` (breast epithelium), `CL:0002327` (mammary epithelial cell), and `EFO:0001203` (MCF-7 breast cancer cell line). Tracks were averaged across the three tissues per output type.

---

## Limitations

**Variant selection bias.** ClinVar over-represents variants that were already clinically suspected. A future extension would apply this pipeline to uncharacterised rare variants from gnomAD to identify novel candidates.

**Computational predictions require experimental validation.** Splicing predictions for the top VUS should be confirmed with a minigene splicing assay or RNA sequencing of patient-derived cell lines before influencing clinical decisions.

**Tissue-specificity of scores.** The highest quantile scores for most variants were observed in non-breast tissues, consistent with BRCA1/2's ubiquitous role in DNA repair. Breast-tissue-specific visualisation was performed separately using validated mammary ontology terms.

**Sample size.** 20 variants were used as a proof-of-concept. The pipeline is designed to scale to full ClinVar VUS datasets.

---

## Requirements

```
Python >= 3.11
snakemake
alphagenome
pandas
matplotlib
numpy
```

Install dependencies:

```bash
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
pip install snakemake alphagenome pandas matplotlib numpy
```

---

## Setup

**1. Clone the repository**

```bash
git clone https://github.com/BRCA Regulatory Variant Scoring Pipeline/brca-variant-scoring.git
cd brca-variant-scoring
```

**2. Get an AlphaGenome API key**

Register at [deepmind.google/science/alphagenome](https://deepmind.google/science/alphagenome). The API is free for non-commercial research use.

**3. Add your API key**

```bash
echo 'api_key: "ag-your-key-here"' > config.yaml
echo "config.yaml" >> .gitignore   # never commit your key
```

**4. Prepare input variants**

Run the data preparation notebook to download ClinVar and generate `data/brca_variants.vcf`:

```bash
jupyter notebook notebook/data_preparation.ipynb
```

---

## Running the pipeline

```bash
# Dry run — check the execution plan without running anything
snakemake --dry-run

# Full run — executes all 4 rules sequentially
snakemake -j1
```

Expected runtime: ~25–30 minutes (20 API calls × 2 scoring scripts + visualisation).

Output files:

```
results/
├── raw_scores.csv           # 430,372 rows — full track-level scores
├── splicing_scores.csv      # 20 rows — merged splicing score per variant
├── ranked_variants.csv      # 20 rows — final ranked report
└── plots/
    ├── summary_bar_chart.png
    ├── REF_vs_ALT_chr13_32394749_GtoA.png
    ├── REF_vs_ALT_chr17_43076605_GtoA.png
    └── REF_vs_ALT_chr17_43092314_CtoT.png
```

---

## Project structure

```
brca-variant-scoring/
├── Snakefile                      # pipeline definition
├── config.yaml                    # API key (gitignored)
├── .gitignore
├── data/
│   └── brca_variants.vcf          # 20 input variants
├── scripts/
│   ├── score_variants.py          # Rule 1 — regulatory scoring
│   ├── compute_splicing.py        # Rule 2 — splicing scoring
│   ├── rank_and_report.py         # Rule 3 — ranking + report
│   └── visualize.py               # Rule 4 — REF vs ALT plots
├── results/                       # generated by pipeline
│   ├── raw_scores.csv
│   ├── splicing_scores.csv
│   ├── ranked_variants.csv
│   └── plots/
└── notebook/
    └── data_preparation.ipynb     # ClinVar download + VCF preparation
```

---

## Tools and data sources

| Tool / Resource | Purpose |
|---|---|
| [AlphaGenome](https://deepmind.google/science/alphagenome/) | Variant effect prediction |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical variant annotations |
| [Snakemake](https://snakemake.readthedocs.io/) | Pipeline management |
| GRCh38 | Human reference genome assembly |

---

## Acknowledgements

AlphaGenome is developed by Google DeepMind. This project uses the AlphaGenome API under non-commercial research terms. Variant data sourced from NCBI ClinVar.
