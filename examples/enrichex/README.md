# EnrichEx Examples

Example scripts for gene set overlap enrichment and burden testing.

For full documentation, see the [EnrichEx docs](https://bigbio.github.io/hvantk/tools/enrichex/) and [EnrichEx examples guide](https://bigbio.github.io/hvantk/examples/enrichex/).

## Directory Structure

```
examples/enrichex/
├── synthetic_gene_sets.json          # Example gene sets for testing
├── overlap_enrichment_example.py     # Python API: overlap enrichment
├── burden_analysis_example.py        # Python API: burden testing with synthetic data
├── create_gene_sets_example.py       # Python API: create gene set collections
├── synthetic_burden_demo.py          # Synthetic cohort burden testing demo
├── chd_burden_run.py                 # CHD burden analysis workflow
└── results/                          # Pre-generated example outputs
    ├── overlap_results.tsv
    ├── enrichex_overlap.png
    ├── enrichex_overlap_report.html
    ├── burden_results.tsv
    ├── enrichex_burden.png
    └── enrichex_burden_report.html
```

## Quick Start

### CLI

```bash
# Overlap enrichment
hvantk enrichex overlap \
  -g my_genes.txt \
  -s synthetic_gene_sets.json \
  -o results/overlap_results.tsv \
  --generate-report

# Burden testing
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s synthetic_gene_sets.json \
  -o results/burden_results.tsv \
  --generate-report
```

### Python API

```bash
python overlap_enrichment_example.py
python burden_analysis_example.py
python create_gene_sets_example.py
```

## Input Requirements

- **Gene list** (overlap): plain text, one gene symbol per line
- **Gene sets**: JSON with `gene_sets` and optional `background_genes`
- **Cohort MatrixTable** (burden): with `SYMBOL`, `gnomad_af`, `cadd_phred` row annotations
- **Phenotype table** (burden): sample IDs + phenotype field + optional covariates

## Expected Outputs

- TSV results with statistics (p-values, odds ratios, confidence intervals)
- PNG plots (dot plot for overlap, forest plot for burden)
- Self-contained HTML reports
