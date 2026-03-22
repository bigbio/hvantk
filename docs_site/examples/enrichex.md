# EnrichEx Examples

This directory contains example scripts and results for EnrichEx gene set enrichment analysis and burden testing workflows.

## Quick Start (CLI)

EnrichEx provides two main analysis workflows via the `hvantk enrichex` command:

### Overlap Enrichment

Test if a gene list is enriched in gene sets using Fisher's exact test:

```bash
# Basic usage
hvantk enrichex overlap \
  -g my_genes.txt \
  -s synthetic_gene_sets.json \
  -o results/overlap_results.tsv

# With report generation (plots + HTML)
hvantk enrichex overlap \
  -g my_genes.txt \
  -s synthetic_gene_sets.json \
  -o results/overlap_results.tsv \
  --generate-report
```

### Burden Testing

Test if cases have excess rare variants in gene set genes using Hail regression:

```bash
# Basic usage
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s synthetic_gene_sets.json \
  -o results/burden_results.tsv

# With covariates and report generation
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s synthetic_gene_sets.json \
  --covariates PC1,PC2,PC3,age,sex \
  --max-af 0.001 \
  --min-cadd 25 \
  -o results/burden_results.tsv \
  --generate-report
```

## Directory Structure

```
examples/enrichex/
├── README.md                         # This file
├── synthetic_gene_sets.json          # Example gene sets for testing
├── overlap_enrichment_example.py     # Python API example for overlap
├── burden_analysis_example.py        # Python API example for burden
├── create_gene_sets_example.py       # How to create gene set collections
├── synthetic_burden_demo.py          # Synthetic cohort burden testing demo
├── chd_burden_run.py                 # CHD burden analysis workflow
└── results/                          # Example output from test runs
    ├── overlap_results.tsv           # Overlap enrichment results
    ├── enrichex_overlap.png          # Overlap enrichment dot plot
    ├── enrichex_overlap_report.html  # HTML report for overlap
    ├── burden_results.tsv            # Burden testing results
    ├── enrichex_burden.png           # Burden forest plot
    └── enrichex_burden_report.html   # HTML report for burden
```

## Example Files

### synthetic_gene_sets.json

Example gene sets for brain cell types (Microglia, Astrocytes, Oligodendrocytes, etc.) ready for use with EnrichEx commands.

### Python API Examples

| Script | Description |
|--------|-------------|
| [`overlap_enrichment_example.py`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/overlap_enrichment_example.py) | Run overlap enrichment using Python API |
| [`burden_analysis_example.py`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/burden_analysis_example.py) | Run burden testing using Python API with synthetic data |
| [`create_gene_sets_example.py`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/create_gene_sets_example.py) | Create gene set collections from various formats |
| [`synthetic_burden_demo.py`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/synthetic_burden_demo.py) | Synthetic cohort burden testing demo |
| [`chd_burden_run.py`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/chd_burden_run.py) | CHD burden analysis workflow |

Run examples:
```bash
cd examples/enrichex
python overlap_enrichment_example.py
python burden_analysis_example.py
python create_gene_sets_example.py
```

### results/

Pre-generated results from successful test runs showing expected output format:
- **TSV files**: Tab-separated results with statistics
- **PNG files**: Publication-ready plots
- **HTML files**: Self-contained reports with embedded plots

## Input Requirements

### Gene List (for overlap)

Plain text file with one gene symbol per line:
```
APOE
TREM2
CD33
MS4A6A
```

### Gene Sets (JSON)

JSON file with gene set definitions:
```json
{
  "gene_sets": {
    "Microglia": {
      "name": "Microglia",
      "genes": ["APOE", "TREM2", "CD33", "MS4A6A"]
    }
  },
  "background_genes": ["APOE", "TREM2", "CD33", ...]
}
```

### Cohort MatrixTable (for burden)

Hail MatrixTable with:
- Row annotations: `SYMBOL` (gene name), `gnomad_af`, `cadd_phred`
- Entry fields: `GT` (genotype)
- Column key: sample ID

### Phenotype Table (for burden)

Hail Table or TSV with:
- Sample IDs matching MatrixTable columns
- Phenotype field (binary or continuous)
- Optional covariates (PCs, age, sex)

See [EnrichEx Docs](../tools/enrichex.md) for full CLI options and output format details.

## Gene Set Sources

**Cell-type markers:**
- Lake et al. 2018: https://www.nature.com/articles/nbt.4038
- PanglaoDB: https://panglaodb.se/
- CellMarker: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

**Pathway databases:**
- MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/
- Gene Ontology: http://geneontology.org/
- KEGG: https://www.genome.jp/kegg/pathway.html

**GWAS genes:**
- GWAS Catalog: https://www.ebi.ac.uk/gwas/

## Documentation

- **Full documentation**: [EnrichEx Docs](../tools/enrichex.md)
- **CLI help**: `hvantk enrichex --help`
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
