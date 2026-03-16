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
| `overlap_enrichment_example.py` | Run overlap enrichment using Python API |
| `burden_analysis_example.py` | Run burden testing using Python API with synthetic data |
| `create_gene_sets_example.py` | Create gene set collections from various formats |

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

## CLI Options

### Overlap Command

```bash
hvantk enrichex overlap --help
```

Key options:
- `-g, --gene-list`: Query gene list file (required)
- `-s, --gene-sets`: Gene sets JSON file (required)
- `-o, --output`: Output TSV path (required)
- `--correction`: Multiple testing correction (benjamini-hochberg, bonferroni, none)
- `--alpha`: Significance threshold (default: 0.05)
- `--generate-report`: Generate plots and HTML report

### Burden Command

```bash
hvantk enrichex burden --help
```

Key options:
- `-m, --cohort-mt`: Cohort MatrixTable (required)
- `-p, --phenotypes`: Phenotype file (required)
- `-s, --gene-sets`: Gene sets JSON file (required)
- `-o, --output`: Output TSV path (required)
- `--phenotype-field`: Phenotype column name (default: is_case)
- `--phenotype-type`: binary or continuous
- `--covariates`: Comma-separated covariate names
- `--max-af`: Max allele frequency filter (default: 0.01)
- `--min-cadd`: Min CADD score filter (default: 20)
- `--genotype-aggregation`: hets, homs, chets, homs_chets
- `--generate-report`: Generate plots and HTML report

## Output Format

### Overlap Results (TSV)

| Column | Description |
|--------|-------------|
| gene_set_name | Name of the gene set |
| n_query | Number of query genes |
| n_gene_set | Number of genes in gene set |
| n_overlap | Number of overlapping genes |
| n_background | Size of background |
| p_value | Fisher's exact test p-value |
| odds_ratio | Effect size |
| ci_lower, ci_upper | 95% confidence interval |
| overlap_genes | Comma-separated overlapping genes |
| p_adjusted | Multiple testing corrected p-value |
| significant | Boolean significance flag |

### Burden Results (TSV)

| Column | Description |
|--------|-------------|
| gene_set_name | Name of the gene set |
| beta | Regression coefficient |
| standard_error | SE of beta |
| z_stat | Z-statistic |
| p_value | Regression p-value |
| odds_ratio | Exponentiated beta (binary only) |
| ci_lower, ci_upper | 95% confidence interval |
| p_adjusted | Multiple testing corrected p-value |
| significant | Boolean significance flag |

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

- **Full documentation**: [docs/tools/enrichex.md](../../docs_site/tools/enrichex.md)
- **CLI help**: `hvantk enrichex --help`
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
