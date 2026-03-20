# hvantk Examples

This directory contains end-to-end workflow examples demonstrating how to use hvantk for genomic variant analysis. Each subdirectory represents a complete analysis pipeline with scripts, results, and documentation.

## Quick Links

- [HGC (Joint Genotyping)](#hgc-joint-genotyping)
- [PSROC (Prediction Score ROC Analysis)](#psroc-prediction-score-roc-analysis)
- [EnrichEx (Gene Set Enrichment)](#enrichex-gene-set-enrichment)
- [ClinVar Data Streaming](#clinvar-data-streaming)
- [Recipe Templates](#recipe-templates)

## Workflow Examples

### HGC (Joint Genotyping)

**Directory:** [`hgc/`](https://github.com/bigbio/hvantk/tree/main/examples/hgc/)

Complete pipeline for joint genotyping of GVCF cohorts with quality control and benchmarking.

**Key scripts:**
- `hgc_qc_example.py` - QC workflow with visualization
- `hgc_cpu_scaling_benchmark.py` - CPU scalability testing
- `hgc_scalability_benchmark.py` - Sample size scalability

**Quick start:**
```bash
python examples/hgc/hgc_qc_example.py
```

**Outputs:** QC reports (HTML), dashboards (PNG), metrics (JSON)

**Documentation:** [HGC Examples](hgc.md) | [HGC Docs](../tools/hgc.md)

---

### PSROC (Prediction Score ROC Analysis)

**Directory:** [`psroc/`](https://github.com/bigbio/hvantk/tree/main/examples/psroc/)

ROC curve analysis for evaluating variant pathogenicity prediction scores using ClinVar labels.

**Key scripts:**
- `run_psroc_example.py` - Complete PSROC workflow with synthetic data

**Quick start:**
```bash
python examples/psroc/run_psroc_example.py
```

**Outputs:** ROC curves (PNG), AUC metrics (JSON), annotated variants (TSV)

**Documentation:** [PSROC Examples](psroc.md) | [PSROC Docs](../tools/psroc.md)

---

### EnrichEx (Gene Set Enrichment)

**Directory:** [`enrichex/`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/)

Gene set enrichment analysis using overlap testing (Fisher's exact) and case-control burden testing (Hail regression).

**Key scripts:**
- `overlap_enrichment_example.py` - Overlap enrichment using Python API
- `burden_analysis_example.py` - Burden testing with synthetic data
- `create_gene_sets_example.py` - Create gene set collections

**Quick start:**
```bash
# CLI usage (recommended)
hvantk enrichex overlap \
  -g my_genes.txt \
  -s examples/enrichex/synthetic_gene_sets.json \
  -o results.tsv \
  --generate-report

# Python API
python examples/enrichex/overlap_enrichment_example.py
```

**Outputs:** Results (TSV), plots (PNG), HTML reports

**Documentation:** [EnrichEx Examples](enrichex.md) | [EnrichEx Docs](../tools/enrichex.md)

---

### ClinVar Data Streaming

**Directory:** [`clinvar/`](https://github.com/bigbio/hvantk/tree/main/examples/clinvar/)

Examples for filtering and processing ClinVar variant annotations.

**Key scripts:**
- `clinvar_streamer_example.py` - ClinVar filtering and export

**Quick start:**
```bash
python examples/clinvar/clinvar_streamer_example.py
```

**Outputs:** Filtered variant tables (TSV/VCF)

**Documentation:** [ClinVar Examples](clinvar.md)

---

## Recipe Templates

**Directory:** [`recipes/`](https://github.com/bigbio/hvantk/tree/main/examples/recipes/)

Ready-to-use recipe templates for batch processing annotation tables and expression matrices.

### Table Recipes

**`tables.example.json`** - Batch-create annotation tables (Hail Tables)

```bash
hvantk mktable-batch --recipe recipes/tables.example.json
```

Demonstrates:
- ClinVar annotation table creation
- INSIDER interactome table creation
- Reference genome specification
- Multiple output formats

### MatrixTable Recipes

**`matrices.example.json`** - Batch-create expression matrices

```bash
hvantk mkmatrix-batch --recipe recipes/matrices.example.json
```

Demonstrates:
- UCSC Cell Browser data conversion
- Expression Atlas data processing
- Multiple matrix creation in one run

**`cptac.example.json`** - CPTAC protein expression data

```bash
hvantk mkmatrix-batch --recipe recipes/cptac.example.json
```

Demonstrates:
- Protein expression matrix creation
- Categorical metadata handling
- Custom column mapping

---

## Running Examples

### Prerequisites

```bash
# Install dependencies
poetry install

# Activate environment
eval "$(poetry env activate)"
```

### Basic Usage

```bash
# Run a workflow example
python examples/<workflow>/script.py

# Run with custom parameters (if supported)
python examples/hgc/hgc_qc_example.py --input my_data.mt --output my_results/
```

### Using Recipes

```bash
# JSON recipe (built-in support)
hvantk mktable-batch --recipe path/to/recipe.json

# YAML recipe (requires PyYAML)
hvantk mkmatrix-batch --recipe path/to/recipe.yaml
```

---

## Example Structure

Each workflow directory follows this structure:

```
<workflow>/
├── README.md              # Workflow-specific documentation
├── *.py                   # Main example scripts
├── scripts/               # Supporting scripts (optional)
│   ├── plot_*.py         # Visualization scripts
│   └── run_*.sh          # Shell automation scripts
└── results/              # Example outputs
    ├── *.json            # Metrics and metadata
    ├── *.tsv             # Tabular results
    └── plots/            # Visualizations (PNG/PDF)
```



See the [Quick Start Guide](../getting-started/quickstart.md) for common workflows, and [Contributing](../contributing.md) for test data information.
