# hvantk Examples

This directory contains end-to-end workflow examples demonstrating how to use hvantk for genomic variant analysis. Each subdirectory represents a complete analysis pipeline with scripts, results, and documentation.

## Quick Links

- [HGC (Joint Genotyping)](#hgc-joint-genotyping)
- [PSROC (Prediction Score ROC Analysis)](#psroc-prediction-score-roc-analysis)
- [EnrichEx (Gene Set Enrichment)](#enrichex-gene-set-enrichment)
- [PTM (Post-Translational Modification)](#ptm-post-translational-modification)
- [Ancestry Inference](#ancestry-inference)
- [ClinVar Data Streaming](#clinvar-data-streaming)
- [ClinGen Streaming](#clingen-streaming)
- [1000 Genomes Reference Build](#1000-genomes-reference-build)
- [Recipe Templates](#recipe-templates)

## Workflow Examples

### HGC (Joint Genotyping)

**Directory:** [`hgc/`](https://github.com/bigbio/hvantk/tree/main/examples/hgc/)

Complete pipeline for joint genotyping of GVCF cohorts with quality control and benchmarking.

**Key scripts:**
- `qc/hgc_qc_example.py` - QC workflow with visualization
- `cpu_scaling/benchmark.py` - CPU scalability testing
- `scalability/benchmark.py` - Sample size scalability

**Quick start:**
```bash
python examples/hgc/qc/hgc_qc_example.py
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

### PTM (Post-Translational Modification)

Map UniProt PTM sites to genomic coordinates, cross-reference with ClinVar and gnomAD variants, and analyze PTM-variant enrichment.

**Quick start:**
```bash
hvantk ptm build --output-dir data/ptm/ --output-ht data/ptm/ptm_sites.ht
hvantk ptm landscape --clinvar-ht clinvar.ht --ptm-ht ptm_sites.ht -o results/landscape/ --save-plots
hvantk ptm report -o report.html --landscape-json results/landscape/landscape_summary.json
```

**Outputs:** PTM sites table (HT), landscape/population JSON summaries, plots (PNG), HTML report

**Documentation:** [PTM Examples](ptm.md) | [PTM Docs](../tools/ptm.md)

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

### Ancestry Inference

**Directory:** [`ancestry/`](https://github.com/bigbio/hvantk/tree/main/examples/ancestry/)

Predict genetic ancestry for samples using PCA and Random Forest classification against a labeled reference panel.

**Key scripts:**
- `basic_inference.py` - Basic ancestry inference workflow
- `custom_reference.py` - Using a custom reference panel
- `run_ancestry_example.py` - End-to-end example with reporting

**Quick start:**
```bash
python examples/ancestry/basic_inference.py
```

**Outputs:** Ancestry predictions (HT/TSV), PCA plots (PNG), HTML reports

**Documentation:** [Ancestry Examples](ancestry.md) | [Ancestry Docs](../tools/ancestry.md)

---

### ClinGen Streaming

**Directory:** [`clingen/`](https://github.com/bigbio/hvantk/tree/main/examples/clingen/)

Examples for querying ClinGen gene-disease validity data and categorizing by ontology.

**Key scripts:**
- `run_with_real_data.py` - ClinGen streamer with real data
- `run_ontology_categorization.py` - Disease ontology categorization

**Quick start:**
```bash
python examples/clingen/run_with_real_data.py
```

**Documentation:** [ClinGen Examples](clingen.md)

---

### 1000 Genomes Reference Build

**Directory:** [`1k_genome/`](https://github.com/bigbio/hvantk/tree/main/examples/1k_genome/)

Scripts for building a Hail MatrixTable from 1000 Genomes NYGC high-coverage data for use as a reference panel.

**Key scripts:**
- `build_1kg_nygc.py` - Build reference MatrixTable from 1000 Genomes
- `build_1kg_nygc_cli.py` - CLI wrapper for the build script

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
