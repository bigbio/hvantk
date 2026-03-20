# Ancestry Inference Examples

This directory contains example scripts demonstrating the ancestry inference pipeline.

## Quick Start

Run the end-to-end example with synthetic data:

```bash
# Activate environment
eval $(poetry env activate)

# Run end-to-end example (generates synthetic data automatically)
python examples/ancestry/run_ancestry_example.py
```

This generates outputs in `examples/ancestry/results/` including:
- `ancestry_report.html` - Comprehensive HTML report
- `predictions.tsv` - Ancestry predictions for all samples
- `plots/` - Visualizations (PCA, confusion matrix, etc.)

## Examples

### 1. End-to-End Example (`run_ancestry_example.py`) ⭐

**Recommended starting point.** Complete runnable workflow:
- Generates synthetic data (5 populations, 300 samples)
- Runs full ancestry inference pipeline
- Validates predictions against true labels
- Generates all visualizations and HTML report
- Demonstrates validation of inference accuracy

```bash
python examples/ancestry/run_ancestry_example.py --output-dir ./ancestry_results
```

### 2. Basic Inference (`basic_inference.py`)

Demonstrates basic ancestry inference workflow:
- Loading query and reference MatrixTables
- Running the inference pipeline
- Examining predictions
- Generating HTML report

### 3. Custom Reference Panel (`custom_reference.py`)

Shows how to use a custom reference panel:
- Preparing a custom reference MatrixTable
- Population label mapping (e.g., CEU → EUR)
- Configuring pipeline parameters
- Custom visualization colors

### 4. CLI Examples (`cli_examples.sh`)

Bash script with CLI usage examples:
- Basic usage
- Custom parameters
- Report generation
- Checkpointing
- Model and loadings export

## Results Directory

The `results/` directory contains example outputs from running the pipeline:

```
results/
├── ancestry_report.html    # Interactive HTML report
├── predictions.tsv         # Ancestry predictions (sample_id, predicted, probability)
├── pipeline_stats.json     # Pipeline execution statistics
└── plots/
    ├── pca_scatter_pc1_pc2.png
    ├── pca_scatter_pc1_pc3.png
    ├── pca_panel.png
    ├── ancestry_proportions.png
    ├── probability_distribution.png
    ├── variance_explained.png
    └── confusion_matrix.png
```

## Prerequisites

- hvantk installed with Hail support
- Reference panel MatrixTable (e.g., 1000 Genomes) for real analysis
- Query cohort MatrixTable

## Data Requirements

The examples expect:
- Query MT: A Hail MatrixTable with genotype calls keyed by `(locus, alleles)`
- Reference MT: A Hail MatrixTable with ancestry labels in a column annotation

For testing, you can generate synthetic data using Hail's `balding_nichols_model()`:

```python
import hail as hl

# Generate combined dataset (reference + query from same simulation)
combined_mt = hl.balding_nichols_model(
    n_populations=5,
    n_samples=300,  # 200 reference + 100 query
    n_variants=10000,
    pop_dist=[0.2, 0.2, 0.2, 0.2, 0.2],
    fst=[0.12, 0.15, 0.10, 0.08, 0.11],
)

# Add population labels
pop_labels = hl.literal(["EUR", "AFR", "EAS", "SAS", "AMR"])
combined_mt = combined_mt.annotate_cols(
    ancestry=pop_labels[combined_mt.pop],
)

# Split into reference (first 200) and query (last 100)
ref_mt = combined_mt.filter_cols(combined_mt.sample_idx < 200)
query_mt = combined_mt.filter_cols(combined_mt.sample_idx >= 200)
```

## Documentation

For full documentation, see [docs/tools/ancestry.md](../../docs_site/tools/ancestry.md).
