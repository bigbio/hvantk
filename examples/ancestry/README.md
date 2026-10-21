# Ancestry Inference Examples

Example scripts for the ancestry inference pipeline (PCA + Random Forest classification).

For full documentation, see the [Ancestry docs](https://bigbio.github.io/hvantk/tools/ancestry/) and [Ancestry examples guide](https://bigbio.github.io/hvantk/examples/ancestry/).

## Contents

| Script | Description |
|--------|-------------|
| `run_ancestry_example.py` | End-to-end example with synthetic data (recommended starting point) |
| `basic_inference.py` | Basic ancestry inference workflow |
| `custom_reference.py` | Using a custom reference panel with population mapping |
| `cli_examples.sh` | CLI usage examples |

## Quick start

```bash
# Run end-to-end example (generates synthetic data automatically)
python examples/ancestry/run_ancestry_example.py

# With custom output directory
python examples/ancestry/run_ancestry_example.py --output-dir ./ancestry_results
```

## Expected outputs

`run_ancestry_example.py` writes (under `--output-dir`, default `./examples/ancestry/results`):

```text
results/
├── data/
│   ├── synthetic_query.mt/
│   └── synthetic_reference.mt/
├── predictions.ht/
├── predictions.tsv
├── pc_scores.ht/
├── rf_model.pkl
├── pca_loadings.ht/
├── pipeline_stats.json
├── ancestry_report.html
├── plots/
│   ├── pca_scatter_pc1_pc2.png
│   ├── pca_scatter_pc1_pc3.png
│   ├── pca_panel.png
│   ├── ancestry_proportions.png
│   ├── probability_distribution.png
│   ├── variance_explained.png
│   └── confusion_matrix.png
└── checkpoints/
```

`basic_inference.py` and `custom_reference.py` write a smaller subset of this
(`predictions.ht`, `predictions.tsv`, `ancestry_report.html`, plus a handful of
individual plot PNGs) directly under their own `--output-dir`, with no `data/`,
`pc_scores.ht`, or `checkpoints/`.

## Requirements

- Hail installed and configured
- Reference panel MatrixTable (e.g., 1000 Genomes) for real analysis
- Query cohort MatrixTable keyed by `(locus, alleles)`
