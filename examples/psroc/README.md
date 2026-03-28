# PSROC Example

End-to-end example of the PSROC pipeline for evaluating variant pathogenicity prediction scores.

For full documentation, see the [PSROC docs](https://bigbio.github.io/hvantk/tools/psroc/) and [PSROC examples guide](https://bigbio.github.io/hvantk/examples/psroc/).

## Contents

- **`run_psroc_example.py`** - Complete workflow using synthetic test data

## Quick Start

```bash
# Run the example (uses synthetic test data from hvantk/tests/testdata/psroc/)
python examples/psroc/run_psroc_example.py

# View results
ls examples/psroc/results/
```

## Test Data

Uses synthetic data from `hvantk/tests/testdata/psroc/`:

- 50 pathogenic + 40 benign + 10 VUS variants across BRCA1, BRCA2, TP53
- 4 prediction scores: REVEL (~0.95 AUC), CADD (~0.90), MetaLR (~0.75), VEST4 (excluded, >30% missing)

## Expected Outputs

```
results/
├── plots/
│   ├── psroc_example_roc_curves.png
│   ├── psroc_example_auc_comparison.png
│   ├── psroc_example_missingness.png
│   └── psroc_example_dashboard.png
├── psroc_example_metrics.json
├── psroc_example_missingness.json
└── psroc_example_annotated.tsv
```
