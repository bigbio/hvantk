# PS-ROC example

Runs the PS-ROC pipeline end to end on synthetic data, so you can see every output without
real ClinVar or dbNSFP tables. Documentation: the
[PS-ROC reference](https://bigbio.github.io/hvantk/tools/psroc/) and the
[example walkthrough](https://bigbio.github.io/hvantk/examples/psroc/).

## Contents

- `run_psroc_example.py`: builds two Hail tables from the synthetic TSVs in
  `hvantk/tests/testdata/psroc/`, runs `PSROCPipeline` on BRCA1, BRCA2 and TP53 with four
  scores and the default settings, and prints the results.

## Run

Requires Hail (Java 11) and the `psroc` extra (`poetry install --extras psroc`):

```bash
python examples/psroc/run_psroc_example.py --output-dir out/psroc_example
```

The script finds the test data relative to its own location, so it runs from any
directory. `--output-dir` defaults to `/tmp/psroc_example`.

## Data

- 100 variants in BRCA1, BRCA2 and TP53: 50 pathogenic, 40 benign and 10 uncertain, all
  with 1- or 2-star review status.
- Four scores, missing for some variants: `CADD_phred` (7), `REVEL_score` (3),
  `MetaLR_score` (10) and `VEST4_score` (42).

## Expected results

The 10 uncertain variants stay in the annotated table as excluded and enter neither the
missingness check nor the ROC, so the missing rates are over the 90 pathogenic and benign
variants (with the number of them missing the score).

| Score | Missing | AUC | Status |
| --- | --- | --- | --- |
| `REVEL_score` | 3.3% (3) | 0.998 | included |
| `CADD_phred` | 7.8% (7) | 0.979 | included |
| `MetaLR_score` | 10.0% (9) | 0.896 | included |
| `VEST4_score` | 42.2% (38) | | excluded (above the 30% threshold) |

```text
out/psroc_example/
├── tables/                  # clinvar.ht and dbnsfp.ht built from the TSVs
└── results/
    ├── psroc_example_metrics.json
    ├── psroc_example_missingness.json
    ├── psroc_example_annotated.ht/
    ├── psroc_example_annotated.tsv
    ├── plots/
    │   ├── psroc_example_roc_curves.png
    │   ├── psroc_example_auc_comparison.png
    │   ├── psroc_example_missingness.png
    │   └── psroc_example_dashboard.png
    ├── logs/                # this run's log file
    └── .pipeline_state.json
```
