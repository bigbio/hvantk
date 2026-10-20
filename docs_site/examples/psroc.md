# PS-ROC example

This example runs PS-ROC end to end on a small synthetic ClinVar and dbNSFP pair that ships
with the tests, so you can see every output without building the real tables. For the
options, the method and the output formats, see the [PS-ROC reference](../tools/psroc.md).

## Run the example

The script needs Hail and the `psroc` extra. It builds two Hail tables from the synthetic
TSVs, runs the pipeline through the Python API and prints the results:

```bash
python examples/psroc/run_psroc_example.py --output-dir out/psroc_example
```

The script finds its data relative to its own location, so you can run it from any
directory. Without `--output-dir` it writes to `/tmp/psroc_example`.

## Data

The two files are in `hvantk/tests/testdata/psroc/`:

- `synthetic_clinvar.tsv`: 100 variants in BRCA1 (35), BRCA2 (33) and TP53 (32), with 50
  pathogenic, 40 benign and 10 `Uncertain_significance` labels. Every review status has 1
  or 2 stars.
- `synthetic_dbnsfp.tsv`: `CADD_phred`, `REVEL_score`, `MetaLR_score` and `VEST4_score` for
  the same variants, missing for 7, 3, 10 and 42 of them.

The script nests `CLNSIG`, `CLNREVSTAT` and `GENEINFO` (plain strings) under an `info`
struct, as in a table built from the ClinVar VCF, and keys both tables by
`(locus, alleles)`.

## Expected results

With the default settings (`min_stars` 1, `max_missingness` 0.3, Youden threshold, 2000
bootstrap resamples), all 100 variants pass the review-status filter. The 10 uncertain
variants stay in the annotated table as excluded but enter neither the missingness check
nor the ROC, so the missing rates below are over the 90 pathogenic and benign variants
(with the number of them missing the score):

| Score | Missing | AUC | Status |
| --- | --- | --- | --- |
| `REVEL_score` | 3.3% (3) | 0.998 | included |
| `CADD_phred` | 7.8% (7) | 0.979 | included |
| `MetaLR_score` | 10.0% (9) | 0.896 | included |
| `VEST4_score` | 42.2% (38) | | excluded (above 30%) |

The script also prints each AUC's 95% confidence interval, the optimal threshold and its
sensitivity and specificity.

## Outputs

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

## The same run from the CLI

On tables built with `hvantk reprocess`, this command runs the same analysis:

```bash
hvantk psroc \
  --genes BRCA1,BRCA2,TP53 \
  --clinvar-ht data/clinvar.ht \
  --dbnsfp-ht data/dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score,VEST4_score" \
  --output-prefix psroc_example \
  --export-tsv \
  --output-dir results/psroc_example
```

The [`examples/psroc/`](https://github.com/bigbio/hvantk/tree/main/examples/psroc/)
directory holds the script and three sample dashboards from a `--gene-sets` run on real
ClinVar and dbNSFP tables.
