# gnomad-metrics

The gnomAD (Genome Aggregation Database) constraint gene metrics table
provides per-gene constraint scores including pLI (probability of loss-of-
function intolerance), oe_lof (observed/expected loss-of-function ratio),
and related statistics derived from the gnomAD v2.1.1 cohort. These metrics
are widely used to assess gene-level intolerance to variation.

Upstream: https://gnomad.broadinstitute.org/downloads

## Dataset

- `gnomad-metrics:metrics` — per-gene constraint metrics, keyed by gene_id

## Download

The constraint tables are small and public, so the plugin ships a downloader
(`hvantk download gnomad-metrics`). Two releases are supported:

- **v2.1.1** (GRCh37, ~4.6 MB) — the default; `by_gene` (keyed by `gene_id`) or
  `by_transcript`. This is what hvantk standardises on.
- **v4.0** (GRCh38, ~82 MB) — `constraint_metrics`; per-transcript rows with
  dotted column names (`lof.pLI`, `lof.oe_ci.upper` = LOEUF) and **no `gene_id`**
  column, so build it with `--plugin-arg key=transcript`. (gnomAD did not
  re-release constraint for v4.1, so v4.0 is the newest.)

```bash
hvantk download gnomad-metrics --output data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
hvantk download gnomad-metrics --version v4.0 --output data/gnomad/gnomad.v4.0.constraint_metrics.tsv
```

## Build

`hvantk reprocess` runs download + build end-to-end (the downloader drops the
file into `--raw-dir`; the builder resolves it):

```bash
hvantk reprocess gnomad-metrics:metrics \
  --raw-dir data/gnomad/ \
  --output gnomad_metrics.ht
```

Building the Hail Table needs Java 11 and adequate driver heap (locally:
`PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"`). Optional field
selection: `--plugin-arg fields='["gene_id", "pLI", "oe_lof"]'`. For v4.0, pass
`--plugin-arg key=transcript` (no `gene_id` column).

The builder is `build_gnomad_metrics_metrics` in
`hvantk/skills/gnomad_metrics/builder.py`, with signature
`(parsed_input, ctx, **params) -> AnnotationTable`. `parsed_input` may be a file
path or a `raw_dir` (the builder resolves the `*.bgz`/`*.tsv` inside it, so
`reprocess` works end-to-end). It imports the table via
`hl.import_table(..., impute=True, key=params.get("key", "gene_id"))`, optionally
selects `params["fields"]`, and wraps the result with provenance
(`schema_id="gnomad-metrics-v1"`). The plugin loader auto-resolves this dataset
from `plugin.yaml`; top-level builds run through `run_builder_for_spec`
(`hvantk/core/plugin/run_builder.py`).

## Phase K notes

This plugin was promoted as part of Phase K of the data-model platform
refactor. The drift probe (`drift_probe.py`, `fetch_fingerprint`) is a
stub; a real probe should be implemented in a follow-up. The downloader
(`cli.py`, `download_dataset` / `download_cmd`) fetches the constraint tables
from the public gnomAD GCS bucket — see the Download section above.

## Schema

Field documentation TBD. The table is keyed by `gene_id` (Ensembl gene ID).
All fields from the upstream TSV are imported with type imputation.

## Tests

A conformance test in
`hvantk/skills/gnomad_metrics/tests/test_gnomad_metrics.py` exercises the
build via `run_builder_for_spec` against the gnomAD lof_metrics fixture
(`gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz`). The `tests:` block in
`plugin.yaml` declares the plugin-relative fixture at
`tests/testdata/raw/gnomad`.

Run with: `pytest hvantk/skills/gnomad_metrics/tests`
