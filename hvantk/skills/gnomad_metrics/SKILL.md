# gnomad-metrics

The gnomAD (Genome Aggregation Database) constraint gene metrics table
provides per-gene constraint scores including pLI (probability of loss-of-
function intolerance), oe_lof (observed/expected loss-of-function ratio),
and related statistics derived from the gnomAD v2.1.1 cohort. These metrics
are widely used to assess gene-level intolerance to variation.

Upstream: https://gnomad.broadinstitute.org/downloads

## Dataset

- `gnomad-metrics:metrics` — per-gene constraint metrics, keyed by gene_id

## Build

```bash
hvantk reprocess gnomad-metrics:metrics \
  --raw-dir <dir-with-lof_metrics-tsv> \
  --output <out.ht>
```

Optional field selection can be passed through as a plugin arg, e.g.
`--plugin-arg fields='["gene_id", "pLI", "oe_lof"]'`.

The builder is `build_gnomad_metrics_metrics` in
`hvantk/skills/gnomad_metrics/builder.py`, with signature
`(parsed_input, ctx, **params) -> AnnotationTable`. It imports the gnomAD
lof_metrics TSV inline via `hl.import_table(..., impute=True, key="gene_id")`,
optionally selects `params["fields"]`, and wraps the result with provenance
(`schema_id="gnomad-metrics-v1"`). The plugin loader auto-resolves this
dataset from `plugin.yaml`; top-level builds run through
`run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).

## Phase K notes

This plugin was promoted as part of Phase K of the data-model platform
refactor. The drift probe (`drift_probe.py`, `fetch_fingerprint`) is a
stub; a real probe should be implemented in a follow-up. The downloader is
not implemented; upstream files are expected to be externally materialized
for now.

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
