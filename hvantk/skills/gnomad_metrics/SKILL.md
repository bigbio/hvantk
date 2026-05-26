# gnomad-metrics

The gnomAD (Genome Aggregation Database) constraint gene metrics table
provides per-gene constraint scores including pLI (probability of loss-of-
function intolerance), oe_lof (observed/expected loss-of-function ratio),
and related statistics derived from the gnomAD v2.1.1 cohort. These metrics
are widely used to assess gene-level intolerance to variation.

Upstream: https://gnomad.broadinstitute.org/downloads

## Dataset

- `gnomad-metrics:metrics` — per-gene constraint metrics, keyed by gene_id

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The drift probe is a stub; a real probe should be implemented in a
follow-up. The downloader is not implemented; upstream files are
expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_id` (Ensembl gene ID).
All fields from the upstream TSV are imported with type imputation.

## Tests

A conformance test in `tests/test_gnomad_metrics.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz`.
