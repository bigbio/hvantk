# gevir

GeVIR (Gene Vulnerability and Intolerance Rank) provides gene-level metrics
quantifying the intolerance of human genes to functional variants. The resource
was published in PMID 31873297 (Chen et al. 2020, Nature Communications) and
provides ranks and scores for ~18,000 protein-coding genes based on the spatial
distribution of de novo mutations.

## Dataset

- `gevir:metrics` — per-gene vulnerability and intolerance rank metrics, keyed by gene_id

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

A conformance test in `tests/test_gevir.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz`.
