# gevir

GeVIR (Gene Vulnerability and Intolerance Rank) provides gene-level metrics
quantifying the intolerance of human genes to functional variants. The resource
was published in PMID 31873297 (Chen et al. 2020, Nature Communications) and
provides ranks and scores for ~18,000 protein-coding genes based on the spatial
distribution of de novo mutations.

## Dataset

- `gevir:metrics` — per-gene vulnerability and intolerance rank metrics, keyed by gene_id

The builder is `build_gevir_metrics` in `hvantk/skills/gevir/builder.py`, with
signature `(parsed_input, ctx, **params) -> AnnotationTable`. It imports the
GeVIR TSV inline via `hl.import_table(..., key="gene_id")`, optionally selects a
`fields` subset, and wraps the result with provenance under schema
`gevir-metrics-v1`. The dataset is resolved from `plugin.yaml` by the plugin
loader (`hvantk/core/plugin/loader.py`) and built through `run_builder_for_spec`
(`hvantk/core/plugin/run_builder.py`); there is no separate builder registry.

## Build

```bash
hvantk reprocess gevir:metrics \
  --raw-dir <dir-containing-gevir-tsv> \
  --output <out.ht>
```

Optional plugin args (e.g. field selection) can be passed with
`--plugin-arg fields=...`.

## Notes

The drift probe (`hvantk/skills/gevir/drift_probe.py`, `fetch_fingerprint`) is a
stub; a real probe should be implemented in a follow-up. No downloader is
implemented; upstream files are expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_id` (Ensembl gene ID).
All fields from the upstream TSV are imported with type imputation.

## Tests

A conformance test in `tests/test_gevir.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz`.
