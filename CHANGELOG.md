# Changelog

## Unreleased

### Added

- Declarative feature selection for `hvantk rerank` (Python API: `Config.selection`). Filters run within each axis — univariate AUC with within-axis BH-FDR, then Spearman redundancy — re-fitted inside every cross-validation fold on the training slice only, so the reported ΔAUC is not inflated by selection that has seen the held-out labels. A third RFECV step is available but **off by default** (`SelectionPolicy(wrapper="rfecv")`): across four real cohorts it eliminated columns almost exclusively in the one with the fewest positives, and pruned the ablation baseline axis, so it needs an out-of-fold outcome comparison before it can be trusted by default. `Config.selection = None` (the default) reproduces the previous code path exactly, and the CLI is unchanged.
- `rerank_arms(config)` runs each analysis as two arms, `clean` and `all`, over identical folds. `clean` (columns with no provenance conflict against the label source) is the headline; `all` adds conflicted and undeclared columns so the circularity channel is a measured number rather than an assumption. `RerankResult.selection` carries the per-fold selection frequency, the global-pass feature list, and both nested and global AUCs.
- Plugin manifests may declare per-predictor training provenance: an optional `scores: {<column>: {trained_on: [...]}}` block per dataset. `hvantk/skills/dbnsfp/plugin.yaml` declares it for 55 of its 57 rankscore predictors. An omitted score means unknown and is never treated as clean.
- Plugin system for data-provider adapters. Each provider now lives in a single folder under `hvantk/skills/<provider>/` with a `plugin.yaml` manifest, builder code, drift probe, downloader CLI, and tests. The loader auto-discovers plugins from the in-tree filesystem and Python entry points.
- `hvantk plugins {list,describe,errors,validate}` commands for inspecting the registry.
- `hvantk drift <provider:dataset>` for upstream-drift detection against committed expected fingerprints.
- `hvantk reprocess <provider:dataset>` for chaining download -> parse -> build -> drift-check from a single command.
- 13 migrated provider plugins: clingen (gene-disease), clinvar, cptac (expression + phospho), expression-atlas, gencc (submissions), gtex-eqtl, gwas-catalog, hgnc, insider, msigdb, peptideatlas (phospho), ucsc-cellbrowser (default / adult-ctx / dev-ctx), uniprot-ptm (sites).
- Scheduled CI workflow (`.github/workflows/drift.yml`) that runs `hvantk drift --all --json` daily and opens a draft PR per drifted plugin with the regenerated fingerprint pre-committed.

### Changed

- Package restructured into 4 purpose-driven roofs: `core/` (platform models, utilities, plugin/tool runtime, streamers, transient builders), `algorithms/` (analytical computation: ptm, psroc, qtlcascade, enrichex, hgc, ancestry, annotation, visualization, expression, statistics, training_sets), `skills/` (data ingestion plugins), `tools/` (CLI surface). Inside `core/` there are now sub-packages `models/`, `utils/`, `streamers/`, `plugin/`, `tool/`, `builders/` so adding a new format helper has one obvious home. One-way dependency rule (`skills/`, `tools/` → `algorithms/` → `core/`) is enforced by `hvantk/tests/test_dependency_directions.py`. `hvantk/data/`, `hvantk/utils/`, `hvantk/tables/`, and 8 top-level algorithm dirs (`hvantk/{ptm,psroc,qtlcascade,enrichex,hgc,ancestry,annotation,visualization}/`) are gone. `ClinVarStreamer` no longer imports from `hvantk.skills.clinvar.builder` — it accepts a pre-built Hail Table via its constructor.
- Registry keys for migrated providers use compound `provider:dataset` form. Recipe JSONs and any custom callers should update from bare names (e.g., `clinvar`) to compound (`clinvar:variants`). The legacy `hvantk mktable` / `hvantk mkmatrix` CLI surfaces have been retired; data builds now go through `hvantk reprocess <provider>:<dataset>` with `--plugin-arg key=value` for builder kwargs.
- Plugin manifests gain an optional `catalog: <path>` field pointing at a per-plugin `catalog/datasets.json`. `unified_registry.HvantkRegistry` now aggregates per-plugin catalogs from the plugin loader in addition to the legacy `resources/registry/genomics/datasets.json`.
- Per-domain catalogs `resources/registry/{transcriptomics,proteomics,epigenomics}/datasets.json` are removed; their entries now live inside each owning plugin's `catalog/datasets.json` (expression-atlas, ucsc-cellbrowser). `registry/genomics/datasets.json` is intentionally retained until orphan entries (dbNSFP, gnomad-metrics, ensembl-gene, gevir, cosmic-cgc) gain owning plugins.
- `hvantk catalog` CLI rewritten to read per-plugin catalogs via `HvantkRegistry`. New subcommands: `list` (with `--omics-type` / `--data-source` / `--organism` filters), `show`, `stats`, `search`. The legacy `catalog build` subcommand is removed; use `hvantk reprocess <provider:dataset>` instead.

### Removed

- Per-provider downloader modules under `hvantk/commands/*_downloader.py` for migrated providers (moved into their plugin folder's `cli.py`).
- Per-provider dataset classes under `hvantk/datasets/*_datasets.py` for migrated providers (moved into `hvantk/skills/<provider>/shared/`).
- Per-provider builder functions in `hvantk/tables/table_builders.py` and `matrix_builders.py` for migrated providers (moved into `hvantk/skills/<provider>/[<dataset>/]builder.py`).
- `hvantk/resources/generate_catalog.py` (regenerated the now-removed per-domain `datasets.json` files). Catalog regeneration is now a per-plugin concern; if a maintainer needs a packaged regenerator in the future it should live alongside each plugin's `catalog/datasets.json`.
- `hvantk/resources/catalog.yaml` (auto-generated summary file pointing at deleted per-domain JSON files). Equivalent information is available on demand via `hvantk catalog stats`.

### Known gaps before first stable release

- The following plugins reference snapshot files (`schema.json`, `sample_rows.json`)
  in their `plugin.yaml` manifests that have not yet been seeded on disk:
  `clingen`, `gencc`, `hgnc`, `uniprot-ptm`, `expression-atlas`, `peptideatlas:phospho`,
  `cptac:expression`, and `cptac:phospho`. The first hail-enabled CI run with
  `--regenerate-snapshots` will bootstrap them. All `ucsc-cellbrowser` variants
  (`default`, `adult-ctx`, `dev-ctx`) already have populated snapshot dirs.
