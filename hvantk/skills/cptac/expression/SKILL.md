---
name: hvantk:resource-cptac-expression
description: Build an AnnData (samples x genes) from long-format CPTAC protein-expression matrices plus sample metadata.
status: provisional
backend: anndata
domain: proteomics
---

# CPTAC protein expression resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. Builder + drift-probe placeholder are migrated under the plugin folder. Builder round-trip snapshots are NOT yet seeded; the fixture directory exists but is empty. Existing builder coverage lives in `hvantk/tests/test_expression_builders_anndata.py` (`TestBuildCptacAd`).
- **In scope:** one long-format CPTAC protein-expression TSV/CSV at a time, paired with a sample metadata file. Output is an AnnData (`.h5ad`) keyed `samples x genes` with `float32` `X`. The user is responsible for staging the inputs -- there is no automated downloader for the expression matrix today (the `cptac` Python package fetches into memory for the phospho path; mirroring that for expression is a follow-up).
- **Out of scope:** phosphoproteomics (sibling dataset `cptac:phospho`); per-cancer-type batch builds (callers loop themselves); any Hail-Table or wide-format representation.

## 2. Source identity

- **Provider:** Clinical Proteomic Tumor Analysis Consortium (CPTAC) (<https://proteomics.cancer.gov/programs/cptac>).
- **Catalog entry:** TODO. No CPTAC entry exists yet in any plugin's `catalog/datasets.json` (verify with `hvantk catalog search cptac`). The cancer-type enum + class map live in `hvantk/ptm/constants.py` (`CPTAC_CANCER_TYPES`, `CPTAC_CANCER_CLASS_MAP`). When a dedicated catalog entry lands under `hvantk/skills/cptac/catalog/datasets.json`, point this plugin's `source.catalog_ref` at it and remove this TODO.

## 3. Backend choice + reasoning

**`backend: anndata`, `domain: proteomics`.** Per `_conventions` § 3, AnnData is the right fit for sparse-ish sample x feature matrices keyed by sample identifier. CPTAC protein expression is `samples x genes`, fits in memory for any single cancer cohort, and pairs naturally with sample-level clinical metadata in `obs`. A Hail Table would force a Spark session for a sub-100 MB matrix and lose the `obs` ergonomics.

## 4. Raw format & gotchas

The builder accepts a long-format table (one row per `(sample, gene)` observation) and pivots it to wide:

- Required columns (defaults): `GeneID`, `SampleID`, `Expression`. Override via `gene_id_col`, `sample_id_col`, `expression_col`.
- Optional: `Gene Name` (column name configurable via `gene_name_col`) -- if present, it is propagated to `var`.
- Delimiter is auto-detected by `pd.read_csv(sep=None, engine="python")`.

Metadata file is read with the same auto-detect logic and indexed by `SampleID` (override via `sample_id_col`). Any extra columns become `obs` columns via a left join, so samples in the metadata that are absent from the expression matrix are silently dropped (the join is on `wide.index`).

`X` is coerced to `np.float32`. Missing values from the pivot become `NaN` and are preserved in `X`.

## 5. Output contract

- **File:** AnnData `.h5ad` at the path passed via `output_path` (callers' choice -- there is no fixed naming).
- **Shape:** `(n_samples, n_genes)` -- samples in `obs`, genes in `var`.
- **`X`:** `float32`, missing values as `NaN`.
- **`obs`:** indexed by sample id; columns mirror the metadata file (after `set_index`).
- **`var`:** indexed by gene id (name = `gene_id_col`); contains `gene_name_col` if it was present in the input.
- **`uns["column_summary"]`:** added by `annotate_column_summary_ad`.
- **Provenance:** stamped on the returned `ExpressionMatrix` via `ctx.provenance(schema_id="cptac-expression-v1")`; persisted by the platform as a sidecar `.provenance.json` next to the `.h5ad`.

## 6. hvantk integration points

- **Builder:** `build_cptac_ad` in `hvantk/skills/cptac/expression/builder.py`.
- **Shared helpers:** `create_anndata_from_cptac_long`, `_parse_site_id` in `hvantk/skills/cptac/shared/cptac.py` (shared with the sibling phospho builder).
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/cptac/expression/drift_probe.py` (fingerprints the installed `cptac` Python package version).
- **Plugin manifest:** `hvantk/skills/cptac/plugin.yaml` (compound dataset key `cptac:expression`).
- **Tests:** parser/helper coverage in `hvantk/skills/cptac/expression/tests/` (drift-probe sanity); builder coverage in `hvantk/tests/test_expression_builders_anndata.py` (`TestBuildCptacAd`).
- **CLI:** end-to-end `hvantk reprocess cptac:expression` is not yet wired — the Phase B `build_cptac_expression` builder expects `parsed_input` to be a `{"expression": <tsv>, "metadata": <tsv>}` dict produced by a `lifecycle.parse` stage that the manifest does not yet declare. Build via the Python API today (`build_cptac_expression(parsed_input, ctx, gene_id_col=…, sample_id_col=…, expression_col=…)`, or the legacy `build_cptac_ad`).

## 7. Workflow steps

When invoked to build a CPTAC protein-expression AnnData:

1. **Stage inputs.** Long-format expression TSV/CSV + sample metadata TSV/CSV. (Today the user produces these manually; an automated downloader for expression is a follow-up -- see `cptac:phospho` for the in-package fetcher pattern.)
2. **Build.** Import `build_cptac_expression` (or the legacy `build_cptac_ad`) directly — see § 6 on why `hvantk reprocess cptac:expression` is not wired yet.
3. **Validate.** `pytest hvantk/skills/cptac/expression/tests` (drift-probe sanity) + `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildCptacAd` (builder round-trip).

## 8. Update playbook

CPTAC ships new cohorts and re-processed runs via the `cptac` Python package. When a new package release lands:

1. Run the drift probe: `python -c "from hvantk.skills.cptac.expression.drift_probe import fetch_fingerprint; print(fetch_fingerprint())"`. A change in `source_version` (the installed `cptac` package version) is the trigger.
2. Re-stage / re-pull the expression matrix for any cohort of interest using the package's `get_proteomics` API (or whatever the new release exposes).
3. Re-build and spot-check sample / gene counts against the previous run.
4. Re-regenerate `tests/drift_fingerprint.json` once the package version is updated.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/cptac/expression/tests/testdata/raw/cptac-expression/` (directory present, not yet seeded -- existing builder tests in `hvantk/tests/test_expression_builders_anndata.py` build mock inputs on the fly).
- **schema_snapshot:** `hvantk/skills/cptac/expression/tests/snapshots/schema.json` (TODO -- created on first `--regenerate-snapshots` run).
- **row_snapshot:** `hvantk/skills/cptac/expression/tests/snapshots/sample_rows.json` (TODO -- same).
- **test_command:** `pytest hvantk/skills/cptac/expression/tests`.
- **drift_fingerprint:** `hvantk/skills/cptac/expression/tests/drift_fingerprint.json` (placeholder shape; refresh via the update playbook).

The plugin manifest declares these paths so the loader contract holds. Drift-probe sanity test passes today; the builder round-trip snapshot is the gap to close in a follow-up.

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run in a hail-enabled environment, use
> `pytest hvantk/skills/cptac/expression/tests/test_builder.py --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
