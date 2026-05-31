---
name: hvantk:resource-ucsc-cellbrowser
description: Onboard, build, or update a UCSC Cell Browser single-cell expression resource for hvantk
status: provisional
backend: anndata
domain: expression
---

# UCSC Cell Browser resource skill

## 1. Status & scope

This skill covers BUILD and UPDATE of a UCSC Cell Browser single-cell AnnData (`.h5ad`) from raw expression + metadata TSVs. It does NOT cover download (see `hvantk/skills/ucsc_cellbrowser/cli.py`), scanpy-based downstream analysis, or atlas-level merging across multiple UCSC collections.

## 2. Source identity

Provider metadata (URLs, collections, organisms, per-collection sample counts) lives in `hvantk/skills/ucsc_cellbrowser/catalog/datasets.json` (per-collection key is `accession`, e.g., `cortex-dev`, `ad-multi-region`). Browse via `hvantk catalog list --data-source UCSC` or `hvantk catalog show cortex-dev`. Read those files; do not restate.

Stable note (not in catalog): UCSC Cell Browser publishes per-collection (one accession at a time), not on a global monthly cadence like ClinVar. Each collection is curated independently and may be updated without warning. There is no global "release version" string — versioning is implicit in the collection's last-modified timestamp captured in the registry.

## 3. Backend choice + reasoning

AnnData (`.h5ad`). Single-cell expression is sparse cells × genes; AnnData's CSR storage and the surrounding scanpy ecosystem are the natural fit. A Hail MatrixTable would force materialization of zeros and lose access to scanpy / cell-level QC tooling. Conventions § 3 ("expression sparse, single-cell → anndata h5ad") applies directly.

## 4. Raw format & gotchas

UCSC publishes two TSVs per collection:

- **Expression matrix TSV**: rows = genes, columns = cells. The first column is a gene identifier (symbol, Ensembl ID, or `<symbol>|<id>` pipe-encoded). The header's first column name is collection-dependent (commonly `gene` or `geneId`) — the builder default is `gene_column="gene"`. Tab-delimited. May be plain `.tsv` or gzipped (`.gz`/`.bgz`); the loader auto-detects via extension.
- **Metadata TSV**: rows = cells, columns = arbitrary cell-level annotations (commonly `celltype`, `cluster`, `donor`). Tab-delimited. Cell IDs must match the expression matrix's column header.

Gotchas:

- **Pipe-encoded gene names.** Some UCSC collections format gene IDs as `<symbol>|<ensembl_id>`. The builder's default `split_gene_field=True` keeps only the part before the `|`. If a collection ships fully-qualified IDs, set `split_gene_field=False` to retain the full string.
- **Cell ID alignment.** The expression matrix header (cell IDs) and the metadata TSV's cell-ID column must align. Mismatches surface only when `obs` is constructed; the builder does not silently re-key.
- **Backed-vs-in-memory auto-selection.** `build_ucsc_cellbrowser(backed=None, ...)` chooses based on input size compared with the threshold defined at `hvantk/skills/ucsc_cellbrowser/builder.py:BACKED_BUILDER_THRESHOLD_BYTES` (default: 1 GiB). Backed mode requires `backed_output_path` (passed via `--plugin-arg backed_output_path=<path>`, same path as `--output`). Backed-mode `obs` writes are guarded by `coerce_obs_for_h5ad()` because anndata's vlen-string HDF5 writer chokes on object-dtype columns mixing strings and `NaN`; the in-memory path applies the same coercion before saving.
- **Backed mode skips column-summary annotation.** `annotate_column_summary_ad` is only called in the in-memory branch — for atlas-scale inputs it would force a full X scan. When auditing a backed build, do not assume the per-column summary fields exist on the returned handle.
- **Streaming row iterator.** `_iter_ucsc_rows` (in `hvantk/skills/ucsc_cellbrowser/shared/ucsc.py`) wraps a file handle that closes on generator exhaustion. Partial iteration leaks the handle until GC; full iteration is the contract.

## 5. Output contract

The builder returns an `ExpressionMatrix` wrapping an `AnnData`; the platform persists it to the `--output` `.h5ad` path via `artifact.save()`. Schema is defined by `hvantk/skills/ucsc_cellbrowser/tests/snapshots/ucsc-cellbrowser/schema.json` (canonical). Human summary: AnnData with `obs` indexed by cell ID and carrying cell metadata (rows = cells), `var` indexed by gene ID (rows = genes), `X` a `float32` sparse matrix (CSR for the in-memory path; CSC for the backed path). Provenance is stamped on the returned `ExpressionMatrix` via `ctx.provenance(schema_id=...)` and persisted by the platform as a sidecar `.provenance.json`.

## 6. hvantk integration points

- Builder entry point: `build_ucsc_cellbrowser` in `hvantk/skills/ucsc_cellbrowser/builder.py`. Signature is `build_ucsc_cellbrowser(parsed_input, ctx, *, gene_column, delimiter, split_gene_field, chunk_size, backed, column_batch, backed_output_path, **params) -> ExpressionMatrix`; `parsed_input` carries the `expression_matrix` and `metadata` paths. One function serves all three datasets (`default`, `adult-ctx`, `dev-ctx`); the per-dataset `schema_id` is resolved from `ctx.dataset`.
- In-memory helper: `create_anndata_from_ucsc_matrix` in `hvantk/skills/ucsc_cellbrowser/shared/ucsc.py`.
- Backed-write helper: `build_ucsc_atlas_backed` in `hvantk/skills/ucsc_cellbrowser/shared/ucsc.py`.
- Metadata loader: `load_ucsc_metadata` in `hvantk/skills/ucsc_cellbrowser/shared/ucsc.py`.
- Provenance: stamped by the platform via `ctx.provenance(...)` (see `hvantk/core/plugin/run_builder.py`). Save / summary helpers: `save_anndata` in `hvantk/core/io/anndata_io.py`, `annotate_column_summary_ad` in `hvantk/core/models/anndata_utils.py`.
- CLI: `hvantk reprocess ucsc-cellbrowser:<dataset> --raw-dir <dir> --output <path>.h5ad` (dataset is one of `default`, `adult-ctx`, `dev-ctx`). Builder kwargs flow through `--plugin-arg key=value`.
- Registry: the plugin loader (`hvantk/core/plugin/loader.py`) auto-resolves each dataset from the plugin manifest at `hvantk/skills/ucsc_cellbrowser/plugin.yaml` via `get_registry().get_dataset("ucsc-cellbrowser:<dataset>")` (one of `default`, `adult-ctx`, `dev-ctx`). Top-level builds run through `run_builder_for_spec` in `hvantk/core/plugin/run_builder.py`. There is no `MATRIX_BUILDERS` registry or `registry.py`.
- Test: `hvantk/skills/ucsc_cellbrowser/tests/test_builder.py`.

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update:

1. Verify `anndata`, `scanpy`, and `scipy` are importable. The current SessionStart hook (`hvantk/skills/_hooks/check-hail.sh`) only checks Hail + Java — anndata/scanpy availability is a known gap; verify manually with `python -c "import anndata, scanpy, scipy"` until the hook is extended.
2. Confirm the expression input is a UCSC matrix: read the first line; the leading column should be a gene identifier (e.g., `gene`/`geneId`/`Symbol`), not a cell ID. The remaining header fields should be cell IDs that match the metadata TSV's cell-ID column.
3. Confirm the metadata input: read the first line of the metadata TSV; the cell-ID column is whatever `UCSC_CELL_ID_COLUMN` resolves to (see `hvantk/skills/ucsc_cellbrowser/shared/constants.py`); other columns become `obs` fields.
4. Invoke `hvantk reprocess ucsc-cellbrowser:<dataset> --raw-dir <dir> --output <path>.h5ad`, passing builder kwargs through `--plugin-arg key=value` (e.g. `--plugin-arg gene_column=geneId --plugin-arg split_gene_field=false`). The platform calls `build_ucsc_cellbrowser(parsed_input, ctx, ...)`; it auto-selects backed mode based on file size; pass `--plugin-arg backed=true` to force.
5. If backed mode is selected (auto or forced), `--plugin-arg backed_output_path=<path>` is required (use the same path as `--output`). Division of responsibility: `build_ucsc_atlas_backed` writes CSC columns directly to disk (no return value of interest); `build_ucsc_cellbrowser` re-opens the written file via `anndata.read_h5ad(backed_output_path, backed='r')` and returns that read-backed handle wrapped in an `ExpressionMatrix`. Per § 4, do not rely on `annotate_column_summary_ad`-derived fields in this branch — the backed path skips that pass.
6. Backed-path operation order is fixed: coerce `obs` via `coerce_obs_for_h5ad` → pre-flight via `_probe_write_elem` (catches obs-serialization failures BEFORE streaming any X columns) → stream X column-batches via `write_elem` + `sparse_dataset(...).append(...)` → finalize by writing `/obs`, `/var`, `/uns`. Reordering risks committing a partial atlas to disk before discovering an obs failure.
7. Run validation: `pytest hvantk/skills/ucsc_cellbrowser/tests/test_builder.py` (this test is NOT marked `@pytest.mark.hail` — it does not require a Hail session).
8. Report: schema diff, head sample-row diff, test pass/fail.

## 8. Update playbook

UCSC datasets are per-collection. To refresh a collection's build:

1. Re-fetch the collection's expression and metadata TSVs via `hvantk/skills/ucsc_cellbrowser/cli.py` (or replace the fixture by hand if regenerating from a stored slice).
2. If updating the round-trip fixture, replace `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser/expression_matrix.tsv` and `metadata.tsv`. The current pilot fixture is a 12 cells × 250 genes slice from `asp_2019_celltype_summary.h5ad`; keep slices small but representative (mix of cell types, sparse rows, at least one gene with pipe-encoded ID if relevant).
3. Run snapshot regeneration: `pytest hvantk/skills/ucsc_cellbrowser/tests/test_builder.py --regenerate-snapshots`. (No `-m hail`; this builder is not Hail-backed.)
4. Inspect the snapshot diff:
   - **Expected diff** — new `obs_columns` (added metadata field), changed `n_obs` / `n_vars` (different fixture shape), new `layers` if a downstream caller starts populating one. Commit the regenerated snapshots with explanation.
   - **Unexpected diff** — `X_dtype` shifts away from `float32` (likely a cast regression), `X_format` changes between `csr_matrix` / `csc_matrix` for the same code path, missing or renamed `obs_columns` (data loss). STOP and investigate before committing.
5. If a parsing gotcha was introduced (e.g., new metadata encoding, different gene-ID separator, new column with mixed-NaN dtype), add it to § 4.
6. Open PR; reviewer focus is the schema + sample-row diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser/expression_matrix.tsv` (primary) and `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser/metadata.tsv` (secondary).
- `schema_snapshot`: `hvantk/skills/ucsc_cellbrowser/tests/snapshots/ucsc-cellbrowser/schema.json`
- `row_snapshot`: `hvantk/skills/ucsc_cellbrowser/tests/snapshots/ucsc-cellbrowser/sample_rows.json`
- `test_command`: `pytest hvantk/skills/ucsc_cellbrowser/tests/test_builder.py`

## 10. Onboarding a new dataset within UCSC Cell Browser

When a new collection appears, follow this triage before deciding what to add:

1. **Inspect metadata TSV header.** Identify the cell-ID column and classify the others as categorical, numeric, or ontology-coded.
2. **Inspect expression TSV header.** Determine the gene-column name (`gene`, `geneId`, `geneSymbol`, etc.). Note pipe-encoded gene names if present.
3. **Sniff a sample of X values.** Integer-typed values (large dynamic range) indicate raw counts; small float values (~0–10) indicate normalized.
4. **Add a catalog entry** to `hvantk/skills/ucsc_cellbrowser/catalog/datasets.json` (and a `datasets:` entry in `hvantk/skills/ucsc_cellbrowser/plugin.yaml` if it is a new compound dataset) with:
   - `accession`, `data_source: "UCSC"`, `url`, `last_updated`
   - Per-dataset build params: `gene_column`, `split_gene_field`, `pre_normalized`, `expected_obs_columns` (for drift detection)
5. **Decide whether to add a snapshot test.** Add one if the dataset's obs columns differ structurally from already-snapshotted datasets. Skip if obs/var schema matches an existing snapshot — registry params + existing snapshot suffice.
6. **Run the builder** with registry-derived params and verify smoke output.
