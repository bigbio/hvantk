---
name: hvantk:conventions
description: Conventions every hvantk resource skill assumes. Read first.
status: provisional
---

# hvantk Resource Conventions

These conventions apply to every per-resource skill. Per-resource skills MAY assume everything below without restating it.

## 1. Repository map

- `hvantk/tables/` — Hail Table builders, including `_create_table_base()` (the boilerplate helper)
- `hvantk/tables/matrix_builders.py` — MatrixTable / anndata builders
- `hvantk/tables/registry.py` — recipe-system registration via `create_table_adapter()`
- `hvantk/commands/` — Click CLI entry points; `make_table_cli.py` and `make_matrix_cli.py` are the dispatch hubs
- `hvantk/datasets/` — provider-specific dataset classes (download + versioning)
- `hvantk/resources/catalog.yaml` — provider catalog (URLs, version cadence, license)
- `hvantk/resources/registry/<domain>/datasets.json` — per-domain JSON registry with one entry per dataset (find your source by `accession` or `title`)
- `hvantk/core/hail_context.py` — Hail initialization (idempotent, thread-safe)
- `hvantk/tests/testdata/<source>/` — fixtures
- `hvantk/tests/snapshots/<source>/` — schema + sample rows snapshots
- `hvantk/skills/<source>/SKILL.md` — per-resource skill

When in doubt, READ existing code under these paths before inferring shape.

## 2. Authoritative spec sources

`resources/catalog.yaml` and `resources/registry/<domain>/datasets.json` are the source of truth for provider metadata: URLs, version strings, license, citation, release cadence. NEVER restate this content in a skill. Reference the catalog instead.

## 3. Keying conventions per data domain

- Variants → key `(locus, alleles)`, Hail Table or MatrixTable
- Genes → key `gene_id`, Hail Table
- Proteins → key `interval` or `protein_id`
- Expression (dense, Hail-friendly) → MatrixTable, rows = genes, cols = samples/cells
- Expression (sparse, single-cell) → anndata h5ad
- Lookup / mapping (e.g., HGNC) → small Hail Table or pandas DataFrame, depending on use

## 4. Required helpers

- `_create_table_base()` — defined in `hvantk/tables/table_builders.py`. Use this for variant/gene Table builders to avoid boilerplate (import, transform, checkpoint, optional TSV export).
- `create_table_adapter()` — `hvantk/tables/registry.py`. Use to register a Hail Table builder in the recipe system via introspection-based parameter mapping.
- `create_matrix_adapter()` — `hvantk/tables/registry.py`. Same role for MatrixTable / anndata builders that take multi-input shapes (e.g., expression matrix + metadata).
- `init_hail()` — `hvantk/core/hail_context.py`. Use to ensure Hail is initialized once. Tests use the session-scoped `hail_session` fixture instead.
- AnnData helpers — `hvantk/core/anndata_utils.py` exposes `build_anndata_metadata`, `save_anndata`, `coerce_obs_for_h5ad`, `annotate_column_summary_ad`. Use for anndata-backed builders.

NEVER paste these helpers' source into a skill. Reference them by path.

## 5. Builder pattern

- Function naming: `create_<source>_tb` (Table) or `build_<source>_ad` (anndata).
- Signature shape: `input_path: str, output_path: str, **kwargs`. Common kwargs: `overwrite: bool`, `export_tsv: bool`, `reference_genome: str`.
- Idempotent: must support `overwrite=True`. Output is checkpointed to disk.
- Returns the built object (`hl.Table`, `hl.MatrixTable`, or `anndata.AnnData`).

## 6. Registry registration

A builder is registered in `hvantk/tables/registry.py` only if it is intended for batch / recipe use. Adapter pattern:

```python
TABLE_BUILDERS["<source>"] = create_table_adapter(
    "hvantk.tables.table_builders", "create_<source>_tb"
)
# Anndata / multi-input matrix builders use the parallel registry:
MATRIX_BUILDERS["<source>"] = create_matrix_adapter(
    "hvantk.tables.matrix_builders", "build_<source>_ad"
)
```

Skip registration for one-off builders.

## 7. CLI command pattern

New commands live under `hvantk/commands/make_table_cli.py` (or `make_matrix_cli.py`). Use Click decorators consistent with existing commands (`@_raw_input_opt`, `@_output_ht_opt`). Naming: `mktable_<source>` for the function, `<source>` for the command name.

## 8. Test pattern

- Round-trip test file: `hvantk/tests/test_<source>_builder.py`.
- Mark with `@pytest.mark.hail` if Hail is required.
- Use the `hail_session` fixture (session-scoped, auto-applied via `conftest.py`).
- Fixtures: `hvantk/tests/testdata/<source>/`. Snapshots: `hvantk/tests/snapshots/<source>/`.
- Assert against snapshots using `hvantk.tests._snapshot_utils`. The `--regenerate-snapshots` pytest flag rewrites snapshots in place.

## 9. Validation contract

Every per-resource SKILL.md MUST declare, by path:

- `fixture`: input file used by the round-trip test
- `schema_snapshot`: `hvantk/tests/snapshots/<source>/schema.json`
- `row_snapshot`: `hvantk/tests/snapshots/<source>/sample_rows.json`
- `test_command`: `pytest hvantk/tests/test_<source>_builder.py -m hail`

## 10. Hard guardrails

- NEVER invent Hail field names. Read the schema from a real run.
- NEVER invent VCF/TSV column names. Read the file header first.
- NEVER assume catalog content. Read `resources/catalog.yaml`.
- NEVER paste code from a builder into the skill. Reference the file path.
- When uncertain, READ existing code (cite which file).

## 11. Out of scope for any skill

- Downloaders. URL drift / version-string handling lives in `hvantk/commands/<source>_downloader.py` and `hvantk/datasets/<source>_datasets.py`.
- Hail context init. Tests use `hail_session`; runtime uses `init_hail()`.
- Cross-resource utilities (gene-ID mapping, locus normalization). Those live in `hvantk/utils/`.
- "How to use the product" — analytical guidance is downstream.
