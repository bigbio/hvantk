---
name: hvantk:resource-cptac-expression
description: Build an AnnData (samples x genes) from long-format CPTAC protein-expression matrices plus sample metadata.
status: provisional
backend: anndata
domain: proteomics
---

# CPTAC protein expression resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

- **Status:** provisional. Downloader + dataset class + builder + drift-probe placeholder live under the plugin folder. Builder round-trip snapshots and the fixture are seeded; `test_builder.py` asserts against both. Plugin tests (downloader unit tests + drift-probe sanity) live in `hvantk/skills/cptac/expression/tests/`.
- **In scope:** per-cancer-type whole-proteome fetches (one of `CPTAC_CANCER_TYPES`) via the `cptac` Python package, written as a long-format matrix CSV (`SampleID`, `GeneID`, `Gene Name`, `Expression`) plus a clinical metadata CSV. The builder produces a `(samples x genes)` AnnData from that pair, one cancer type at a time. A long-format TSV/CSV staged by hand (the pre-#386 workflow) still works via the `{"expression","metadata"}` mapping form of `parsed_input`.
- **Out of scope:** phosphoproteomics (sibling dataset `cptac:phospho`); pan-cancer joint analysis (`--all` downloads every cancer type but does not merge them -- gene sets and sample columns differ per cancer type, unlike phospho's site-level TSV which does have a pan-cancer merge); any Hail-Table or wide-format representation.

## 2. Source identity

- **Provider:** Clinical Proteomic Tumor Analysis Consortium (CPTAC) (<https://proteomics.cancer.gov/programs/cptac>).
- **Catalog entry:** TODO. No CPTAC entry exists yet in any plugin's `catalog/datasets.json` (verify with `hvantk catalog search cptac`). The cancer-type enum + class map live in `hvantk/skills/cptac/shared/constants.py` (`CPTAC_CANCER_TYPES`, `CPTAC_CANCER_CLASS_MAP`). When a dedicated catalog entry lands, point this plugin's `source.catalog_ref` at it and remove this TODO.

## 3. Backend choice + reasoning

**`backend: anndata`, `domain: proteomics`.** Per `_conventions` § 3, AnnData is the right fit for sparse-ish sample x feature matrices keyed by sample identifier. CPTAC protein expression is `samples x genes`, fits in memory for any single cancer cohort, and pairs naturally with sample-level clinical metadata in `obs`. A Hail Table would force a Spark session for a sub-100 MB matrix and lose the `obs` ergonomics.

## 4. Raw format & gotchas

The `cptac` package's `get_proteomics()` returns a `pd.DataFrame` (samples x genes) with a 2-level `MultiIndex` on columns: `Name` (gene symbol, e.g. `TP53`) and `Database_ID` (an Ensembl protein accession, e.g. `ENSP00000000233.5` -- confirmed live, see below). `Database_ID` is what keeps columns unique when a protein group maps several accessions to one symbol (isoforms) -- the same shape `get_phosphoproteomics()` documents for the sibling phospho dataset (see its SKILL.md § 4), minus the `Site` level.

`melt_proteomics_wide_to_long` (in `hvantk/skills/cptac/shared/datasets.py`) flattens this to the long format the builder reads:

- `GeneID` = `Database_ID` (falls back to `Name` when a column has no accession, so a blank ID never collapses two different genes together).
- `Gene Name` = `Name`.
- Any residual duplicate `(SampleID, GeneID)` rows are collapsed by averaging `Expression` -- mirroring `write_matrix_csv`'s duplicate-site handling for phospho.
- Sample-level tissue type is inferred from the sample-id suffix (`*.N` -> `normal`, everything else -> `tumor`) by the shared `_tag_tissue_type` helper, then joined onto the clinical metadata as a `tissue_type` column -- identical convention to `cptac:phospho`.
- **Live-verification status (2026-10-12, cptac 1.5.14):** `get_proteomics(source="umich")` is confirmed live for `brca` (125 samples x 12,922 proteins; the exact `(Name, Database_ID)` column shape above) and round-trips through `build_cptac_expression` into a `(125, 12922)` AnnData. The clinical-metadata half (`get_clinical(source="mssm")`) is **not yet verified live** -- Zenodo returned 504 for the whole verification session, so `source="mssm"`'s validity for `get_clinical()` remains an informed guess (mirroring the phospho downloader's choice), not a confirmed one.
- **No default-source retry.** The downloader always passes an explicit `source` (`"umich"` for proteomics, `"mssm"` for clinical) and never retries a getter with no `source` on failure: `cptac` 1.5.14's no-source path hits `len(sources_for_data)` on a generator inside `cptac/cancers/cancer.py` `get_dataframe()` (~line 682) and raises an opaque `TypeError: object of type 'generator' has no len()` instead of a clear "no source available" error. A fetch failure instead raises a `RuntimeError` naming the cancer type, data type (`proteomics`/`clinical`), and source that failed. `cptac:phospho`'s pre-existing `_get_phospho`/`get_clinical` fallbacks (`shared/datasets.py` ~line 481 and ~528) still retry with no source and so still carry this bug; fixing them is out of scope here.

The builder itself then accepts that long-format table (one row per `(sample, gene)` observation) and pivots it to wide:

- Required columns (defaults): `GeneID`, `SampleID`, `Expression`. Override via `gene_id_col`, `sample_id_col`, `expression_col`.
- Optional: `Gene Name` (column name configurable via `gene_name_col`) -- if present, it is propagated to `var`.
- Delimiter is auto-detected by `pd.read_csv(sep=None, engine="python")`.

Metadata file is read with the same auto-detect logic and indexed by `SampleID` (override via `sample_id_col`). Any extra columns become `obs` columns via a left join, so samples in the metadata that are absent from the expression matrix are silently dropped (the join is on `wide.index`).

`X` is coerced to `np.float32`. Missing values from the pivot become `NaN` and are preserved in `X`.

## 5. Output contract

Per-cancer download (in `output_dir/` or `raw_dir/`):

- `cptac-expression-<ct>-matrix.csv` -- long-format table: `SampleID`, `GeneID`, `Gene Name`, `Expression`.
- `cptac-expression-<ct>-metadata.csv` -- clinical metadata with `tissue_type` + `cancer_type` columns (same convention as `cptac:phospho`).

AnnData (`.h5ad`) builder output:

- **Artifact:** an `ExpressionMatrix` (AnnData-backed) returned by the builder; the platform persists it to `.h5ad` at the `--output` path the caller passes to `hvantk reprocess`.
- **Shape:** `(n_samples, n_genes)` -- samples in `obs`, genes in `var`.
- **`X`:** `float32`, missing values as `NaN`.
- **`obs`:** indexed by sample id; columns mirror the metadata file (after `set_index`).
- **`var`:** indexed by gene id (name = `gene_id_col`); contains `gene_name_col` if it was present in the input.
- **`uns["column_summary"]`:** added by `annotate_column_summary_ad`.
- **Provenance:** stamped on the returned `ExpressionMatrix` via `ctx.provenance(schema_id="cptac-expression-v1")`; persisted by the platform as a sidecar `.provenance.json` next to the `.h5ad`.

## 6. hvantk integration points

- **Dataset class + parser helpers:** `CPTACExpressionDataset`, `melt_proteomics_wide_to_long` in `hvantk/skills/cptac/shared/datasets.py`. `resolve_raw_dir_inputs`, `_tag_tissue_type`, and `download_cancer_types` are shared with `cptac:phospho` in the same module.
- **Builder:** `build_cptac_expression(parsed_input, ctx, *, cancer_type=None, gene_id_col="GeneID", gene_name_col="Gene Name", sample_id_col="SampleID", expression_col="Expression", **_)` in `hvantk/skills/cptac/expression/builder.py`. `parsed_input` is either a `{"expression","metadata"}` dict (direct API) OR a raw download directory, in which case `cancer_type` selects the per-cancer `cptac-expression-<ct>-matrix.csv` / `-metadata.csv` pair (this is the `reprocess` path -- see § 7). Returns an `ExpressionMatrix` (shares `create_anndata_from_cptac_long` in `hvantk/skills/cptac/shared/cptac.py`).
- **Downloader CLI:** `download_cmd` in `hvantk/skills/cptac/expression/cli.py`, declared as command `cptac-expression-download` in the manifest's `cli:` block. The plugin loader strips the `-download` suffix and wires it as `hvantk download cptac-expression`.
- **Lifecycle entry points:** `download_dataset` (in `expression/cli.py`), wired via `lifecycle.download` in `plugin.yaml`. There is **no** `lifecycle.parse`: the download already writes the builder's matrix+metadata inputs, so `reprocess` hands the raw dir straight to the builder (mirrors `cptac:phospho`).
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/cptac/expression/drift_probe.py` (fingerprints the installed `cptac` Python package version).
- **Plugin manifest:** `hvantk/skills/cptac/plugin.yaml` (compound dataset key `cptac:expression`). The loader auto-resolves the dataset via `get_registry().get_dataset("cptac:expression")`; top-level builds run through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).
- **Tests:** `hvantk/skills/cptac/expression/tests/` (downloader unit tests -- including the fetch-both-then-atomically-publish guarantee -- in `test_downloader.py` + builder round-trip in `test_builder.py` + drift-probe sanity in `test_drift_probe.py`).

## 7. Workflow steps

When invoked to build or update CPTAC expression data:

**One AnnData per cancer type.** CPTAC expression matrices differ in genes and
samples across cancer types, so a single `reprocess` builds one cancer's AnnData;
run it once per cancer type. End-to-end via the documented CLI:

```bash
hvantk reprocess cptac:expression \
  --raw-dir /data/cptac \
  --output cptac_expression_brca.h5ad \
  --plugin-arg cancer_type=brca
```

The one `cancer_type` flows to both the download (fetch only that cancer) and the
build (select its matrix+metadata pair). Steps:

1. **Download.** Via the download CLI (`hvantk download cptac-expression --cancer-type brca -o /data/cptac`) or the lifecycle entry point (`reprocess` calls `download_dataset(raw_dir=..., cancer_type=brca)`). Pass `--all` to fetch every cancer type (no pan-cancer merge -- see § 1); the downloader skips-and-continues past any cancer that fails upstream.
2. **Builder.** `build_cptac_expression(<raw_dir>, ctx, cancer_type="brca")` (or a `{"expression","metadata"}` dict for the direct API / a hand-staged TSV) produces the `(samples x genes)` `ExpressionMatrix`.
3. **Validate.** `pytest hvantk/skills/cptac/expression/tests` (downloader + builder + drift-probe).

## 8. Update playbook

CPTAC publishes refreshed datasets through the `cptac` Python package. When a new package release lands:

1. Run the drift probe: `python -c "from hvantk.skills.cptac.expression.drift_probe import fetch_fingerprint; print(fetch_fingerprint())"`. A change in `source_version` (the installed `cptac` package version) is the trigger.
2. Re-download every cancer type of interest with `--overwrite` and spot-check sample / gene counts and `tissue_type` distributions against the previous build.
3. If any new column shape appears from `get_proteomics()` (e.g. a different `Database_ID` convention), document it in § 4 and adjust `melt_proteomics_wide_to_long`.
4. Regenerate `tests/drift_fingerprint.json` with `hvantk drift --regenerate cptac:expression` once the package version is updated, and commit the rebuild-ledger row it writes with it.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/cptac/expression/tests/testdata/raw/cptac-expression/` — seeded (`expression.tsv` + `metadata.tsv`).
- **schema_snapshot:** `hvantk/skills/cptac/expression/tests/snapshots/schema.json` — seeded.
- **row_snapshot:** `hvantk/skills/cptac/expression/tests/snapshots/sample_rows.json` — seeded.
- **command:** `pytest hvantk/skills/cptac/expression/tests`.
- **drift_fingerprint:** `hvantk/skills/cptac/expression/tests/drift_fingerprint.json` (placeholder shape; refresh via the update playbook).

The plugin manifest declares these paths so the loader contract holds. Downloader, drift-probe sanity, and builder round-trip tests all pass.

> **Snapshot status:** seeded. `test_builder.py` asserts the build against both
> snapshots. Regenerate after an intentional schema change with
> `pytest hvantk/skills/cptac/expression/tests/test_builder.py --regenerate-snapshots`, then commit.
