---
name: hvantk:resource-cptac-phospho
description: Download CPTAC phosphoproteomics data via the cptac Python package and build a sites x samples AnnData plus a wide intermediate TSV for the hvantk PTM pipeline.
status: provisional
backend: anndata
domain: proteomics
---

# CPTAC phosphoproteomics resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. Downloader + dataset class + parser helpers + builder + drift-probe placeholder + parser unit tests live under the plugin folder. Builder round-trip snapshots are NOT yet seeded; the fixture directory exists but is empty. Plugin tests (parser + drift-probe sanity) live in `hvantk/skills/cptac/phospho/tests/`.
- **In scope:** per-cancer-type phospho fetches (one of `CPTAC_CANCER_TYPES`) via the `cptac` Python package, written as a 13-column intermediate TSV consumed by the PTM pipeline (`hvantk/algorithms/ptm/pipeline.py`) AND a wide-format matrix CSV + metadata CSV consumed by the AnnData builder. The builder produces a `(samples x sites)` AnnData with parsed site annotations in `var`.
- **Out of scope:** non-phospho PTM atlases on CPTAC (those would land as sibling datasets under `hvantk/skills/cptac/<ptm-type>/`); pan-cancer joint analysis (the downloader merges per-cancer TSVs but does not produce a joint AnnData).

## 2. Source identity

- **Provider:** Clinical Proteomic Tumor Analysis Consortium (CPTAC) (<https://proteomics.cancer.gov/programs/cptac>), accessed via the [`cptac` Python package](https://pypi.org/project/cptac/) (NOT direct HTTP).
- **Catalog entry:** TODO. No CPTAC entry exists yet in any plugin's `catalog/datasets.json` (verify with `hvantk catalog search cptac`). Cancer-type enum + class map live in `hvantk/skills/cptac/shared/constants.py` (`CPTAC_CANCER_TYPES`, `CPTAC_CANCER_CLASS_MAP`). The installed `cptac` package version is the de-facto release pin (fingerprinted by the drift probe).

## 3. Backend choice + reasoning

**`backend: anndata`, `domain: proteomics`.** Per `_conventions` § 3, AnnData is the right fit for sparse-ish sample x feature matrices. CPTAC phospho is `samples x sites` with parsed site annotations (`gene_symbol`, `amino_acid`, `residue_pos`) in `var` and clinical metadata in `obs`. The intermediate site-level TSV stays alongside the AnnData because the downstream PTM pipeline (`hvantk/algorithms/ptm/pipeline.py`) reads the TSV directly -- the AnnData is for analysis, the TSV is for the PTM-site joining workflow.

## 4. Raw format & gotchas

The `cptac` package returns a `pd.DataFrame` with a `MultiIndex` on columns whose first two levels are `Gene` and `Site`. `Site` tokens come in two shapes:

- Single site: `S65`, `T185`, `Y243` (one of S/T/Y followed by digits).
- Multi-site: `T185_Y187` (multiple sites separated by underscores). The parser expands each into its own row.

Aggregation gotchas (see `extract_phospho_sites` and `write_matrix_csv` in `shared/datasets.py`):

- The same `(gene, amino_acid, position)` can appear via multiple peptide columns -- `extract_phospho_sites` aggregates `n_observations` and weighted-averages `mean_intensity` across columns. `write_matrix_csv` groups duplicate site labels by mean intensity per sample.
- Sample-level tissue type is inferred from the sample-id suffix: `*.N` -> `normal`, everything else -> `tumor`. This is then joined onto the clinical metadata as a new `tissue_type` column.
- Some cancer types have no normal samples; the downloader logs and continues without a `*-normal.tsv`.
- The package fetches data from upstream **at runtime** (no static URLs). The `cptac:phospho` drift probe therefore fingerprints the installed package version, not an HTTP `Last-Modified` header.

For the AnnData builder, the input is the *wide* matrix CSV (sites in rows / first column, samples in remaining columns). Site IDs follow the `<gene>_<aa><pos>` convention (e.g. `TP53_S15`). The builder parses these via `_parse_site_id` and writes `gene_symbol`, `amino_acid`, `residue_pos` into `var`.

## 5. Output contract

Per-cancer download (in `output_dir/` or `raw_dir/`):

- `cptac-phospho-<ct>.tsv` -- combined site table (tumor + normal), 13-column.
- `cptac-phospho-<ct>-tumor.tsv` -- tumor-only sites.
- `cptac-phospho-<ct>-normal.tsv` -- normal-only sites (if any).
- `cptac-phospho-<ct>-matrix.csv` -- wide intensity matrix (sites x samples).
- `cptac-phospho-<ct>-metadata.csv` -- clinical metadata with `tissue_type` + `cancer_type`.

TSV column order (`_TSV_COLUMNS`): `accession`, `gene_symbol`, `position`, `description`, `amino_acid`, `ensembl_xrefs`, `sequence_length`, `n_observations`, `source_db`, `evidence_type`, `tissue_type`, `cancer_type`, `mean_intensity`. `source_db` is always `"CPTAC"`; `evidence_type` is always `"mass_spectrometry"`.

AnnData (`.h5ad`) builder output:

- **Shape:** `(n_samples, n_sites)`.
- **`X`:** `float32`, missing intensities as `NaN`.
- **`obs`:** indexed by sample id; columns mirror the metadata CSV.
- **`var`:** indexed by site id (e.g. `TP53_S15`); columns `gene_symbol`, `amino_acid`, `residue_pos`.
- **`uns["column_summary"]`** as for `cptac:expression`.
- **Provenance:** stamped via `ctx.provenance(schema_id="cptac-phospho-v1")` and persisted in the platform sidecar `.provenance.json`.

## 6. hvantk integration points

- **Dataset class + parser helpers:** `CPTACPhosphoDataset`, `parse_phospho_site`, `extract_phospho_sites`, `write_intermediate_tsv`, `write_matrix_csv`, `write_metadata_csv`, `parse_raw_dir` in `hvantk/skills/cptac/shared/datasets.py`.
- **Builder:** `build_cptac_phospho(parsed_input, ctx, *, cancer_type=None, site_id_col="Site", sample_id_col="SampleID", **_)` in `hvantk/skills/cptac/phospho/builder.py`. `parsed_input` is either a `{"expression": <matrix.csv>, "metadata": <metadata.csv>}` dict (direct API) OR a raw download directory, in which case `cancer_type` selects the per-cancer `cptac-phospho-<ct>-matrix.csv` / `-metadata.csv` pair (this is the `reprocess` path — see § 7). `site_id_col` defaults to `"Site"` to match the matrix CSV header. Returns an `ExpressionMatrix` (shares `create_anndata_from_cptac_phospho` in `hvantk/skills/cptac/shared/cptac.py`).
- **Downloader CLI:** `download_cmd` in `hvantk/skills/cptac/phospho/cli.py`, declared as command `cptac-phospho-download` in the manifest's `cli:` block. The plugin loader strips the `-download` suffix and wires it as `hvantk download cptac-phospho` (via `hvantk/tools/plugins/download_cli.py`).
- **Lifecycle entry points:** `download_dataset` (in `phospho/cli.py`), wired via `lifecycle.download` in `plugin.yaml`. There is **no** `lifecycle.parse` for the AnnData build: the download already writes the builder's matrix+metadata inputs, so `reprocess` hands the raw dir straight to the builder (mirrors `cptac:expression`). `parse_raw_dir` (in `shared/datasets.py`) remains a public helper that consolidates the site-level TSVs for the PTM pipeline; it is not wired into the AnnData lifecycle.
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/cptac/phospho/drift_probe.py` (fingerprints the installed `cptac` package version).
- **Downstream consumer:** `hvantk/algorithms/ptm/pipeline.py` (`PTMBuildConfig.cptac_tsv`) -- reads the per-cancer intermediate TSV and maps PTM sites to genomic coordinates. Exposed at the user-facing level by `hvantk/algorithms/ptm/atlas.py` (`PTMAtlasConfig.cptac_tsv`).
- **Plugin manifest:** `hvantk/skills/cptac/plugin.yaml` (compound dataset key `cptac:phospho`). The loader auto-resolves the dataset via `get_registry().get_dataset("cptac:phospho")`; top-level builds run through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).
- **Tests:** `hvantk/skills/cptac/phospho/tests/` (parser unit tests in `test_phospho.py` + drift-probe sanity in `test_drift_probe.py`).

## 7. Workflow steps

When invoked to build or update CPTAC phospho data:

**One AnnData per cancer type.** CPTAC phospho matrices differ in sites and samples
across cancer types, so a single `reprocess` builds one cancer's AnnData; run it once
per cancer type. End-to-end via the documented CLI:

```bash
hvantk reprocess cptac:phospho \
  --raw-dir /data/cptac \
  --output cptac_phospho_brca.h5ad \
  --plugin-arg cancer_type=brca
```

The one `cancer_type` flows to both the download (fetch only that cancer) and the
build (select its matrix+metadata pair). Steps:

1. **Download.** Via the download CLI (`hvantk download cptac-phospho --cancer-type brca -o /data/cptac`) or the lifecycle entry point (`reprocess` calls `download_dataset(raw_dir=..., cancer_type=brca)`). Pass `--all` to fetch every cancer type; the downloader then also produces a pan-cancer merge TSV and skips-and-continues past any cancer that fails upstream.
2. **Builder.** `build_cptac_phospho(<raw_dir>, ctx, cancer_type="brca")` (or a `{"expression","metadata"}` dict for the direct API) produces the `(samples x sites)` `ExpressionMatrix`.
3. **Validate.** `pytest hvantk/skills/cptac/phospho/tests` (parser + builder + drift-probe).

## 8. Update playbook

CPTAC publishes refreshed datasets through the `cptac` Python package. When a new package release lands:

1. Run the drift probe: `python -c "from hvantk.skills.cptac.phospho.drift_probe import fetch_fingerprint; print(fetch_fingerprint())"`. A change in `source_version` (the installed `cptac` package version) is the trigger.
2. Update the Poetry dep pin in `pyproject.toml` (`cptac = { version = "*", optional = true }` -- pin to the new release if reproducibility matters).
3. Re-download every cancer type of interest with `--overwrite` and spot-check row counts and `tissue_type` distributions against the previous build.
4. If any new modification notation or column shape appears, document it in § 4.
5. Re-regenerate `tests/drift_fingerprint.json` to record the new package version.
6. If a builder snapshot exists (TODO), run `pytest hvantk/skills/cptac/phospho/tests --regenerate-snapshots` and inspect the diff.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/cptac/phospho/tests/testdata/raw/cptac-phospho/` (directory present, not yet seeded -- parser unit tests build mock multi-index DataFrames on the fly).
- **schema_snapshot:** `hvantk/skills/cptac/phospho/tests/snapshots/schema.json` (TODO -- created on first `--regenerate-snapshots` run).
- **row_snapshot:** `hvantk/skills/cptac/phospho/tests/snapshots/sample_rows.json` (TODO -- same).
- **test_command:** `pytest hvantk/skills/cptac/phospho/tests`.
- **drift_fingerprint:** `hvantk/skills/cptac/phospho/tests/drift_fingerprint.json` (placeholder shape; refresh via the update playbook).

The plugin manifest already declares these paths so the loader contract holds. Parser unit tests + drift-probe sanity test pass today; the builder round-trip snapshot is the gap to close in a follow-up.

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run, use
> `pytest hvantk/skills/cptac/phospho/tests --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
