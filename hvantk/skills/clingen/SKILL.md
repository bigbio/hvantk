---
name: hvantk:resource-clingen
description: Onboard, build, or update the ClinGen Gene-Disease Validity resource for hvantk
status: provisional
backend: hail
domain: genomics
---

# ClinGen Gene-Disease Validity resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

This skill covers BUILD and UPDATE of the ClinGen Gene-Disease Validity Hail Table. It does NOT cover download (see `hvantk/skills/clingen/cli.py`) or downstream streamer logic. The per-plugin streamer subclass lives in `hvantk/skills/clingen/streamers.py` (`ClinGenGeneDiseaseTableStreamer`), built on the shared base `GeneDiseaseTableStreamer` in `hvantk/core/streamers/gene_disease_table.py`.

## 2. Source identity

Provider metadata (URL, license, citation) lives in the plugin's `catalog/datasets.json` (or query it with `hvantk catalog show <accession>` once a ClinGen catalog entry exists). The URL/version constants live in `hvantk/skills/clingen/shared/constants.py` (`CLINGEN_BASE_URL`, `CLINGEN_DOWNLOADS_URL`, `CLINGEN_FILE_PREFIX`, `CLINGEN_HEADER_SKIP_LINES`). Read those files; do not restate.

Stable provider notes the catalog will not capture:
- ClinGen serves a single rolling Gene-Disease Validity CSV; there are no dated archives. Freshness is determined by the file's HTTP `Last-Modified` header (when present) and by the values in the leading metadata block.
- The CSV is generated on demand from the live curation database, so two downloads minutes apart may differ.

## 3. Backend choice + reasoning

`hail`. ClinGen is small (~5k rows) and `_conventions` § 3 allows pandas for mapping tables, but every current consumer (`ClinGenStreamer`, gene-disease join points in PSROC, recipe-driven multi-table joins) reads it as a Hail Table — producing one avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: comma-separated, double-quoted fields, **6-line metadata header** before the column-header line (`CLINGEN_HEADER_SKIP_LINES = 6`). The column header begins with `"GENE SYMBOL"`. Separator rows containing `++++++` are interleaved with the metadata block.
- Preprocessing: the builder streams the file via `hl.hadoop_open` and writes a cleaned temp CSV containing only the column header + data rows. This is required because `hl.import_table` cannot skip arbitrary leading metadata.
- Import: `hl.import_table(delimiter=",", quote='"', impute=False, min_partitions=10)`. All fields stay as strings; no type inference.
- Field renaming is driven by `CLINGEN_GENE_DISEASE_FIELDS` (`hvantk/skills/clingen/shared/constants.py`). Notable renames: `GENE SYMBOL → gene_symbol`, `GENE ID (HGNC) → hgnc_id`, `DISEASE LABEL → disease_label`, `DISEASE ID (MONDO) → mondo_id`, `MOI → mode_of_inheritance`, `CLASSIFICATION → classification`, `GCEP → gene_curation_expert_panel`.
- ID prefix stripping: the builder strips the `HGNC:` and `MONDO:` prefixes from `hgnc_id` and `mondo_id`. This is asymmetric with the HGNC table (which keeps the prefix); downstream joins (e.g., `clingen_streamer`) account for this.
- Classification levels (`CLINGEN_CLASSIFICATION_LEVELS`, strongest first): `Definitive`, `Strong`, `Moderate`, `Limited`, `Disputed`, `Refuted`. The builder annotates `classification_level` as the numeric position (lower is stronger). Unknown values get `len(CLINGEN_CLASSIFICATION_LEVELS)`.
- Keying: keyed by `(hgnc_id, mondo_id)`.

## 5. Output contract

Hail Table at `<output_path>.ht`, keyed by `(hgnc_id, mondo_id)`.

Row schema includes: `hgnc_id`, `gene_symbol`, `disease_label`, `mondo_id`, `mode_of_inheritance`, `sop`, `classification`, `classification_level`, `classification_date`, `gene_curation_expert_panel`, `report_url`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/clingen/plugin.yaml` (the loader resolves this dataset via `get_registry().get_dataset("clingen:gene-disease")`; top-level builds run through `run_builder_for_spec` in `hvantk/core/plugin/run_builder.py`).
- Builder: `build_clingen_gene_disease` in `hvantk/skills/clingen/builder.py`. Signature `(parsed_input, ctx, *, min_classification=None, fields=None) -> AnnotationTable`. The table is built inline (`hl.import_table` + renames/transforms); the only shared helper is `cleanup_temp_file` from `hvantk/core/utils/hail_helpers.py`.
- Downloader CLI: `download_cmd` (Click `clingen-download`) in `hvantk/skills/clingen/cli.py`; lifecycle entry-point `download_dataset(raw_dir=...)` (declared under `lifecycle.download` in `plugin.yaml`).
- Dataset class: `ClinGenGeneDiseaseDataset` in `hvantk/skills/clingen/shared/datasets.py`.
- Build CLI: `hvantk reprocess clingen:gene-disease --raw-dir <dir> --output <path>.ht` (pass builder kwargs via `--plugin-arg key=value`, e.g. `--plugin-arg min_classification=Strong`).
- Streamer (per-plugin): `ClinGenGeneDiseaseTableStreamer` in `hvantk/skills/clingen/streamers.py`, subclass of `GeneDiseaseTableStreamer` (`hvantk/core/streamers/gene_disease_table.py`).
- Constants: `CLINGEN_BASE_URL`, `CLINGEN_DOWNLOADS_URL`, `CLINGEN_FILE_PREFIX`, `CLINGEN_HEADER_SKIP_LINES`, `CLINGEN_GENE_DISEASE_FIELDS`, `CLINGEN_CLASSIFICATION_LEVELS` in `hvantk/skills/clingen/shared/constants.py`.
- Tests: `hvantk/skills/clingen/tests/test_downloader.py`, `test_drift_probe.py`, `test_streamer.py`.

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Confirm the raw CSV is present at `<raw_dir>/Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv`. If absent, run the plugin downloader (Click `clingen-download`, e.g. `hvantk download clingen --output-dir <raw_dir>`).
3. The builder keys the output by `(hgnc_id, mondo_id)`. Optionally pass `min_classification` to drop weaker associations, and `fields` to select a column subset.
4. Build via CLI: `hvantk reprocess clingen:gene-disease --raw-dir <dir> --output <path>.ht [--plugin-arg min_classification=Strong]`.
5. Sanity-check the output: row count plausible (~5k associations live, 11 in fixture); key fields present; HGNC/MONDO prefixes stripped on at least one known row (e.g., BRCA1 → `hgnc_id == "1100"`, not `"HGNC:1100"`).
6. Run validation: `pytest hvantk/skills/clingen/tests -m hail`.

## 8. Update playbook

When ClinGen publishes an updated snapshot (any download is effectively a new snapshot):

1. Re-download: `hvantk download clingen --output-dir <raw_dir> --overwrite`. Capture the `Last-Modified` header (or the snapshot date label) in the PR description.
2. Diff the new CSV header against the previous fixture header. New columns alone are non-breaking — they will not appear in the built table unless added to `CLINGEN_GENE_DISEASE_FIELDS`. Removed/renamed columns require updating `CLINGEN_GENE_DISEASE_FIELDS`.
3. Drift-probe diff: regenerate `hvantk/skills/clingen/tests/drift_fingerprint.json` and inspect for column-list changes or hash changes.
4. If the fixture (`hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv`) is no longer representative (new classification value, new GCEP under test), refresh it from a curated sub-sample.
5. Re-run the plugin tests (`pytest hvantk/skills/clingen/tests -m hail`); expected diffs: new optional columns, refreshed `classification_date` values. Unexpected diffs (changed key semantics, stripped prefixes leaking back in) require investigation before commit.
6. Open PR; reviewer checks the diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv`
- `drift_fingerprint`: `hvantk/skills/clingen/tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/clingen/tests -m hail`

> **Snapshot status:** the `tests/snapshots/schema.json` and
> `tests/snapshots/sample_rows.json` paths declared in `plugin.yaml` have NOT yet
> been seeded for this plugin (the `tests/snapshots/` directory does not exist
> yet). Until seeded, there is no fixed-schema round-trip check; the live tests
> are `test_downloader.py`, `test_drift_probe.py`, and `test_streamer.py`.
