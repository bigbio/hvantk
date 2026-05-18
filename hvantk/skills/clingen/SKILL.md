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

This skill covers BUILD and UPDATE of the ClinGen Gene-Disease Validity Hail Table. It does NOT cover download (see `hvantk/skills/clingen/cli.py`) or downstream streamer logic (see `hvantk/data/clingen_streamer.py`, which is a cross-cutting utility shared with the GenCC/gene-disease streamer family and intentionally stays outside this plugin folder).

## 2. Source identity

Provider metadata (URL, license, citation) lives in `hvantk/resources/catalog.yaml` (entry `clingen`). The URL/version constants live in `hvantk/core/constants.py` (`CLINGEN_BASE_URL`, `CLINGEN_DOWNLOADS_URL`, `CLINGEN_FILE_PREFIX`, `CLINGEN_HEADER_SKIP_LINES`). Read those files; do not restate.

Stable provider notes the catalog will not capture:
- ClinGen serves a single rolling Gene-Disease Validity CSV; there are no dated archives. Freshness is determined by the file's HTTP `Last-Modified` header (when present) and by the values in the leading metadata block.
- The CSV is generated on demand from the live curation database, so two downloads minutes apart may differ.

## 3. Backend choice + reasoning

`hail`. ClinGen is small (~5k rows) and `_conventions` § 3 allows pandas for mapping tables, but every current consumer (`ClinGenStreamer`, gene-disease join points in PSROC, recipe-driven multi-table joins) reads it as a Hail Table — producing one avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: comma-separated, double-quoted fields, **6-line metadata header** before the column-header line (`CLINGEN_HEADER_SKIP_LINES = 6`). The column header begins with `"GENE SYMBOL"`. Separator rows containing `++++++` are interleaved with the metadata block.
- Preprocessing: the builder streams the file via `hl.hadoop_open` and writes a cleaned temp CSV containing only the column header + data rows. This is required because `hl.import_table` cannot skip arbitrary leading metadata.
- Import: `hl.import_table(delimiter=",", quote='"', impute=False, min_partitions=10)`. All fields stay as strings; no type inference.
- Field renaming is driven by `CLINGEN_GENE_DISEASE_FIELDS` (`hvantk/core/constants.py`). Notable renames: `GENE SYMBOL → gene_symbol`, `GENE ID (HGNC) → hgnc_id`, `DISEASE LABEL → disease_label`, `DISEASE ID (MONDO) → mondo_id`, `MOI → mode_of_inheritance`, `CLASSIFICATION → classification`, `GCEP → gene_curation_expert_panel`.
- ID prefix stripping: the builder strips the `HGNC:` and `MONDO:` prefixes from `hgnc_id` and `mondo_id`. This is asymmetric with the HGNC table (which keeps the prefix); downstream joins (e.g., `clingen_streamer`) account for this.
- Classification levels (`CLINGEN_CLASSIFICATION_LEVELS`, strongest first): `Definitive`, `Strong`, `Moderate`, `Limited`, `Disputed`, `Refuted`. The builder annotates `classification_level` as the numeric position (lower is stronger). Unknown values get `len(CLINGEN_CLASSIFICATION_LEVELS)` and are clamped to the last valid index in the gene-aggregation branch.
- Keying: default `key_by="gene_disease"` keys by `(hgnc_id, mondo_id)`. `key_by="gene"` aggregates per `hgnc_id` (collects diseases, classifications, modes-of-inheritance into sets; counts diseases). The gene branch also derives `max_classification_label` from `max_classification_level` via a clamped lookup.
- TSV export is forwarded to `_create_table_base` (no nested struct to flatten); ClinGen is not special in this regard.

## 5. Output contract

Hail Table at `<output_path>.ht`. Default keying: `(hgnc_id, mondo_id)`. With `key_by="gene"`, keyed by `hgnc_id` only.

Default row schema includes: `hgnc_id`, `gene_symbol`, `disease_label`, `mondo_id`, `mode_of_inheritance`, `sop`, `classification`, `classification_level`, `classification_date`, `gene_curation_expert_panel`, `report_url`. Gene-aggregation rows replace the disease/classification scalars with `disease_labels`, `disease_mondo_pairs`, `mondo_ids`, `classifications`, `modes_of_inheritance`, `max_classification_level`, `max_classification_label`, `n_diseases`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/clingen/plugin.yaml` (drives loader registration; compound dataset key `clingen:gene-disease`).
- Builder: `create_clingen_gene_disease_tb` in `hvantk/skills/clingen/builder.py` (uses `_create_table_base()` per `_conventions` § 4).
- Downloader CLI: `download_cmd` (Click `clingen-download`) in `hvantk/skills/clingen/cli.py`; lifecycle entry-point `download_dataset(raw_dir=...)`. Wired into the umbrella `hvantk download clingen` group in `hvantk/commands/download_cli.py`.
- Dataset class: `ClinGenGeneDiseaseDataset` in `hvantk/skills/clingen/shared/datasets.py`.
- Build CLI: `hvantk mktable clingen-gene-disease` in `hvantk/commands/make_table_cli.py`.
- Streamer (out-of-plugin, intentionally): `hvantk/data/clingen_streamer.py` (`ClinGenStreamer`, subclass of `GeneDiseaseValidityStreamer`); shared with the GenCC/gene-disease streamer family.
- Constants: `CLINGEN_BASE_URL`, `CLINGEN_DOWNLOADS_URL`, `CLINGEN_FILE_PREFIX`, `CLINGEN_HEADER_SKIP_LINES`, `CLINGEN_GENE_DISEASE_FIELDS`, `CLINGEN_CLASSIFICATION_LEVELS` in `hvantk/core/constants.py`.
- Tests: `hvantk/skills/clingen/tests/test_builder.py`, `test_downloader.py`, `test_drift_probe.py`. Streamer tests stay at `hvantk/tests/test_clingen_streamer.py` (test the unmoved streamer).

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Confirm the raw CSV is present at `<raw_dir>/Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv`. If absent, run `hvantk download clingen --output-dir <raw_dir>`.
3. Decide keying: `gene_disease` for full granularity (joins on `(hgnc_id, mondo_id)`); `gene` for gene-level aggregation (downstream geneset extraction).
4. Build via Python (`create_clingen_gene_disease_tb(input_path, output_path, key_by=..., min_classification=..., overwrite=True)`) or CLI (`hvantk mktable clingen-gene-disease --raw-input ... --output-ht ...`).
5. Sanity-check the output: row count plausible (~5k associations live, 11 in fixture); key fields present; HGNC/MONDO prefixes stripped on at least one known row (e.g., BRCA1 → `hgnc_id == "1100"`, not `"HGNC:1100"`).
6. Run validation: `pytest hvantk/skills/clingen/tests -m hail`.

## 8. Update playbook

When ClinGen publishes an updated snapshot (any download is effectively a new snapshot):

1. Re-download: `hvantk download clingen --output-dir <raw_dir> --overwrite`. Capture the `Last-Modified` header (or the snapshot date label) in the PR description.
2. Diff the new CSV header against the previous fixture header. New columns alone are non-breaking — they will not appear in the built table unless added to `CLINGEN_GENE_DISEASE_FIELDS`. Removed/renamed columns require updating `CLINGEN_GENE_DISEASE_FIELDS`.
3. Drift-probe diff: regenerate `hvantk/skills/clingen/tests/drift_fingerprint.json` and inspect for column-list changes or hash changes.
4. If the fixture (`hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv`) is no longer representative (new classification value, new GCEP under test), refresh it from a curated sub-sample.
5. Re-run the snapshot round-trip with `--regenerate-snapshots`; expected diffs: new optional columns, refreshed `classification_date` values. Unexpected diffs (changed key semantics, stripped prefixes leaking back in) require investigation before commit.
6. Open PR; reviewer checks the snapshot diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv`
- `drift_fingerprint`: `hvantk/skills/clingen/tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/clingen/tests -m hail`

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run in a hail-enabled environment, use
> `pytest hvantk/skills/clingen/tests/test_builder.py --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
