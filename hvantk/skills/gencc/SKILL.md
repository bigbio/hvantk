---
name: hvantk:resource-gencc
description: Onboard, build, or update the GenCC submissions resource for hvantk
status: provisional
backend: hail
domain: genomics
---

# GenCC (Gene Curation Coalition) submissions resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

This skill covers BUILD and UPDATE of the GenCC Submissions Hail Table. It does NOT cover download (see `hvantk/skills/gencc/cli.py`) or downstream streamer logic (see `hvantk/data/gencc_streamer.py`, which is a cross-cutting utility shared with the ClinGen/gene-disease streamer family and intentionally stays outside this plugin folder).

## 2. Source identity

Provider metadata (URL, license, citation) lives in `hvantk/resources/catalog.yaml` (entry `gencc`). The URL/version constants live in `hvantk/core/constants.py` (`GENCC_BASE_URL`, `GENCC_FILE_PREFIX`). Read those files; do not restate.

Stable provider notes the catalog will not capture:
- GenCC serves a single rolling submissions TSV; there are no dated archives. Freshness is determined by the file's HTTP `Last-Modified` header (when present) and by the `submitted_as_date` values in the body.
- The TSV is generated on demand by the GenCC backend, so two downloads minutes apart may differ. The body may contain multiline quoted fields and blank lines that `hl.import_table` cannot handle natively, so the dataset's `download()` runs `sanitize_tsv()` in-place before the file is consumed.

## 3. Backend choice + reasoning

`hail`. GenCC is small (~5-10k rows) and `_conventions` § 3 allows pandas for mapping tables, but every current consumer (`GenCCStreamer`, gene-disease join points in PSROC, recipe-driven multi-table joins) reads it as a Hail Table — producing one avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: tab-separated, optional double-quoted fields. The column header is the **first** non-blank line and begins with `sgc_id`. There is no metadata preamble (unlike ClinGen).
- Preprocessing: handled inside the downloader (`sanitize_tsv` strips multiline quoted fields and blank lines). The builder itself calls `hl.import_table` directly on the sanitized TSV.
- Import: `hl.import_table(delimiter="\t", impute=False, min_partitions=10)`. All fields stay as strings; no type inference.
- Field renaming is driven by `GENCC_SUBMISSION_FIELDS` (`hvantk/core/constants.py`). Notable renames: `gene_curie → hgnc_id`, `gene_symbol → gene_symbol`, `disease_curie → mondo_id`, `disease_title → disease_label`, `classification_title → classification`, `moi_title → mode_of_inheritance`, `submitter_title → submitter`, `submitted_as_date → submission_date`, `submitted_as_public_report_url → report_url`, `submitted_as_pmids → pmids`.
- ID prefix stripping: the builder strips the `HGNC:` and `MONDO:` prefixes from `hgnc_id` and `mondo_id` (asymmetric with HGNC, symmetric with ClinGen).
- Classification levels (`GENCC_CLASSIFICATION_LEVELS`, strongest first): `Definitive`, `Strong`, `Moderate`, `Supportive`, `Limited`, `Disputed Evidence`, `Refuted Evidence`, `No Known Disease Relationship`. The builder annotates `classification_level` as the numeric position (lower is stronger); unknown values get `len(GENCC_CLASSIFICATION_LEVELS)` and are clamped to the last valid index in the aggregation branches.
- Keying: default `key_by="gene_disease_submitter"` keys by `(hgnc_id, mondo_id, submitter)` (one row per submitter assertion). `key_by="gene_disease"` aggregates across submitters and re-derives `classification`/`classification_level` from `max_classification_level`. `key_by="gene"` aggregates all diseases per gene.
- TSV export is forwarded to `_create_table_base` (no nested struct to flatten).

## 5. Output contract

Hail Table at `<output_path>.ht`. Default keying: `(hgnc_id, mondo_id, submitter)`.

Default row schema includes: `sgc_id`, `hgnc_id`, `gene_symbol`, `mondo_id`, `disease_label`, `disease_original_id`, `disease_original_label`, `classification`, `classification_level`, `mode_of_inheritance`, `submitter`, `submission_date`, `report_url`, `pmids`. Gene-disease aggregation rows replace per-submitter scalars with `submitters`, `classifications`, `modes_of_inheritance`, `max_classification_level`, `max_classification_label`, `n_submitters`. Gene aggregation additionally collapses `disease_labels`, `disease_mondo_pairs`, `mondo_ids`, `n_diseases`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/gencc/plugin.yaml` (drives loader registration; compound dataset key `gencc:submissions`).
- Builder: `create_gencc_submissions_tb` in `hvantk/skills/gencc/builder.py` (uses `_create_table_base()` per `_conventions` § 4).
- Downloader CLI: `download_cmd` (Click `gencc-download`) in `hvantk/skills/gencc/cli.py`; lifecycle entry-point `download_dataset(raw_dir=...)`. Wired into the umbrella `hvantk download gencc` group in `hvantk/tools/plugins/download_cli.py`.
- Dataset class: `GenCCSubmissionsDataset` in `hvantk/skills/gencc/shared/datasets.py`.
- Build CLI: `hvantk mktable gencc-submissions` in `hvantk/tools/build/make_table_cli.py`.
- Streamer (out-of-plugin, intentionally): `hvantk/data/gencc_streamer.py` (`GenCCStreamer`, subclass of `GeneDiseaseValidityStreamer`); shared with the ClinGen/gene-disease streamer family. Re-exported through `hvantk/data/__init__.py`.
- Constants: `GENCC_BASE_URL`, `GENCC_FILE_PREFIX`, `GENCC_SUBMISSION_FIELDS`, `GENCC_CLASSIFICATION_LEVELS` in `hvantk/core/constants.py`.
- Tests: `hvantk/skills/gencc/tests/test_gencc.py` (dataset class + builder-driven streamer tests), `test_drift_probe.py`.

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Confirm the raw TSV is present at `<raw_dir>/gencc-submissions-<YYYY-MM-DD>.tsv`. If absent, run `hvantk download gencc --output-dir <raw_dir>`.
3. Decide keying: `gene_disease_submitter` for full granularity (per-submitter assertions); `gene_disease` for consensus joins; `gene` for gene-level rollups.
4. Build via Python (`create_gencc_submissions_tb(input_path, output_path, key_by=..., min_classification=..., overwrite=True)`) or CLI (`hvantk mktable gencc-submissions --raw-input ... --output-ht ...`).
5. Sanity-check the output: row count plausible (per-submitter granularity gives more rows than ClinGen); key fields present; HGNC/MONDO prefixes stripped on at least one known row (e.g., BRCA1 → `hgnc_id == "1100"`, not `"HGNC:1100"`).
6. Run validation: `pytest hvantk/skills/gencc/tests -m hail`.

## 8. Update playbook

When GenCC publishes an updated snapshot (any download is effectively a new snapshot):

1. Re-download: `hvantk download gencc --output-dir <raw_dir> --overwrite`. Capture the `Last-Modified` header (or the snapshot date label) in the PR description.
2. Diff the new TSV header against the previous fixture header. New columns alone are non-breaking — they will not appear in the built table unless added to `GENCC_SUBMISSION_FIELDS`. Removed/renamed columns require updating `GENCC_SUBMISSION_FIELDS`.
3. Drift-probe diff: regenerate `hvantk/skills/gencc/tests/drift_fingerprint.json` and inspect for column-list changes or hash changes.
4. If the fixture (`hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv`) is no longer representative (new classification value, new submitter under test), refresh it from a curated sub-sample.
5. Re-run the snapshot round-trip with `--regenerate-snapshots`; expected diffs: new optional columns, refreshed `submission_date` values. Unexpected diffs (changed key semantics, stripped prefixes leaking back in) require investigation before commit.
6. Open PR; reviewer checks the snapshot diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv`
- `drift_fingerprint`: `hvantk/skills/gencc/tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/gencc/tests -m hail`

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run in a hail-enabled environment, use
> `pytest hvantk/skills/gencc/tests/test_builder.py --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
