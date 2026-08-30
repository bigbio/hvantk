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

This skill covers BUILD and UPDATE of the GenCC Submissions Hail Table. It does NOT cover download (see `hvantk/skills/gencc/cli.py`) or downstream streamer logic (see `hvantk/skills/gencc/streamers.py`, `GenCCGeneDiseaseTableStreamer`, a thin subclass of the shared `GeneDiseaseTableStreamer` base in `hvantk/core/streamers/gene_disease_table.py`).

## 2. Source identity

Provider metadata (URL, license, citation) lives in the plugin's `catalog/datasets.json` (or query it with `hvantk catalog show <accession>` once a GenCC catalog entry exists). The URL/version constants live in `hvantk/skills/gencc/shared/constants.py` (`GENCC_BASE_URL`, `GENCC_FILE_PREFIX`). Read those files; do not restate.

Stable provider notes the catalog will not capture:
- GenCC serves a single rolling submissions TSV; there are no dated archives. Freshness is **not** determined by `Last-Modified`: the drift probe's content signal is the HTTP `Content-Length` header (`extras.content_length`) instead. Across the 4 fingerprint commits from 2026-07-27 to 2026-08-23 the checksum (`6f07ac79f9…`) never moved while `Last-Modified` moved on every regeneration. `Last-Modified` is still recorded, under `informational` (excluded from drift comparison) alongside the `submitted_as_date` values in the body; `source_version` is `null`, and the probe (`probe_version` 2) fails closed — raises `DriftProbeError` — if the server omits `Content-Length`.
- The TSV is generated on demand by the GenCC backend, so two downloads minutes apart may differ. The body may contain multiline quoted fields and blank lines that `hl.import_table` cannot handle natively, so the dataset's `download()` runs `sanitize_tsv()` in-place before the file is consumed.

## 3. Backend choice + reasoning

`hail`. GenCC is small (~5-10k rows) and `_conventions` § 3 allows pandas for mapping tables, but every current consumer (`GenCCGeneDiseaseTableStreamer`, gene-disease join points in PSROC, multi-table joins) reads it as a Hail Table — producing one avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: tab-separated, optional double-quoted fields. The column header is the **first** non-blank line and begins with `sgc_id`. There is no metadata preamble (unlike ClinGen).
- Preprocessing: handled inside the downloader (`sanitize_tsv` strips multiline quoted fields and blank lines). The builder itself calls `hl.import_table` directly on the sanitized TSV.
- GenCC rate-limits: a scheduled drift regeneration on 2026-08-04 got back HTTP 429 and failed the run outright. The drift probe uses `request_with_retry` (`hvantk/core/utils/http.py`) instead of bare `requests` for its HEAD and GET calls — it retries transient statuses (429 plus the 5xx family) and connection errors with exponential backoff, honouring `Retry-After` but clamped to a max sleep so a long backoff can't stall CI. Still load-bearing: the 2026-08-29 baseline regeneration for `probe_version: 2` hit 429 again and the retry recovered it on attempt 3/4, so the probe completed rather than failing the run.
- Import: `hl.import_table(delimiter="\t", impute=False, min_partitions=10)`. All fields stay as strings; no type inference.
- Field renaming is driven by `GENCC_SUBMISSION_FIELDS` (`hvantk/skills/gencc/shared/constants.py`). Notable renames: `gene_curie → hgnc_id`, `gene_symbol → gene_symbol`, `disease_curie → mondo_id`, `disease_title → disease_label`, `classification_title → classification`, `moi_title → mode_of_inheritance`, `submitter_title → submitter`, `submitted_as_date → submission_date`, `submitted_as_public_report_url → report_url`, `submitted_as_pmids → pmids`.
- ID prefix stripping: the builder strips the `HGNC:` and `MONDO:` prefixes from `hgnc_id` and `mondo_id` (asymmetric with HGNC, symmetric with ClinGen).
- Classification levels (`GENCC_CLASSIFICATION_LEVELS`, strongest first): `Definitive`, `Strong`, `Moderate`, `Supportive`, `Limited`, `Disputed Evidence`, `Refuted Evidence`, `No Known Disease Relationship`. The builder annotates `classification_level` as the numeric position (lower is stronger); unknown values get `len(GENCC_CLASSIFICATION_LEVELS)`.
- Keying: the builder always keys by `(hgnc_id, mondo_id, submitter)` (one row per submitter assertion). There is no `key_by`/aggregation parameter.
- Builder params (passed as `--plugin-arg key=value`): `min_classification` filters to assertions at or above the given confidence level (must be one of `GENCC_CLASSIFICATION_LEVELS`); `fields` selects a subset of row fields.
- The builder constructs the table inline with `hl.import_table` plus rename/annotate/key_by transforms, then wraps it via `AnnotationTable.from_hail`. There is no `_create_table_base` helper and no `output_path`/`overwrite`/`export_tsv` plumbing in the builder signature; the real signature is `build_gencc_submissions(parsed_input, ctx, *, min_classification=None, fields=None) -> AnnotationTable`.

## 5. Output contract

Hail Table at `<output_path>.ht`, keyed by `(hgnc_id, mondo_id, submitter)`.

Row schema includes: `sgc_id`, `hgnc_id`, `gene_symbol`, `mondo_id`, `disease_label`, `disease_original_id`, `disease_original_label`, `classification`, `classification_level`, `mode_of_inheritance`, `submitter`, `submission_date`, `report_url`, `pmids`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/gencc/plugin.yaml` (drives loader registration; compound dataset key `gencc:submissions`).
- Builder: `build_gencc_submissions` in `hvantk/skills/gencc/builder.py`.
- Downloader CLI: `download_cmd` (Click `gencc-download`) in `hvantk/skills/gencc/cli.py`; lifecycle entry-point `download_dataset(raw_dir=..., overwrite=...)` (declared under `lifecycle.download` in `plugin.yaml`). The umbrella download group lives in `hvantk/tools/plugins/download_cli.py`.
- Dataset class: `GenCCSubmissionsDataset` in `hvantk/skills/gencc/shared/datasets.py`.
- Build CLI: `hvantk reprocess gencc:submissions --raw-dir <dir> --output <path>.ht` (pass builder kwargs via `--plugin-arg key=value`, e.g. `--plugin-arg min_classification=Strong`). The plugin loader (`hvantk/core/plugin/loader.py`) auto-resolves the dataset from `plugin.yaml` via `get_registry().get_dataset("gencc:submissions")`; the build runs through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).
- Streamer: `GenCCGeneDiseaseTableStreamer` in `hvantk/skills/gencc/streamers.py`, a thin subclass of `GeneDiseaseTableStreamer` (base in `hvantk/core/streamers/gene_disease_table.py`); shared with the ClinGen/gene-disease streamer family.
- Constants: `GENCC_BASE_URL`, `GENCC_FILE_PREFIX`, `GENCC_SUBMISSION_FIELDS`, `GENCC_CLASSIFICATION_LEVELS` in `hvantk/skills/gencc/shared/constants.py`.
- Tests: `hvantk/skills/gencc/tests/test_gencc.py` (dataset class + builder-driven streamer tests), `test_drift_probe.py`.

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Confirm the raw TSV is present at `<raw_dir>/gencc-submissions-<YYYY-MM-DD>.tsv`. If absent, run `hvantk download gencc --output-dir <raw_dir>`.
3. Build via CLI: `hvantk reprocess gencc:submissions --raw-dir <dir> --output <path>.ht` (pass `--plugin-arg min_classification=Strong` to filter).
5. Sanity-check the output: row count plausible (per-submitter granularity gives more rows than ClinGen); key fields present; HGNC/MONDO prefixes stripped on at least one known row (e.g., BRCA1 → `hgnc_id == "1100"`, not `"HGNC:1100"`).
6. Run validation: `pytest hvantk/skills/gencc/tests -m hail`.

## 8. Update playbook

When GenCC publishes an updated snapshot (any download is effectively a new snapshot):

1. Re-download: `hvantk download gencc --output-dir <raw_dir> --overwrite`. Capture the new `Content-Length` in the PR description — that is what the drift probe treats as the content signal (see § 2). `Last-Modified` is recorded for reference only.
2. Diff the new TSV header against the previous fixture header. New columns alone are non-breaking — they will not appear in the built table unless added to `GENCC_SUBMISSION_FIELDS`. Removed/renamed columns require updating `GENCC_SUBMISSION_FIELDS`.
3. Drift-probe diff: regenerate `hvantk/skills/gencc/tests/drift_fingerprint.json` and inspect for column-list changes or hash changes.
4. If the fixture (`hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv`) is no longer representative (new classification value, new submitter under test), refresh it from a curated sub-sample.
5. Re-run the snapshot round-trip with `--regenerate-snapshots`; expected diffs: new optional columns, refreshed `submission_date` values. Unexpected diffs (changed key semantics, stripped prefixes leaking back in) require investigation before commit.
6. Open PR; reviewer checks the snapshot diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv`
- `drift_fingerprint`: `hvantk/skills/gencc/tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/gencc/tests -m hail`

> **Snapshot status:** seeded. `tests/snapshots/schema.json` and
> `tests/snapshots/sample_rows.json` are committed, and `test_builder.py` asserts the
> build against them. Regenerate after an intentional schema change with
> `pytest hvantk/skills/gencc/tests/test_builder.py --regenerate-snapshots`, then
> commit the result.
