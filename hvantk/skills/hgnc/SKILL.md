---
name: hvantk:resource-hgnc
description: HGNC gene nomenclature lookup table — authoritative human gene symbols, IDs, and cross-references keyed by hgnc_id.
status: provisional
backend: hail
domain: mapping
---

# HGNC

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the HGNC complete-set TSV → Hail Table builder used as the canonical human-gene lookup across hvantk (gene-symbol ↔ HGNC ID ↔ Ensembl/Entrez/UniProt mapping; symbol-history resolution).

Out of scope for this skill (per `_conventions` § 11):
- Downloading the raw file. See `hvantk/skills/hgnc/cli.py` (`download_dataset` lifecycle entry point; CLI `hvantk hgnc-download`).
- Cross-resource mapping logic. Lives in `hvantk/skills/hgnc/streamers.py` (`HGNCGeneCatalogStreamer`, which absorbed the retired `GeneMapper`).
- Downstream consumers (ClinGen/GenCC streamers, gene-symbol resolution in PSROC). Those reference the built table by path.

## 2. Source identity

HGNC = HUGO Gene Nomenclature Committee. The complete-set TSV is the authoritative reference for current approved human gene symbols, IDs, and curated cross-references.

> Catalog gap: HGNC is **not yet registered** in any plugin's `catalog/datasets.json` or in `hvantk/resources/registry/genomics/datasets.json` (verified 2026-05-10; re-check with `hvantk catalog search HGNC`). Until it is, the URL/version constants live in `hvantk/skills/hgnc/shared/constants.py` (`HGNC_DOWNLOAD_URL`, `HGNC_INFO_URL`). Do **not** restate them here. When HGNC is added to the catalog, drop this paragraph and reference the catalog entry.

Stable provider notes the catalog will not capture:
- HGNC publishes a single rolling "complete set" (no dated versions in the URL). Freshness is **not** determined by `Last-Modified`: HGNC re-publishes byte-identical content under a fresh timestamp, so the drift probe's content signal is the HTTP `Content-Length` header (`extras.content_length`) instead. Across all 8 committed fingerprints from 2026-05-16 to 2026-08-27 the checksum (`b8cfde58…`) never moved while `Last-Modified` moved on every regeneration — and the checksum only hashes the column-header line (a SCHEMA signal), so it could not have caught a content-only change either. `Last-Modified` is still recorded, under `informational` (excluded from drift comparison), and `source_version` is `null`. The probe (`probe_version` 2) fails closed — raises `DriftProbeError` — if the server omits `Content-Length`.
- The TSV uses `\N` for missing in some columns and empty strings in others; both are treated as missing on import (see § 4).

## 3. Backend choice + reasoning

`hail`. Although HGNC is small (~43k rows) and `_conventions` § 3 allows pandas for mapping tables, every current consumer (`HGNCGeneCatalogStreamer`, the ClinGen/GenCC gene-disease streamers, PSROC pipeline) joins it against Hail Tables, so producing a Hail Table avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: tab-separated, single header row, ~50 columns. Header confirmed in `hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`.
- Imported with `hl.import_table(impute=False, missing="")` — all fields stay as strings; no type inference is attempted.
- Field renaming is driven by `HGNC_GENE_FIELDS` (`hvantk/skills/hgnc/shared/constants.py`). Notable renames: `symbol → gene_symbol`, `name → gene_name`, `refseq_accession → refseq_id`, `orphanet → orphanet_id`, `date_approved_reserved → date_approved`. The builder only renames fields that are present in the input; columns absent from the upstream file are silently skipped.
- Pipe-separated multi-value fields are split into `array<str>` *after* renaming. The list of fields treated this way is `HGNC_PIPE_SEPARATED_FIELDS`. Empty/absent values become `[]`, not `missing`.
- `status` filter: by default the builder keeps only rows where `status == "Approved"`. Pass `include_withdrawn=True` to keep symbols, entry-withdrawn rows, etc. Withdrawn rows often have a populated `hgnc_id` but missing cross-references — joining to them silently produces nulls.
- Symbol history: `prev_symbols` and `alias_symbols` are needed to resolve legacy gene symbols. Downstream resolvers (e.g., `HGNCGeneCatalogStreamer.resolve_symbol`) walk these arrays — do not strip them when selecting `--fields`.
- The builder does **not** strip the `HGNC:` prefix from `hgnc_id` (e.g., key is `"HGNC:1100"`, not `"1100"`). Streamers that ingest other sources (ClinGen, GenCC) explicitly strip the `HGNC:` prefix on *their* side before joining, so this asymmetry matters — see `hvantk/core/streamers/gene_disease_table.py` (`GeneDiseaseTableStreamer`) and the per-plugin subclasses in `hvantk/skills/clingen/streamers.py` / `hvantk/skills/gencc/streamers.py`.

## 5. Output contract

Hail Table keyed by `hgnc_id` (string, with `HGNC:` prefix preserved). Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) wrapping the Hail Table with provenance (`schema_id="hgnc-lookup-v1"`). Schema is the source of truth — see `hvantk/skills/hgnc/tests/snapshots/schema.json` (declared in `plugin.yaml`; not yet seeded — regenerate via `--regenerate-snapshots` on first round-trip run).

Summary: one row per approved gene (≈43k in the live release; 5 in the fixture). Row fields fall into core identifiers (`hgnc_id`, `gene_symbol`, `gene_name`, `status`), symbol history arrays (`alias_symbols`, `alias_names`, `prev_symbols`, `prev_names`), cross-reference IDs (`ensembl_gene_id`, `entrez_id`, `uniprot_ids`, `refseq_id`, `ucsc_id`, `ccds_id`), classification (`locus_group`, `locus_type`, `gene_group`), location (`location`, `location_sortable`), clinical links (`omim_id`, `orphanet_id`, `gencc`, `mane_select`), and audit dates. Pipe-separated multi-value fields are arrays; the rest are scalars (mostly `tstr`).

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/hgnc/plugin.yaml` (the loader auto-resolves the dataset via `get_registry().get_dataset("hgnc:lookup")`; drives `hvantk drift hgnc:lookup`).
- Builder: `build_hgnc_gene_lookup` in `hvantk/skills/hgnc/builder.py`. Signature is `build_hgnc_gene_lookup(parsed_input, ctx, *, include_withdrawn=False, fields=None) -> AnnotationTable`. It builds the table inline with `hl.import_table` + renames + pipe-field splitting (no `_create_table_base`, no `output_path`/`overwrite` kwargs).
- CLI: `hvantk reprocess hgnc:lookup --raw-dir <dir> --output <path>.ht`. Builder kwargs (`include_withdrawn`, `fields`) flow through `--plugin-arg key=value`. Top-level builds run through `run_builder_for_spec` in `hvantk/core/plugin/run_builder.py`.
- Constants: `HGNC_GENE_FIELDS`, `HGNC_PIPE_SEPARATED_FIELDS`, `HGNC_DOWNLOAD_URL`, `HGNC_INFO_URL` in `hvantk/skills/hgnc/shared/constants.py`.
- Downloader: `hvantk/skills/hgnc/cli.py` (`download_dataset` lifecycle entry point, wired via `plugin.yaml` `lifecycle.download`; CLI command `hvantk hgnc-download` via the `cli:` block).
- Streamer: `HGNCGeneCatalogStreamer` in `hvantk/skills/hgnc/streamers.py` (subclass of `GeneCatalogStreamer` in `hvantk/core/streamers/gene_catalog.py`; validates the table is keyed by `hgnc_id` and absorbs the retired `GeneMapper`/`gene_aliases` logic).
- Existing tests: `hvantk/skills/hgnc/tests/test_downloader.py` and `test_drift_probe.py`. The snapshot round-trip test (see § 9) is declared in `plugin.yaml` but not yet seeded; create it on first round-trip run.
- Downstream consumers (read-only): the ClinGen/GenCC streamers in `hvantk/skills/clingen/streamers.py` and `hvantk/skills/gencc/streamers.py` (over `GeneDiseaseTableStreamer`), and `hvantk/algorithms/psroc/pipeline.py`.

## 7. Workflow steps

When invoked to build, refresh, or extend the HGNC table:

1. **Confirm the raw file is present.** If absent, point the user at `hvantk hgnc-download --output <path>`; do not attempt to download from inside this workflow.
2. **Verify the header.** `head -1` the TSV and confirm every key in `HGNC_GENE_FIELDS` either exists or is acceptably missing. New upstream columns are non-breaking; *missing* expected columns mean the upstream schema drifted — stop and surface the diff.
3. **Build the table** by calling `build_hgnc_gene_lookup(parsed_input, ctx, include_withdrawn=…, fields=…)` (Python API) or `hvantk reprocess hgnc:lookup --raw-dir <dir> --output <path>.ht [--plugin-arg include_withdrawn=true] [--plugin-arg fields=…]` (CLI). The builder constructs the table inline and returns an `AnnotationTable`; checkpointing/output is handled by the reprocess pipeline, not the builder.
4. **Sanity-check the output.** Confirm the table is keyed by `hgnc_id`, row count is in the expected range (~43k approved; +~5k if `--include-withdrawn`), and pipe-separated fields are arrays — not strings — for at least one known multi-value gene (e.g., BRCA1 → `alias_symbols` contains `BRCC1`).
5. **Run the snapshot round-trip test** (§ 9). If snapshots do not yet exist, create them with `pytest … --regenerate-snapshots`, review the diff, and commit alongside the builder change.
6. **Do not modify** the `HGNC:` prefix on `hgnc_id` keys. Downstream code relies on the prefix being preserved here and stripped at the join site.

## 8. Update playbook

Triggered when HGNC publishes an updated complete-set file or when an upstream schema change surfaces.

1. Re-download the raw TSV (`hvantk hgnc-download --overwrite`). Capture the new `Content-Length` in the PR description — that is what the drift probe treats as the content signal (see § 2). `Last-Modified` is recorded for reference only and does not indicate whether content actually changed.
2. Diff the new TSV header against the previous fixture header (`diff <(head -1 old.tsv) <(head -1 new.tsv)`). New columns alone are non-breaking — they will not appear in the built table unless added to `HGNC_GENE_FIELDS`. Removed/renamed columns require updating `HGNC_GENE_FIELDS` (and possibly `HGNC_PIPE_SEPARATED_FIELDS`).
3. If the test fixture (`hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`) is no longer representative (e.g., a tested gene was withdrawn, a new pipe-separated field was added), regenerate it from the live file by sub-sampling the same gene set (`HGNC:1100`, `HGNC:1101`, `HGNC:4641`, plus a withdrawn row to exercise `include_withdrawn`).
4. Re-run the round-trip test with `--regenerate-snapshots`. Expected diffs: new optional columns added to the schema; widened pipe-separated arrays; refreshed `date_modified` values in `sample_rows.json`. Unexpected diffs: changed `hgnc_id` keys, missing core fields (`gene_symbol`, `ensembl_gene_id`), changed `status` semantics — investigate before committing.
5. Update the catalog entry once it exists (see § 2). Until then, no version-string update is required because the URL is stable.
6. Re-run the plugin test suite (`pytest hvantk/skills/hgnc/tests -m hail`) to confirm streamer-side invariants still hold.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (all paths plugin-relative under `hvantk/skills/hgnc/`):

- `fixture`: `tests/testdata/raw/hgnc/hgnc_test_sample.tsv`
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/hgnc/tests -m hail`

The snapshot directory and round-trip test file are declared but not yet created — initial run uses `pytest … --regenerate-snapshots` to seed them, per `_conventions` § 8.

> **Snapshot status:** seeded. `tests/snapshots/schema.json` and
> `tests/snapshots/sample_rows.json` are committed, and `test_builder.py` asserts the
> build against them. Regenerate after an intentional schema change with
> `pytest hvantk/skills/hgnc/tests/test_builder.py --regenerate-snapshots`, then
> commit the result.
