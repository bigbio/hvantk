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
- Downloading the raw file. See `hvantk/commands/hgnc_downloader.py`.
- Cross-resource mapping logic. Lives in `hvantk/data/gene_mapper.py` (`GeneMapper`).
- Downstream consumers (ClinGen/GenCC streamers, gene-symbol resolution in PSROC). Those reference the built table by path.

## 2. Source identity

HGNC = HUGO Gene Nomenclature Committee. The complete-set TSV is the authoritative reference for current approved human gene symbols, IDs, and curated cross-references.

> Catalog gap: HGNC is **not yet registered** in `hvantk/resources/catalog.yaml` or any `hvantk/resources/registry/<domain>/datasets.json` (verified 2026-05-10). Until it is, the URL/version constants live in `hvantk/core/constants.py` (`HGNC_DOWNLOAD_URL`, `HGNC_INFO_URL`). Do **not** restate them here. When HGNC is added to the catalog, drop this paragraph and reference the catalog entry.

Stable provider notes the catalog will not capture:
- HGNC publishes a single rolling "complete set" (no dated versions in the URL); freshness is determined by the file's HTTP `Last-Modified` header.
- The TSV uses `\N` for missing in some columns and empty strings in others; both are treated as missing on import (see § 4).

## 3. Backend choice + reasoning

`hail`. Although HGNC is small (~43k rows) and `_conventions` § 3 allows pandas for mapping tables, every current consumer (`GeneMapper`, `clingen_streamer`, `gencc_streamer`, PSROC pipeline) joins it against Hail Tables, so producing a Hail Table avoids redundant materialization at every join site.

## 4. Raw format & gotchas

- File: tab-separated, single header row, ~50 columns. Header confirmed in `hvantk/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`.
- Imported with `hl.import_table(impute=False, missing="")` — all fields stay as strings; no type inference is attempted.
- Field renaming is driven by `HGNC_GENE_FIELDS` (`hvantk/core/constants.py`). Notable renames: `symbol → gene_symbol`, `name → gene_name`, `refseq_accession → refseq_id`, `orphanet → orphanet_id`, `date_approved_reserved → date_approved`. The builder only renames fields that are present in the input; columns absent from the upstream file are silently skipped.
- Pipe-separated multi-value fields are split into `array<str>` *after* renaming. The list of fields treated this way is `HGNC_PIPE_SEPARATED_FIELDS`. Empty/absent values become `[]`, not `missing`.
- `status` filter: by default the builder keeps only rows where `status == "Approved"`. Pass `include_withdrawn=True` to keep symbols, entry-withdrawn rows, etc. Withdrawn rows often have a populated `hgnc_id` but missing cross-references — joining to them silently produces nulls.
- Symbol history: `prev_symbols` and `alias_symbols` are needed to resolve legacy gene symbols. Downstream resolvers (e.g., `GeneMapper.resolve_symbol`) walk these arrays — do not strip them when selecting `--fields`.
- The builder does **not** strip the `HGNC:` prefix from `hgnc_id` (e.g., key is `"HGNC:1100"`, not `"1100"`). Streamers that ingest other sources (ClinGen, GenCC) explicitly strip the `HGNC:` prefix on *their* side before joining, so this asymmetry matters — see `hvantk/data/gene_disease_streamer.py`.

## 5. Output contract

Hail Table keyed by `hgnc_id` (string, with `HGNC:` prefix preserved). Schema is the source of truth — see `hvantk/tests/snapshots/hgnc/schema.json` (declared; regenerate via `--regenerate-snapshots` on first round-trip run).

Summary: one row per approved gene (≈43k in the live release; 5 in the fixture). Row fields fall into core identifiers (`hgnc_id`, `gene_symbol`, `gene_name`, `status`), symbol history arrays (`alias_symbols`, `alias_names`, `prev_symbols`, `prev_names`), cross-reference IDs (`ensembl_gene_id`, `entrez_id`, `uniprot_ids`, `refseq_id`, `ucsc_id`, `ccds_id`), classification (`locus_group`, `locus_type`, `gene_group`), location (`location`, `location_sortable`), clinical links (`omim_id`, `orphanet_id`, `gencc`, `mane_select`), and audit dates. Pipe-separated multi-value fields are arrays; the rest are scalars (mostly `tstr`).

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/hgnc/plugin.yaml` (drives loader registration and `hvantk drift hgnc:lookup`).
- Builder: `create_hgnc_gene_tb` in `hvantk/tables/table_builders.py` (uses `_create_table_base()` per `_conventions` § 4).
- CLI: `hvantk mktable hgnc` defined in `hvantk/commands/make_table_cli.py` (`mktable_hgnc`). Supports `--include-withdrawn`, `--fields`, `--overwrite`, `--export-tsv`.
- Constants: `HGNC_GENE_FIELDS`, `HGNC_PIPE_SEPARATED_FIELDS`, `HGNC_DOWNLOAD_URL`, `HGNC_INFO_URL` in `hvantk/core/constants.py`.
- Downloader: `hvantk/commands/hgnc_downloader.py` (CLI: `hvantk download hgnc`, wired in `hvantk/commands/download_cli.py`).
- Registry: **not registered** in `hvantk/tables/registry.py`. HGNC is built as a one-off lookup ahead of recipe runs, not as part of a batch recipe — register only if a real recipe-driven workflow demands it.
- Existing tests: assertion-based unit + GeneMapper tests in `hvantk/tests/test_hgnc_table_hail.py` (marked `hail` and `slow`). The snapshot round-trip test (see § 9) is a *new* file the agent should create on first run; do not extend `test_hgnc_table_hail.py` to do snapshot work — keep concerns separated.
- Downstream consumers (read-only): `hvantk/data/gene_mapper.py` (`GeneMapper` validates the table is keyed by `hgnc_id`), `hvantk/data/gene_disease_streamer.py` (and its `clingen_streamer` / `gencc_streamer` subclasses), `hvantk/psroc/pipeline.py`, `hvantk/utils/gene_aliases.py`, `hvantk/utils/geneset_io.py`.

## 7. Workflow steps

When invoked to build, refresh, or extend the HGNC table:

1. **Confirm the raw file is present.** If absent, point the user at `hvantk download hgnc --output <path>`; do not attempt to download from inside this workflow.
2. **Verify the header.** `head -1` the TSV and confirm every key in `HGNC_GENE_FIELDS` either exists or is acceptably missing. New upstream columns are non-breaking; *missing* expected columns mean the upstream schema drifted — stop and surface the diff.
3. **Build the Hail Table** by calling `create_hgnc_gene_tb(input_path, output_path, overwrite=…)` (Python API) or `hvantk mktable hgnc --raw-input … --output-ht … [--include-withdrawn] [--fields …]` (CLI). Both go through `_create_table_base`, so checkpointing and optional TSV export are handled.
4. **Sanity-check the output.** Confirm the table is keyed by `hgnc_id`, row count is in the expected range (~43k approved; +~5k if `--include-withdrawn`), and pipe-separated fields are arrays — not strings — for at least one known multi-value gene (e.g., BRCA1 → `alias_symbols` contains `BRCC1`).
5. **Run the snapshot round-trip test** (§ 9). If snapshots do not yet exist, create them with `pytest … --regenerate-snapshots`, review the diff, and commit alongside the builder change.
6. **Do not modify** the `HGNC:` prefix on `hgnc_id` keys. Downstream code relies on the prefix being preserved here and stripped at the join site.

## 8. Update playbook

Triggered when HGNC publishes an updated complete-set file or when an upstream schema change surfaces.

1. Re-download the raw TSV (`hvantk download hgnc --overwrite`). Capture the new `Last-Modified` header in the PR description — that is the de-facto version handle.
2. Diff the new TSV header against the previous fixture header (`diff <(head -1 old.tsv) <(head -1 new.tsv)`). New columns alone are non-breaking — they will not appear in the built table unless added to `HGNC_GENE_FIELDS`. Removed/renamed columns require updating `HGNC_GENE_FIELDS` (and possibly `HGNC_PIPE_SEPARATED_FIELDS`).
3. If the test fixture (`hvantk/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`) is no longer representative (e.g., a tested gene was withdrawn, a new pipe-separated field was added), regenerate it from the live file by sub-sampling the same gene set (`HGNC:1100`, `HGNC:1101`, `HGNC:4641`, plus the withdrawn row used by `test_create_hgnc_gene_tb_include_withdrawn`).
4. Re-run the round-trip test with `--regenerate-snapshots`. Expected diffs: new optional columns added to the schema; widened pipe-separated arrays; refreshed `date_modified` values in `sample_rows.json`. Unexpected diffs: changed `hgnc_id` keys, missing core fields (`gene_symbol`, `ensembl_gene_id`), changed `status` semantics — investigate before committing.
5. Update the catalog entry once it exists (see § 2). Until then, no version-string update is required because the URL is stable.
6. Re-run the assertion-based suite (`pytest hvantk/tests/test_hgnc_table_hail.py -m hail`) to confirm GeneMapper-side invariants still hold.

## 9. Validation contract

- `fixture`: `hvantk/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`
- `schema_snapshot`: `hvantk/tests/snapshots/hgnc/schema.json`
- `row_snapshot`: `hvantk/tests/snapshots/hgnc/sample_rows.json`
- `test_command`: `pytest hvantk/tests/test_hgnc_builder.py -m hail`

The snapshot directory and round-trip test file are declared but not yet created — initial run uses `pytest … --regenerate-snapshots` to seed them, per `_conventions` § 8.
