---
name: hvantk:resource-dbnsfp
description: dbNSFP variant functional annotation database — per-variant prediction scores, rankscores and population frequencies keyed by (locus, alleles).
status: provisional
backend: hail
domain: genomics
---

# dbNSFP

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the dbNSFP academic release TSV/BGZ → Hail Table builder (`dbnsfp:variants`): the per-variant functional-prediction table used as the score substrate across hvantk (PSROC, rerank, annotation).

Out of scope for this skill (per `_conventions` § 11):
- Acquisition. `plugin.yaml` declares no `lifecycle.download` — see § 2 and § 8; raw files are materialised externally.
- Score interpretation and circularity handling downstream. The `scores:` declaration in `plugin.yaml` is the input; the consumer is `hvantk.algorithms.rerank.provenance.resolve_arms`.

## 2. Source identity

dbNSFP (database of **N**on**S**ynonymous **F**unctional **P**redictions) aggregates prediction scores from many tools (SIFT, PolyPhen-2, CADD, REVEL, …) and population frequencies (gnomAD, ExAC, 1000Gp3, ESP6500) for essentially every possible non-synonymous SNV in the human genome. Landing page: <https://sites.google.com/site/jpopgen/dbNSFP>.

Catalog entry: `hvantk/skills/dbnsfp/catalog/datasets.json` (accession `dbNSFP_v4.9a`, GRCh38, academic-use licence, quarterly cadence). Do not restate URLs or versions here — query with `hvantk catalog show dbNSFP_v4.9a`.

> **The documented download path is broken.** Every `dbNSFP*.zip` the landing page links returns 404 — the S3 bucket answers `NoSuchBucket` — and the `database.liulab.science` mirror does not resolve. Tracked as issue #321. Acquisition currently has no working public route, so the raw file must be obtained out of band.

Releases are distinguished by an `a` (academic) or `c` (commercial) suffix; hvantk builds from the academic one.

## 3. Backend choice + reasoning

`hail`. The full academic release is ~45 GB compressed over ~84M rows (`catalog/datasets.json`), far past what pandas can hold, and every consumer joins it against variant-keyed Hail Tables. The builder imports with `min_partitions=200` by default for that reason.

## 4. Raw format & gotchas

These are the points that actually break a first-pass implementation. Read them before writing any import code.

- **The chromosome column is literally named `#chr`, with the leading hash.** `hl.import_table` takes the header line verbatim, so the field is `ht["#chr"]` — not `ht.chr`, and not `ht["chr"]`. The builder renames it (`builder.py:71-72`) and only then keys off `chr`. Skipping the rename yields `LookupError: Table instance has no field '#chr'` or a missing-field error downstream. The builder accepts an already-renamed `chr` column too, and raises `ValueError` if neither is present.
- **Contigs carry no `chr` prefix and must be prefixed for GRCh38.** The fixture's first data row starts `10\t47057\tC\tA…`; the default `reference_genome="GRCh38"` has contigs named `chr10`. Passing the bare value to `hl.parse_variant` raises `Invalid locus '10:83406' … Contig '10' is not in the reference genome 'GRCh38'`. The builder prefixes conditionally (`builder.py:77-82`) — `hl.if_else(_chr_str.lower().startswith("chr"), _chr_str, hl.str("chr") + _chr_str)` — so an already-prefixed file is left alone. **Note the prefix is added regardless of `reference_genome`**, and Hail's `GRCh37` uses bare contigs (`1`, `2`, …), so passing `reference_genome=GRCh37` produces `chr10` and fails to parse. The parameter exists, but the builder is effectively GRCh38-only until that is fixed.
- **Several column names contain characters that block attribute access**, notably `pos(1-based)` (and its `hg19_pos(1-based)` / `hg18_pos(1-based)` siblings). Always subscript: `ht["pos(1-based)"]`. The builder validates that `pos(1-based)`, `ref` and `alt` are all present before constructing the key.
- **The missing sentinel is `.`, not the empty string**, and `impute=False` — every column arrives as `tstr`. `hl.parse_float` maps `.` to missing, which is why an unscored variant stays missing rather than becoming `0.0` (which would read as "confidently benign").
- **The table is very wide** — 458 columns in the committed fixture, wider still in the full academic release. The builder groups work into a few `annotate` passes by *kind* — transcript-keyed `*_score` dicts, then per-variant rankscores, then prefix structs — rather than one pass over all 458 columns. Note it does issue a single `ht.annotate(**ann)` covering all ~44 score columns at once (`builder.py:136-145`), so that width is demonstrably fine; group by kind for readability, not to dodge a limit. Relatedly, do **not** pass `sep=" "` to `hl.import_table` — the file is tab-separated, and a space separator turns each row into hundreds of bogus fields and runs the test suite into a multi-minute timeout.
- **Multi-value fields are `;`-delimited**, aligned positionally with `Ensembl_transcriptid` (e.g. `aapos` = `411;408;445`). The builder splits `Ensembl_transcriptid` on `;` and zips it against each `*_score` / `CADD_phred` column to build a `dict<str, float64>` keyed by transcript. A single-valued score is broadcast to every transcript rather than left scalar.
- **`*_rankscore` columns are per-variant, not per-transcript** — one value each, cast to `float64`, never dicts. They are dbNSFP's only mutually comparable surface: all on 0–1, and the `_converted_` ones (SIFT, FATHMM, PROVEAN, LRT, MutationTaster, bStatistic) are inverted at source so that for *every* rankscore higher = more damaging. Raw scores are **not** direction-consistent — SIFT and FATHMM run the other way — so a feature axis built on raw scores is silently backwards for those tools.
- **Compression is ambiguous.** Files are distributed `.gz` but are often BGZF; `resolve_compression` (`hvantk/core/utils/file_utils.py`) detects which and sets `force_bgz` accordingly, with `auto_convert_bgz` available to rewrite a plain gzip in place.

## 5. Output contract

Hail Table keyed by `(locus, alleles)` — `hl.tlocus(reference_genome)` and `hl.tarray(hl.tstr)`, built via `hl.parse_variant` over a `chr:pos:ref:alt` string. Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) stamped with `schema_id="dbnsfp-v1"`.

The staging columns `variant_key`, `chr`, `pos(1-based)`, `ref`, `alt` are dropped after keying. Remaining fields fall into: transcript-keyed score dicts (`*_score`, `CADD_phred` → `dict<str, float64>`), per-variant `*_rankscore` floats (57 in the fixture), and population-frequency structs — `gnomAD`, `ExAC`, `1000Gp3`, `ESP6500`, `clinvar` (the default `group_prefixes`), each collecting every column sharing that prefix. Everything else stays `tstr`.

The schema snapshot is the source of truth: `hvantk/skills/dbnsfp/tests/snapshots/schema.json`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/dbnsfp/plugin.yaml`. Resolved as `get_registry().get_dataset("dbnsfp:variants")`; drives `hvantk drift dbnsfp:variants`.
- Builder: `build_dbnsfp_variants` in `hvantk/skills/dbnsfp/builder.py`, signature `(parsed_input, ctx, **params) -> AnnotationTable`. Params: `reference_genome` (default `"GRCh38"`), `min_partitions` (200), `force_bgz` (True), `parse_transcript_scores` (True), `group_prefixes` (list of str), `auto_convert_bgz` (False).
- CLI: `hvantk reprocess dbnsfp:variants --raw-dir <dir> --output <out>.ht [--plugin-arg reference_genome=GRCh38]`. Builds run through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).
- Helpers used: `resolve_compression` (`hvantk/core/utils/file_utils.py`), `get_row_fields` (`hvantk/core/utils/table_utils.py`).
- Drift probe: `fetch_fingerprint` in `hvantk/skills/dbnsfp/drift_probe.py` (see § 8).
- Training provenance: the `scores:` block in `plugin.yaml` (see § 10).
- Downstream consumers (read-only): `hvantk/algorithms/psroc/`, `hvantk/algorithms/rerank/`. The shared fixture is also read by `hvantk/tests/test_plugin_conformance.py` and a psroc integration test.

## 7. Workflow steps

When invoked to build, refresh, or extend the dbNSFP table:

1. **Confirm the raw file is present.** There is no downloader (§ 2); if absent, stop and say so rather than attempting a fetch — the advertised URLs are dead.
2. **Read the header before writing any field name** (`_conventions` § 10). `zcat <file> | head -1 | tr '\t' '\n' | head -20`. Confirm `#chr`, `pos(1-based)`, `ref`, `alt`.
3. **Build** via `build_dbnsfp_variants(parsed_input, ctx, **params)` (Python) or `hvantk reprocess dbnsfp:variants --raw-dir <dir> --output <path>.ht` (CLI). The builder constructs the table inline and returns an `AnnotationTable`; output and checkpointing belong to the reprocess pipeline.
4. **Sanity-check.** Key is `(locus, alleles)`; a locus renders as `chr10:47057`, not `10:47057`; at least one `*_rankscore` is `float64` and missing (not `0.0`) for an unscored variant; `gnomAD` is a struct, not a flat set of columns.
5. **Run the snapshot round-trip test** (§ 9). Regenerate only for an intentional change, and review the diff before committing.
6. **Do not** flatten the prefix structs or scalarise the transcript dicts — downstream selectors address `gnomAD.*` and index score dicts by transcript ID.

## 8. Update playbook

Triggered when dbNSFP advertises a new release, or when `hvantk drift dbnsfp:variants` reports a change.

1. The drift probe scrapes the **release list advertised on the landing page** (`drift_probe.py`, `probe_version` 2), not any data file. It records the release set under `headers` and the highest academic release under `source_version`. It **never hashes the page body**: Google Sites re-renders per request (352,830 vs 352,716 bytes on two consecutive fetches), so a body hash would flag drift on every run. It fails closed — `DriftProbeError` — if no release matches, rather than baking an empty baseline.
   - Detects: a new release being advertised. Cannot detect: an in-place change to an archive's contents, or the download links being repaired (the markup names the same archives either way).
2. Obtain the new academic release out of band (§ 2) and re-run step 2 of § 7 against it. A header diff is the real schema signal.
3. New score columns are non-breaking — they flow through as `tstr`, or join a prefix struct if they match one. Removed or renamed `pos(1-based)` / `ref` / `alt` / `#chr` are breaking and raise from the builder.
4. If the fixture (`hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz`) is no longer representative, regenerate it by sub-sampling the same chr10 variant window the snapshot keys use (`chr10:47057`–`chr10:47059`). It is a **shared** fixture — changing it moves `hvantk/tests/test_plugin_conformance.py` and a psroc test too.
5. Re-run the round-trip test with `--regenerate-snapshots`, review, commit alongside the builder change. Expected diffs: new columns, wider score dicts. Unexpected: changed keys, a rankscore that stopped being `float64`, a struct that flattened.
6. Update the `scores:` block if a predictor was added or retrained (§ 10), and refresh `catalog/datasets.json` for the new accession.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (paths resolved relative to `hvantk/skills/dbnsfp/`):

- `fixture`: `../../tests/testdata/raw/dbnsfp` — the **shared** form (`_conventions` § 9). The file is `hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz`, referenced in place rather than duplicated because `hvantk/tests/test_plugin_conformance.py` and a psroc integration test read the same bytes.
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/dbnsfp/tests -m hail`

Snapshot status: **seeded**. `tests/test_builder.py` asserts schema and six sample rows (real keys taken from a build of the fixture — `collect_sample_rows` raises `KeyError` for an invented one). `tests/test_dbnsfp.py` exercises the build through `run_builder_for_spec`; `tests/test_drift_probe.py` covers the probe. Regenerate with `pytest hvantk/skills/dbnsfp/tests/test_builder.py -m hail --regenerate-snapshots`, then commit.

## 10. Cross-reference notes

`plugin.yaml` declares `trained_on` for 55 of dbNSFP's 57 rankscore predictors — which curated database, simulated allele set, or population resource each was fit on. dbNSFP is the reason the mechanism exists: roughly half its predictors are supervised on ClinVar or HGMD, so against a curated-database label they are partly circular, and a purely statistical filter would *reward* that circularity rather than catch it.

Three states: `[]` (nothing label-derived, e.g. `phyloP100way_vertebrate_rankscore`), an explicit source list (`REVEL_rankscore: [HGMD, ClinVar]`), and **omitted**, which means unknown and is never treated as clean. `MutFormer_rankscore` and `PHACTboost_rankscore` are deliberately omitted — both are plausibly supervised on curated disease variants, but their training sets could not be confirmed, and a wrong "clean" declaration is more damaging than an honest "unknown".

See `_conventions` § 6 for the contract and `hvantk.algorithms.rerank.provenance.resolve_arms` for the consumer. The list is a claim about training data, not about quality; correct it as tools are retrained.
