---
name: hvantk:resource-peptideatlas-phospho
description: Download and parse a PeptideAtlas Human Phospho build into a wide intermediate TSV consumed by the hvantk PTM pipeline.
status: provisional
backend: pandas
domain: proteomics
---

# PeptideAtlas Human Phospho resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

- **Status:** provisional. Downloader + parser + downloader tests + drift-probe placeholder are migrated under the plugin folder. Builder round-trip snapshots (`schema.json`, `sample_rows.json`) and the fixture are seeded; `test_builder.py` asserts against both.
- **In scope:** one PeptideAtlas human phospho build at a time (default: `(202512, 606)` — see `PEPTIDEATLAS_LATEST_BUILD_DATE` / `PEPTIDEATLAS_LATEST_BUILD_ID` in `hvantk/skills/peptideatlas/phospho/shared/constants.py`). Output is the wide intermediate TSV consumed by `hvantk/algorithms/ptm/pipeline.py` via its `peptideatlas_tsv` argument.
- **Out of scope:** non-human PeptideAtlas builds; non-phospho PTM atlases on PeptideAtlas (those would land as sibling datasets under `hvantk/skills/peptideatlas/<ptm-type>/`); any Hail-Table or AnnData representation — the only downstream consumer today reads the TSV directly.

## 2. Source identity

- **Provider:** PeptideAtlas (<https://peptideatlas.org/>).
- **Catalog entry:** TODO. No PeptideAtlas entry exists yet in any plugin's `catalog/datasets.json` (verify with `hvantk catalog search peptideatlas`). The pinned build coordinates live in `hvantk/skills/peptideatlas/phospho/shared/constants.py` as `PEPTIDEATLAS_PHOSPHO_BASE_URL`, `PEPTIDEATLAS_LATEST_BUILD_DATE`, and `PEPTIDEATLAS_LATEST_BUILD_ID`. When a dedicated catalog entry lands, point this plugin's `source.catalog_ref` at it and remove this TODO.
- Download URL is composed by `PeptideAtlasPhosphoDataset.from_build(build_date, build_id)` as `{base}/{build_date}/atlas_build_{build_id}.tsv.zip`. The dataset class enforces HTTPS + `*.peptideatlas.org` host validation before any GET.

## 3. Backend choice + reasoning

**`backend: pandas`, `domain: proteomics`.** Per `_conventions` § 3, Hail Tables / MatrixTables are for variant- or expression-keyed data joined into the wider hvantk genomics graph; AnnData is for sparse single-cell or sample x gene matrices. PeptideAtlas phospho output is neither: it is a narrow site-level table (one row per `(accession, position)`) keyed by UniProt accession and is consumed downstream by the PTM pipeline (`hvantk/algorithms/ptm/pipeline.py`) as a TSV file. A Hail Table would force a Spark session for a sub-100 MB lookup table. A pandas DataFrame is the natural fit; the plugin builder `build_peptideatlas_phospho(parsed_input, ctx, **params)` loads the intermediate TSV via `pd.read_csv(... sep="\t" ...)` and wraps it as an `AnnotationTable`.

## 4. Raw format & gotchas

The raw upstream artefact is a single TSV zip (e.g. `atlas_build_606.tsv.zip`, multiple GB) containing PeptideAtlas's relational dump:

- `biosequence.tsv` — proteins (`biosequence_id`, accession, gene name, sequence).
- `protein_identification.tsv` — protein confidence; `presence_level_id == 1` means *canonical*. The parser filters by this and drops `DECOY_` / `CONTAM_` prefixed accessions in `biosequence.tsv`. This table is **optional** to the parser — it is not in `parse_peptideatlas_zip`'s `required_tables`. Three cases:
  1. Missing or misnamed: the parse runs without the canonical filter, so non-canonical isoforms pass through. The only sign is a `logger.warning` (tracked in #416). Real builds ship the table.
  2. Present but with no row where `presence_level_id == "1"`: `parse_peptideatlas_zip` raises `ValueError` instead of silently returning an unfiltered or empty table.
  3. Present and populated: only its canonical `biosequence_id`s pass the filter.
- `peptide_instance.tsv` — distinct peptides + `n_observations` (sum of PSMs across experiments, over every form of the peptide, unmodified included).
- `peptide_mapping.tsv` — peptide-to-protein coordinates (`start_in_biosequence`); 15M+ rows; must be streamed, not slurped.
- `modified_peptide_instance.tsv` — modified peptide sequences in bracket notation, one row per modified form (sequence × charge), each with its own `n_observations`; 3M+ rows.
- **Required columns.** Every table above also has required columns; `_iter_tsv_from_zip` raises `ValueError` if the header is missing any of them. As the code has them today: `protein_identification.tsv` needs `biosequence_id`, `presence_level_id`; `biosequence.tsv` needs `biosequence_id`, `biosequence_accession`, `biosequence_gene_name`, `biosequence_seq`; `peptide_instance.tsv` needs `peptide_instance_id`; `peptide_mapping.tsv` needs `peptide_instance_id`, `matched_biosequence_id`, `start_in_biosequence`; `modified_peptide_instance.tsv` needs `peptide_instance_id`, `modified_peptide_sequence`, `n_observations`.

Modification notation gotchas (see `_extract_phospho_offsets` in `shared/datasets.py`):

- Text form: `S[Phospho]`, `T[Phospho]`, `Y[Phospho]`, or `S[Phospho:80]` (substring match handles the latter).
- Numeric form: `S[167]`, `T[181]`, `Y[243]` — the total modified-residue mass (residue + HPO3 ≈ 80). The parser only treats numeric brackets as phospho when the preceding residue is S/T/Y and the mass matches within ±1.0 Da.
- ProForma delta-mass form: `S[+79.966]`, `T[+79.966]`, `Y[+79.966]` (±0.005 Da around `_PHOSPHO_MOD_MASS` = 79.96633 Da: wide enough for the 2-decimal `+79.97`, tight enough to exclude sulfation, 79.9568 Da) and `S[UNIMOD:21]` / `T[UNIMOD:21]` / `Y[UNIMOD:21]`. Accepted defensively: a 1.55M-row sample of build 606, the live build this parser targets, showed only `[Phospho]` (see the live-build check below).
- N-terminal labels like `[TMT6plex]-` precede the first residue and must be skipped (no preceding amino acid).
- Lowercase letters, digits, and dashes in the sequence are ignored.

Aggregation gotcha: the same `(accession, position)` can be observed via multiple distinct peptides and via several modified forms of one peptide. `parse_peptideatlas_zip` sums each phospho form's own `modified_peptide_instance.n_observations` over all of them. It must not use the parent `peptide_instance.n_observations`: that count also covers the unmodified and other forms, and adding it once per form inflated 91% of the real build's 259,932 sites (82× in total, median 16× per site; #425). A phospho form that maps to a kept protein and has no integer count fails the parse. In build 606 the forms' counts sum exactly to the peptide's for all 381,252 peptides (42,499 forms are unmodified). Tables built before the fix (plugin version 0.1.0) carry the inflated counts; `hvantk download peptideatlas-phospho` returns an existing parsed TSV unchanged, so delete it (or pass `--overwrite`, which also downloads the zip again) before rebuilding. `peptide_mapping` repeats a (peptide, protein, start) row 508,098 times in build 606, but only for proteins the parser drops (decoys, contaminants, non-canonical isoforms): the 389,497 kept mappings are all distinct, so no site is counted twice through a repeated mapping. Re-check that on a new build, since `mapping_by_pi` keeps repeats.

**Confirmed against a live build (202512 / 606, checked 2026-10-06).** Downloaded and parsed the full ~549 MB zip on the cluster to check the facts above against real data (no rows retained — see the licence note in § 9):

- Zip members in this build are named *exactly* `biosequence.tsv`, `peptide_instance.tsv`, `peptide_mapping.tsv`, `modified_peptide_instance.tsv`, `protein_identification.tsv` — `_find_table_in_zip`'s exact-match branch fires; the substring-match fallback (`atlas_build_<id>_biosequence.tsv`-style names) is not exercised by this build but may be by others.
- Every header line *and* every data line ends with one extra delimiter, producing a spurious trailing empty-named column in every table. `csv.DictReader` absorbs this harmlessly, but a reader adding a column by position rather than by name would be off by one.
- Missing values are written as the literal string `\N`, not an empty field, throughout every table the builder touches; `_iter_tsv_from_zip` normalizes every such value to `""` before it reaches the parser.
- A phospho form that maps to a kept protein (non-decoy/non-contaminant, and canonical when the canonical filter is on) must have an integer `n_observations`; `\N` (i.e. `""`) or any other non-integer value raises `ValueError`. A form that is not phospho, or maps to no kept protein, is never read for its count, so it may carry `\N` freely.
- `peptide_instance.tsv` and `peptide_mapping.tsv` each carry two near-duplicate columns in this build: `lowest_n_missed_cleavages ` (trailing space) and `lowest_n_missed_cleavages` (no trailing space) are both present as distinct header names. Neither is read by the parser today, but a future column addition must match the exact spelling, trailing space included.
- The canonical flag is exactly what the parser assumes: `protein_presence_level.tsv` labels `presence_level_id == 1` as `"canonical"`, and in the live build 11,852 of ~223,070 `protein_identification.tsv` rows carried that value.
- Only the text notation `<residue>[Phospho]` was observed for phospho sites — a 1.55M-row / 300 MB sample of `modified_peptide_sequence` contained zero numeric-mass brackets (`[167]`, `[181]`, `[243]`) and zero `[Phospho:N]`-style matches. The numeric-mass branch in `_extract_phospho_offsets` is therefore untested by this build's data, though it remains in the parser for older/other builds. Many *non*-phospho bracket modifications co-occur in the same column and are correctly ignored: `TMT6plex`, `iTRAQ4plex`, `iTRAQ8plex`, `TMTpro`, `Carbamidomethyl`, `Deamidated`, `Label:…`, `Oxidation`, `Dimethyl:…`, `Gln->pyro-Glu`, `Glu->pyro-Glu`, `Acetyl`, `Pyro-carbamidomethyl` — including a trailing (C-terminal/side-chain) form such as `K[TMT6plex]`, confirming the "only when `last_aa` is S/T/Y" guard is exercised on internal/trailing brackets, not just the N-terminal ones.
- End-to-end `parse_raw_dir` over the live build: 78 s wall time, ~450 MB peak RSS (`/usr/bin/time -v` `Maximum resident set size`), 259,932 distinct phospho sites from 597,904 biosequences (11,852 canonical) / 381,252 peptide instances / 15,489,152 peptide mappings (373,403 resolved to a kept protein) / 3,097,865 modified peptide instances (2,607,316 carrying a phospho site). Output `amino_acid` values were exactly `{S, T, Y}`, no `DECOY_`/`CONTAM_` accession leaked through, and the well-known p53 CDK-phosphorylation site `P04637` Ser315 was present — confirming the streaming joins and filters hold at full scale, not just against the mock tables in `test_phospho.py`.

## 5. Output contract

- **File:** wide TSV at the path returned by `download_dataset` / `parse_raw_dir`, named `peptideatlas-phospho-<build_date>-<build_id>.tsv` (e.g. `peptideatlas-phospho-202512-606.tsv`).
- **Shape:** one row per `(accession, position)` phospho site on a canonical protein.
- **Columns** (fixed order, see `_TSV_COLUMNS`): `accession`, `gene_symbol`, `position`, `description`, `amino_acid`, `ensembl_xrefs`, `sequence_length`, `n_observations`, `source_db`, `evidence_type`. `source_db` is always `"PeptideAtlas"`; `evidence_type` is always `"mass_spectrometry"`; `ensembl_xrefs` is currently always empty.
- **Builder return value:** `build_peptideatlas_phospho` returns an `AnnotationTable` (via `AnnotationTable.from_pandas(...)`) wrapping `pd.read_csv(<tsv>, sep="\t", dtype=str)` — all columns as strings so downstream consumers don't get silent numeric coercion on `position` / `n_observations`. Provenance is stamped from `ctx.provenance(schema_id="peptideatlas-phospho-v1")`.

## 6. hvantk integration points

- **Dataset class + parser:** `PeptideAtlasPhosphoDataset`, `parse_peptideatlas_zip`, `write_intermediate_tsv`, `parse_raw_dir` in `hvantk/skills/peptideatlas/phospho/shared/datasets.py`.
- **Builder:** `build_peptideatlas_phospho(parsed_input, ctx, **params) -> AnnotationTable` in `hvantk/skills/peptideatlas/phospho/builder.py` (declared in `plugin.yaml` under `datasets[].builder`). A legacy `build_peptideatlas_phospho_tb(input_path, output_path, ...)` helper also lives in that module but is NOT the loader entry point.
- **Downloader CLI:** `download_cmd` in `hvantk/skills/peptideatlas/phospho/cli.py` (declared in the manifest's `cli:` block as `peptideatlas-phospho-download`; the loader strips the `-download` suffix and binds it under the `download` group, so the invocation is `hvantk download peptideatlas-phospho`).
- **Lifecycle entry points:** `download_dataset` and `parse_raw_dir` (loader-wired via `lifecycle.download` + `lifecycle.parse` in `plugin.yaml`).
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/peptideatlas/phospho/drift_probe.py` (HEAD against the pinned build's zip URL).
- **Downstream consumer:** `hvantk/algorithms/ptm/pipeline.py` (`PTMBuildConfig.peptideatlas_tsv`) — reads the intermediate TSV produced here and maps PTM sites to genomic coordinates. Exposed at the user-facing level by `hvantk ptm build` (`--peptideatlas-tsv`).
- **Plugin manifest:** `hvantk/skills/peptideatlas/plugin.yaml` (compound dataset key `peptideatlas:phospho`).
- **Tests:** `hvantk/skills/peptideatlas/phospho/tests/` (parser unit tests + drift-probe sanity test).

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update the PeptideAtlas phospho intermediate:

1. **End-to-end (recommended).** Run the full download -> parse -> build chain through the plugin loader:

   ```bash
   hvantk reprocess peptideatlas:phospho --raw-dir /data/peptideatlas --output /out/peptideatlas-phospho.parquet
   ```

   The loader auto-resolves the dataset from `plugin.yaml` (`get_registry().get_dataset("peptideatlas:phospho")`) and runs the build through `run_builder_for_spec`. The `lifecycle.download` (`download_dataset`) and `lifecycle.parse` (`parse_raw_dir`) entry points run first, then the plugin builder `build_peptideatlas_phospho`.
2. **Download only.** Either via the standalone CLI (`hvantk download peptideatlas-phospho -o /data/peptideatlas`) or the lifecycle entry point `download_dataset(raw_dir=...)`. Both produce `<raw_dir>/atlas_build_<id>.tsv.zip` *and* the parsed `<raw_dir>/peptideatlas-phospho-<date>-<id>.tsv`.
3. **(Lifecycle) parse-only step.** `parse_raw_dir(raw_dir=..., output_path=...)` re-parses an existing zip from `raw_dir` into a fresh intermediate TSV — used when downstream code wants the TSV at a different path than the dataset class's default.
4. **Builder.** `build_peptideatlas_phospho(parsed_input, ctx, **params)` loads the intermediate TSV (the path produced by `parse_raw_dir`) as a pandas DataFrame and returns an `AnnotationTable`.
5. **Validate.** `pytest hvantk/skills/peptideatlas/phospho/tests` — parser unit tests + drift-probe sanity. Builder round-trip snapshot test is TODO (see § 9).

## 8. Update playbook

PeptideAtlas releases a new human phospho build once or twice per year. When a new build is announced:

1. Run the drift probe: `python -c "from hvantk.skills.peptideatlas.phospho.drift_probe import fetch_fingerprint; print(fetch_fingerprint())"`. A change in `source_version` (the `Last-Modified` header) is the trigger.
2. Update the pinned build coordinates in `hvantk/skills/peptideatlas/phospho/shared/constants.py` (`PEPTIDEATLAS_LATEST_BUILD_DATE`, `PEPTIDEATLAS_LATEST_BUILD_ID`).
3. Re-download (`hvantk download peptideatlas-phospho -o /data/peptideatlas --overwrite`) and spot-check the row count against the previous build.
4. If any new modification notation or table appears in the dump, document it in § 4.
5. Re-regenerate `tests/drift_fingerprint.json` with the new build's filename + headers.
6. If a builder snapshot exists (TODO), run `pytest hvantk/skills/peptideatlas/phospho/tests --regenerate-snapshots` and inspect the diff.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/peptideatlas/phospho/tests/testdata/raw/peptideatlas-phospho/` (seeded) — a directory holding `atlas_build_606-synthetic.tsv.zip` plus a `README.md` describing it. This is a **synthetic miniature raw build**, not a truncation of a real PeptideAtlas build — the real `atlas_build_*.tsv.zip` is ~549 MB (`content_length` in `tests/drift_fingerprint.json`) and is not vendored into this repo. Its five tables (`biosequence.tsv`, `peptide_instance.tsv`, `peptide_mapping.tsv`, `modified_peptide_instance.tsv`, `protein_identification.tsv`) carry the real column headers and real `<residue>[Phospho]` modification notation confirmed against a live build (§ 4); every row is fabricated. The round-trip test runs the real `parse_raw_dir` over this zip, so it now covers upstream zip parsing (table joins, offset extraction, canonical/DECOY_/CONTAM_ filtering, observation aggregation) in addition to the intermediate-TSV -> `AnnotationTable` contract — `test_phospho.py`'s own hand-written mock-zip tests remain the place for notation edge cases (numeric-mass brackets, N-terminal labels, canonical filtering in isolation) this fixture does not need to re-cover.
- **schema_snapshot:** `hvantk/skills/peptideatlas/phospho/tests/snapshots/schema.json` (seeded).
- **row_snapshot:** `hvantk/skills/peptideatlas/phospho/tests/snapshots/sample_rows.json` (seeded, keyed on `(accession, position)` for four of the fixture's six phospho sites — an aggregated two-peptide site, a two-site-in-one-peptide site, one arm of a multi-mapping peptide, and a site seen in two phospho forms of one peptide).
- **command:** `pytest hvantk/skills/peptideatlas/phospho/tests`.
- **drift_fingerprint:** `hvantk/skills/peptideatlas/phospho/tests/drift_fingerprint.json` — a real baseline captured from a live probe run (build 202512 / 606, fetched 2026-07-27), not a placeholder. Refresh via the update playbook.

The plugin manifest declares these paths and every one now resolves to a real committed file — `peptideatlas:phospho` is off the `KNOWN_INCOMPLETE` ledger in `hvantk/tests/test_plugin_contract_artifacts.py`.

**Round-trip test:** `hvantk/skills/peptideatlas/phospho/tests/test_builder.py`. It first runs `parse_raw_dir` over the fixture directory, then `build_peptideatlas_phospho` over the resulting intermediate TSV, so one test now grades both plugin stages against one committed fixture. It deliberately does NOT go through `hvantk.tests._snapshot_utils.regenerate_snapshots` — that helper only dispatches on `anndata.AnnData` or a Hail-table fallback, and this builder returns a pandas-backed `AnnotationTable` with no native Hail/AnnData representation (see § 3). Routing through the Hail branch would only work via the `to_hail()` escape hatch, at the cost of a full Hail/Spark session (`@pytest.mark.hail`) for a builder that never touches Hail, and it also distorts dtypes in the process (an all-empty `ensembl_xrefs` column round-trips as Hail `float64` instead of the pandas `str`/`NaN` the builder actually produces). The test instead asserts directly against `AnnotationTable.schema` and `AnnotationTable.to_pandas()`, so it runs Hail-free under the default (non-`-m hail`) pytest selection — consistent with § 3's reasoning for choosing the pandas backend in the first place.

Regenerate after an intentional builder or parser change:

```bash
pytest hvantk/skills/peptideatlas/phospho/tests/test_builder.py --regenerate-snapshots
pytest hvantk/skills/peptideatlas/phospho/tests -q
```
