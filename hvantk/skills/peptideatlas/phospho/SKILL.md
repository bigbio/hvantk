---
name: hvantk:resource-peptideatlas-phospho
description: Download and parse a PeptideAtlas Human Phospho build into a wide intermediate TSV consumed by the hvantk PTM pipeline.
status: provisional
backend: pandas
domain: proteomics
---

# PeptideAtlas Human Phospho resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. Downloader + parser + downloader tests + drift-probe placeholder are migrated under the plugin folder. Builder round-trip snapshots (`schema.json`, `sample_rows.json`) are NOT yet seeded; the fixture directory exists but is empty.
- **In scope:** one PeptideAtlas human phospho build at a time (default: `(202512, 606)` — see `PEPTIDEATLAS_LATEST_BUILD_DATE` / `PEPTIDEATLAS_LATEST_BUILD_ID` in `hvantk/skills/peptideatlas/phospho/shared/constants.py`). Output is the wide intermediate TSV consumed by `hvantk/algorithms/ptm/pipeline.py` via its `peptideatlas_tsv` argument.
- **Out of scope:** non-human PeptideAtlas builds; non-phospho PTM atlases on PeptideAtlas (those would land as sibling datasets under `hvantk/skills/peptideatlas/<ptm-type>/`); any Hail-Table or AnnData representation — the only downstream consumer today reads the TSV directly.

## 2. Source identity

- **Provider:** PeptideAtlas (<https://peptideatlas.org/>).
- **Catalog entry:** TODO. No PeptideAtlas entry exists yet in any plugin's `catalog/datasets.json` (verify with `hvantk catalog search peptideatlas`). The pinned build coordinates live in `hvantk/skills/peptideatlas/phospho/shared/constants.py` as `PEPTIDEATLAS_PHOSPHO_BASE_URL`, `PEPTIDEATLAS_LATEST_BUILD_DATE`, and `PEPTIDEATLAS_LATEST_BUILD_ID`. When a dedicated catalog entry lands under `hvantk/skills/peptideatlas/catalog/datasets.json`, point this plugin's `source.catalog_ref` at it and remove this TODO.
- Download URL is composed by `PeptideAtlasPhosphoDataset.from_build(build_date, build_id)` as `{base}/{build_date}/atlas_build_{build_id}.tsv.zip`. The dataset class enforces HTTPS + `*.peptideatlas.org` host validation before any GET.

## 3. Backend choice + reasoning

**`backend: pandas`, `domain: proteomics`.** Per `_conventions` § 3, Hail Tables / MatrixTables are for variant- or expression-keyed data joined into the wider hvantk genomics graph; AnnData is for sparse single-cell or sample x gene matrices. PeptideAtlas phospho output is neither: it is a narrow site-level table (one row per `(accession, position)`) keyed by UniProt accession and is consumed downstream by the PTM pipeline (`hvantk/algorithms/ptm/pipeline.py`) as a TSV file. A Hail Table would force a Spark session for a sub-100 MB lookup table. A pandas DataFrame is the natural fit; the Phase B builder `build_peptideatlas_phospho(parsed_input, ctx, **params)` loads the intermediate TSV via `pd.read_csv(... sep="\t" ...)` and wraps it as an `AnnotationTable`.

## 4. Raw format & gotchas

The raw upstream artefact is a single TSV zip (e.g. `atlas_build_606.tsv.zip`, multiple GB) containing PeptideAtlas's relational dump:

- `biosequence.tsv` — proteins (`biosequence_id`, accession, gene name, sequence).
- `protein_identification.tsv` — protein confidence; `presence_level_id == 1` means *canonical*. The parser filters by this and drops `DECOY_` / `CONTAM_` prefixed accessions in `biosequence.tsv`.
- `peptide_instance.tsv` — distinct peptides + `n_observations` (sum of PSMs across experiments).
- `peptide_mapping.tsv` — peptide-to-protein coordinates (`start_in_biosequence`); 15M+ rows; must be streamed, not slurped.
- `modified_peptide_instance.tsv` — modified peptide sequences in bracket notation; 3M+ rows.

Modification notation gotchas (see `_extract_phospho_offsets` in `shared/datasets.py`):

- Text form: `S[Phospho]`, `T[Phospho]`, `Y[Phospho]`, or `S[Phospho:80]` (substring match handles the latter).
- Numeric form: `S[167]`, `T[181]`, `Y[243]` — the total modified-residue mass (residue + HPO3 ≈ 80). The parser only treats numeric brackets as phospho when the preceding residue is S/T/Y and the mass matches within ±1.0 Da.
- N-terminal labels like `[TMT6plex]-` precede the first residue and must be skipped (no preceding amino acid).
- Lowercase letters, digits, and dashes in the sequence are ignored.

Aggregation gotcha: the same `(accession, position)` can be observed via multiple distinct peptides — `parse_peptideatlas_zip` sums `n_observations` across all contributing peptide instances.

## 5. Output contract

- **File:** wide TSV at the path returned by `download_dataset` / `parse_raw_dir`, named `peptideatlas-phospho-<build_date>-<build_id>.tsv` (e.g. `peptideatlas-phospho-202512-606.tsv`).
- **Shape:** one row per `(accession, position)` phospho site on a canonical protein.
- **Columns** (fixed order, see `_TSV_COLUMNS`): `accession`, `gene_symbol`, `position`, `description`, `amino_acid`, `ensembl_xrefs`, `sequence_length`, `n_observations`, `source_db`, `evidence_type`. `source_db` is always `"PeptideAtlas"`; `evidence_type` is always `"mass_spectrometry"`; `ensembl_xrefs` is currently always empty.
- **Builder return value:** `build_peptideatlas_phospho` returns an `AnnotationTable` (via `AnnotationTable.from_pandas(...)`) wrapping `pd.read_csv(<tsv>, sep="\t", dtype=str)` — all columns as strings so downstream consumers don't get silent numeric coercion on `position` / `n_observations`. Provenance is stamped from `ctx.provenance(schema_id="peptideatlas-phospho-v1")`.

## 6. hvantk integration points

- **Dataset class + parser:** `PeptideAtlasPhosphoDataset`, `parse_peptideatlas_zip`, `write_intermediate_tsv`, `parse_raw_dir` in `hvantk/skills/peptideatlas/phospho/shared/datasets.py`.
- **Builder:** `build_peptideatlas_phospho(parsed_input, ctx, **params) -> AnnotationTable` in `hvantk/skills/peptideatlas/phospho/builder.py` (declared in `plugin.yaml` under `datasets[].builder`). A legacy `build_peptideatlas_phospho_tb(input_path, output_path, ...)` helper also lives in that module but is NOT the loader entry point.
- **Downloader CLI:** `download_cmd` in `hvantk/skills/peptideatlas/phospho/cli.py` (registered as `hvantk peptideatlas-phospho-download` and re-bound under `hvantk download peptideatlas-phospho` via `hvantk/tools/plugins/download_cli.py`).
- **Lifecycle entry points:** `download_dataset` and `parse_raw_dir` (loader-wired via `lifecycle.download` + `lifecycle.parse` in `plugin.yaml`).
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/peptideatlas/phospho/drift_probe.py` (HEAD against the pinned build's zip URL).
- **Downstream consumer:** `hvantk/algorithms/ptm/pipeline.py` (`PTMBuildConfig.peptideatlas_tsv`) — reads the intermediate TSV produced here and maps PTM sites to genomic coordinates. Exposed at the user-facing level by `hvantk/algorithms/ptm/atlas.py` (`PTMAtlasConfig.peptideatlas_tsv`).
- **Plugin manifest:** `hvantk/skills/peptideatlas/plugin.yaml` (compound dataset key `peptideatlas:phospho`).
- **Tests:** `hvantk/skills/peptideatlas/phospho/tests/` (parser unit tests + drift-probe sanity test).

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update the PeptideAtlas phospho intermediate:

1. **End-to-end (recommended).** Run the full download -> parse -> build chain through the plugin loader:

   ```bash
   hvantk reprocess peptideatlas:phospho --raw-dir /data/peptideatlas --output /out/peptideatlas-phospho.ht
   ```

   The loader auto-resolves the dataset from `plugin.yaml` (`get_registry().get_dataset("peptideatlas:phospho")`) and runs the build through `run_builder_for_spec`. The `lifecycle.download` (`download_dataset`) and `lifecycle.parse` (`parse_raw_dir`) entry points run first, then the Phase B `build_peptideatlas_phospho`.
2. **Download only.** Either via the standalone CLI (`hvantk peptideatlas-phospho-download -o /data/peptideatlas`) or the lifecycle entry point `download_dataset(raw_dir=...)`. Both produce `<raw_dir>/atlas_build_<id>.tsv.zip` *and* the parsed `<raw_dir>/peptideatlas-phospho-<date>-<id>.tsv`.
3. **(Lifecycle) parse-only step.** `parse_raw_dir(raw_dir=..., output_path=...)` re-parses an existing zip from `raw_dir` into a fresh intermediate TSV — used when downstream code wants the TSV at a different path than the dataset class's default.
4. **Builder.** `build_peptideatlas_phospho(parsed_input, ctx, **params)` loads the intermediate TSV (the path produced by `parse_raw_dir`) as a pandas DataFrame and returns an `AnnotationTable`.
5. **Validate.** `pytest hvantk/skills/peptideatlas/phospho/tests` — parser unit tests + drift-probe sanity. Builder round-trip snapshot test is TODO (see § 9).

## 8. Update playbook

PeptideAtlas releases a new human phospho build once or twice per year. When a new build is announced:

1. Run the drift probe: `python -c "from hvantk.skills.peptideatlas.phospho.drift_probe import fetch_fingerprint; print(fetch_fingerprint())"`. A change in `source_version` (the `Last-Modified` header) is the trigger.
2. Update the pinned build coordinates in `hvantk/skills/peptideatlas/phospho/shared/constants.py` (`PEPTIDEATLAS_LATEST_BUILD_DATE`, `PEPTIDEATLAS_LATEST_BUILD_ID`).
3. Re-download (`hvantk peptideatlas-phospho-download -o /data/peptideatlas --overwrite`) and spot-check the row count against the previous build.
4. If any new modification notation or table appears in the dump, document it in § 4.
5. Re-regenerate `tests/drift_fingerprint.json` with the new build's filename + headers.
6. If a builder snapshot exists (TODO), run `pytest hvantk/skills/peptideatlas/phospho/tests --regenerate-snapshots` and inspect the diff.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/peptideatlas/phospho/tests/testdata/raw/peptideatlas-phospho/` (directory present, not yet seeded — the existing parser unit tests build mock zips on-the-fly via the `mock_pa_zip` fixture in `test_phospho.py`). A future round-trip test will seed a slim multi-table zip here.
- **schema_snapshot:** `hvantk/skills/peptideatlas/phospho/tests/snapshots/schema.json` (TODO — created on first `--regenerate-snapshots` run).
- **row_snapshot:** `hvantk/skills/peptideatlas/phospho/tests/snapshots/sample_rows.json` (TODO — same).
- **test_command:** `pytest hvantk/skills/peptideatlas/phospho/tests`.
- **drift_fingerprint:** `hvantk/skills/peptideatlas/phospho/tests/drift_fingerprint.json` (placeholder shape; refresh via the update playbook).

The plugin manifest already declares these paths so the loader contract holds. Parser unit tests + drift-probe sanity test pass today; the builder round-trip is the gap to close in a follow-up.

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. This is a pure-pandas builder (no Hail/Spark required). On the
> first round-trip run, use
> `pytest hvantk/skills/peptideatlas/phospho/tests/test_builder.py --regenerate-snapshots`
> (`test_builder.py` to be added) to bootstrap them, then commit. Until seeded, the
> round-trip test cannot verify output against a fixed schema.
