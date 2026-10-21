---
name: hvantk:resource-uniprot-ptm
description: Onboard, build, or update the UniProt PTM sites resource for hvantk
status: provisional
backend: hail
domain: proteomics
---

# UniProt PTM sites resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

This skill covers BUILD and UPDATE of the UniProt PTM sites Hail Table. It does NOT cover the internals of the PTM coordinate mapper (`hvantk/algorithms/ptm/pipeline.py`, `hvantk/algorithms/ptm/mapper.py`), which this plugin's parse stage calls, nor variant annotation against PTM regions (`hvantk/algorithms/ptm/annotate.py`) or the constraint LMM (`hvantk/algorithms/ptm/lmm.py`) — PTM analysis tooling that *consumes* the table built here.

## 2. Source identity

UniProt REST API endpoint, query, fields, and batch size live in `hvantk/skills/uniprot_ptm/shared/constants.py` (`UNIPROT_API_URL`, `UNIPROT_HUMAN_PTM_QUERY`, `UNIPROT_API_FIELDS`, `UNIPROT_BATCH_SIZE`). Read those; do not restate.

Stable provider notes:
- UniProt provides a live JSON search endpoint at `https://rest.uniprot.org/uniprotkb/search`. There is no archival versioning at the query level — `Last-Modified` may be absent. The version date on the downloaded TSV is the date of download.
- The drift probe (`hvantk/skills/uniprot_ptm/drift_probe.py`, probe version 2) fingerprints `X-UniProt-Release` (the release tag, e.g. `2026_03`) as `source_version`, `X-Total-Results` (the result count of the exact PTM query) as `extras.total_results`, and the sorted JSON key list of the first result as `headers.uniprot_entry_keys`; `X-UniProt-Release-Date` is recorded under `informational` and never signals drift. It fails closed — raises `DriftProbeError` rather than recording a partial fingerprint — when the release or count header is absent/unusable, the query returns zero results, or the HTTP request fails, because the scheduled drift bot regenerates a drifted baseline unattended and one transient omission would otherwise be committed as the new baseline (#271). What it cannot see: an in-place revision to a PTM annotation within the same release that leaves the result count unchanged — the count is a content proxy, not a full-body hash.
- The downloader (`UniProtPTMDataset.download`) pages by cursor (HTTP `Link: <...>; rel="next"`) and writes one TSV row per MOD_RES feature. It does NOT produce mapped genomic coordinates.
- The builder consumes the *mapped* TSV produced by `hvantk.algorithms.ptm.pipeline.map_ptm_sites`, not the raw UniProt download. The mapping step joins the UniProt TSV against an Ensembl GTF (and optionally CPTAC/PeptideAtlas phospho sources) to derive `chrom`, `codon_start`, `codon_end`, `strand`, `residue_pos`, `tissue_type`, etc. The plugin's `lifecycle.parse` runs that step for `hvantk reprocess`, on UniProt sites only and against the pinned Ensembl GRCh38 GTF (`hvantk/resources/ensembl_release.py`); `hvantk ptm build` runs the same mapper and can add the CPTAC/PeptideAtlas sites.

## 3. Backend choice + reasoning

`hail`. The mapped PTM TSV is small (~250k rows across all sources) but every downstream consumer joins it locus-wise against gnomAD-scale variant tables, so the table is materialised as a Hail Table keyed by `locus` to keep those joins cheap.

## 4. Raw format & gotchas

- File: tab-separated, single header row. Columns: `chrom`, `codon_start`, `codon_end`, `strand`, `uniprot_id`, `gene_symbol`, `residue_pos`, `amino_acid`, `ptm_type`, `ptm_category`, `source_db`, `evidence_type`, `n_observations`, `tissue_type` (see `PTM_OUTPUT_COLUMNS` in `hvantk/algorithms/ptm/constants.py`).
- The UniProt raw TSV (`uniprot-ptm-human-<YYYY-MM-DD>.tsv`) and the mapped TSV are different files. The raw TSV has columns `accession`, `gene_symbol`, `position`, `description`, `amino_acid`, `ensembl_xrefs`, `sequence_length` (`_TSV_COLUMNS` in `hvantk/skills/uniprot_ptm/shared/datasets.py`) and feeds the mapping step.
- Under `hvantk reprocess`, `--raw-dir` holds the dated raw TSVs, the pinned GTF (`Homo_sapiens.GRCh38.<release>.gtf.gz`, downloaded once and reused; delete it to fetch it again) and, unless `--intermediate` says otherwise, the mapped TSV (`uniprot-ptm_sites.intermediate`). The parse stage maps the newest dated raw TSV when there are several, writes block-gzip when the intermediate's name ends in `.bgz` and plain text otherwise (Hail picks the codec from the name), refuses a `.gz` name before doing any work (the builder's `hl.import_table` will not read a `.gz` file in parallel without `force_bgz`), and fails without writing anything if no site maps. It rejects any `reference_genome` but `GRCh38`, the build of the GTF.
- Import: `hl.import_table(impute=False, min_partitions=16, types={"codon_start": tstr, "codon_end": tstr, "residue_pos": tstr, "n_observations": tstr})`. Numeric fields stay as strings on import and are cast inside the builder so that empty cells do not abort `hl.import_table`.
- Contig remap: `MT → M`, then `chr` is **unconditionally** prepended (`hl.str("chr") + contig_remap.get(ht.chrom, ht.chrom)`) — there is no check for an already-`chr`-prefixed value. The mapped TSV's `chrom` column is always bare today (`17`, `X`, `MT`, ...; confirmed in the fixture), so this produces `chr17`/`chrX`/`chrM`; if `chrom` ever arrived pre-prefixed, this would double-prefix to `chrchr17` and the row would be silently dropped by the contig-validity filter below rather than erroring. Rows whose remapped contig is not in the reference genome contig list are filtered out (this is also how the non-reference `GL000009.2` fixture row is dropped — see § 9).
- `tissue_type` is provenance-bearing: CPTAC contributes `"normal"`/`"tumor"`, curated UniProt contributes `""`. The field is kept as-is by the builder so downstream stratification works.
- Keying: `key_by="locus"`. A `flanking_interval` field (`codon_start - flanking_codons*3` to `codon_end + flanking_codons*3`, clamped to chromosome bounds) is added. No hvantk code reads it today: `hvantk ptm annotate` builds its own window from its `--flanking-codons` option.
- `flanking_codons` defaults to 5; the builder raises `ValueError` if a negative value is passed.

## 5. Output contract

Hail Table at `<output_path>.ht`. Keyed by `locus`.

Row schema includes: `locus`, `chrom`, `codon_start`, `codon_end`, `strand`, `uniprot_id`, `gene_symbol`, `residue_pos`, `amino_acid`, `ptm_type`, `ptm_category`, `source_db`, `evidence_type`, `n_observations`, `tissue_type`, `flanking_interval`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/uniprot_ptm/plugin.yaml` (drives loader registration; compound dataset key `uniprot-ptm:sites`). The loader auto-resolves the dataset via `get_registry().get_dataset("uniprot-ptm:sites")`; top-level builds run through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).
- Builder: `build_uniprot_ptm_sites` in `hvantk/skills/uniprot_ptm/builder.py`. Signature `build_uniprot_ptm_sites(parsed_input, ctx, *, reference_genome="GRCh38", flanking_codons=5, fields=None, overwrite=False)` returning an `AnnotationTable`; `overwrite` is the download stage's argument and is ignored (the platform writes the output). It builds the table inline (`hl.import_table` + contig remap + locus/flanking-interval transforms); there is no `_create_table_base` helper and no `output_path` kwarg.
- Downloader CLI: `download_cmd` (Click `uniprot-ptm-download`) in `hvantk/skills/uniprot_ptm/cli.py`; lifecycle entry-point `download_dataset(raw_dir=..., overwrite=...)` (wired via `plugin.yaml` `lifecycle.download` and the `cli:` block). The umbrella `hvantk download` group (`hvantk/tools/plugins/download_cli.py`) strips the `-download` suffix, exposing it as `hvantk download uniprot-ptm`.
- Dataset class: `UniProtPTMDataset` in `hvantk/skills/uniprot_ptm/shared/datasets.py` (`from_latest()`, `from_date()`, `download()`).
- Build CLI: `hvantk reprocess uniprot-ptm:sites --raw-dir <dir> --output <path>.ht` runs download, parse (GTF + coordinate mapping) and build. Builder kwargs (`flanking_codons`, `fields`) flow through `--plugin-arg key=value`; `reference_genome` must stay `GRCh38`; `overwrite=true` downloads today's UniProt TSV again. `reprocess` passes every `--plugin-arg` to all three stages, so each stage names all four of these arguments and ignores the ones that are not its own; a misspelt one still fails.
- Parse stage: `map_uniprot_raw_dir` in `hvantk/algorithms/ptm/pipeline.py`, declared as `lifecycle.parse` in `plugin.yaml`. It reuses `download_ensembl_gtf`, `parse_ensembl_gtf` and `map_ptm_sites`, and lives in `algorithms/` because no module under `skills/` may import that layer (`hvantk/tests/test_dependency_directions.py`); the manifest names it instead.
- Multi-source build: `hvantk ptm build` (`ptm_build_pipeline` in `hvantk/tools/ptm/pipeline.py`) downloads UniProt itself, maps it with the same mapper, can add PeptideAtlas and CPTAC sites, and runs this plugin's builder on the combined TSV.
- Constants: `UNIPROT_API_URL`, `UNIPROT_HUMAN_PTM_QUERY`, `UNIPROT_API_FIELDS`, `UNIPROT_BATCH_SIZE` in `hvantk/skills/uniprot_ptm/shared/constants.py`; `PTM_OUTPUT_COLUMNS` in `hvantk/algorithms/ptm/constants.py`.
- Tests: `hvantk/skills/uniprot_ptm/tests/test_builder.py` (`test_uniprot_ptm_snapshot_round_trip`, the `@pytest.mark.hail @pytest.mark.round_trip("uniprot-ptm:sites")` test `plugin.yaml`'s `tests.command` selects — see § 9) and `test_drift_probe.py`. Broader mapping + builder integration is additionally exercised via `hvantk/tests/test_ptm.py`.

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Build end to end: `hvantk reprocess uniprot-ptm:sites --raw-dir <dir> --output <path>.ht`, adding `--plugin-arg flanking_codons=...` to change the stored `flanking_interval` (not the window `hvantk ptm annotate` uses). It downloads today's UniProt TSV and the pinned Ensembl GTF into `<dir>`, maps the sites (`lifecycle.parse`, see § 6) and runs the builder, which returns an `AnnotationTable` rather than writing to an `output_path`.
3. To build from files already in `<dir>` (offline, or a dated download kept for reproducibility), add `--skip-download`; the parse stage maps the newest `uniprot-ptm-human-<YYYY-MM-DD>.tsv` there and fetches the GTF only if it is missing. For a table that also carries PeptideAtlas or CPTAC sites, use `hvantk ptm build` instead (§ 6).
4. Reuse a mapped TSV without mapping again: `hvantk reprocess uniprot-ptm:sites --skip-download --skip-parse --raw-dir <dir> --intermediate <mapped.tsv> --output <path>.ht`.
5. Sanity-check the output: row count plausible (55,549 UniProt sites from the October 2026 release via `reprocess`; hundreds of thousands once `ptm build` adds the other sources); `locus` populated; `flanking_interval` length is roughly `(codon_end - codon_start) + 2 * flanking_codons * 3`.
6. Run validation: `pytest hvantk/skills/uniprot_ptm/tests -m hail`.

## 8. Update playbook

When UniProt releases a new snapshot (roughly every eight weeks):

1. Re-download: `hvantk download uniprot-ptm --output-dir <raw_dir> --overwrite`. Note the date label.
2. Diff the new raw TSV header against the previous fixture. New `_TSV_COLUMNS` entries imply a downloader schema change and require updating `shared/datasets.py`.
3. Regenerate the drift fingerprint with `hvantk drift --regenerate uniprot-ptm:sites --accepted-in "PR #N"`, never by writing `tests/drift_fingerprint.json` by hand: the command refuses a failing or placeholder-shaped probe result, replaces the file atomically, and records the dataset's row in the rebuild ledger (`hvantk/resources/drift_ledger.json`) when upstream changed (none for a `probe_version` change; `_conventions` § 12); commit both files. Inspect `source_version` and `extras.total_results` for the new release tag and result count, and `headers.uniprot_entry_keys` for additions or removals.
4. Re-run the mapping pipeline (`hvantk.algorithms.ptm.pipeline.map_ptm_sites`) to refresh the mapped TSV input.
5. Re-run `pytest hvantk/skills/uniprot_ptm/tests -m hail` and the broader PTM test suite (`pytest hvantk/tests/test_ptm.py -m hail`).
6. Open PR; reviewer checks the drift fingerprint diff and the row-count delta in the PR description.

## 9. Validation contract

- `fixture`: `hvantk/skills/uniprot_ptm/tests/testdata/raw/uniprot-ptm/ptm_sites.tsv` — synthetic, not a real slice (#426): hand-authored to carry the same content `test_plugin_conformance.py`'s round-trip builds inline (#213), using real gene symbols and UniProt accessions (TP53/`P04637`, BRCA2/`P51587`, NEMO-IKBKG/`Q9Y6K9`, MT-ND1/`P03886`) but fabricated residue positions, amino acids and genomic coordinates that do not match the real UniProt features at those accessions.
- `schema_snapshot`: `hvantk/skills/uniprot_ptm/tests/snapshots/schema.json`
- `row_snapshot`: `hvantk/skills/uniprot_ptm/tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `hvantk/skills/uniprot_ptm/tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/uniprot_ptm/tests -m hail`

> **Snapshot status:** seeded. `test_builder.py`'s `test_uniprot_ptm_snapshot_round_trip`
> drives the builder through `phase_b_snapshot_adapter` and asserts both snapshots.
> Regenerate after an intentional schema change with
> `pytest hvantk/skills/uniprot_ptm/tests/test_builder.py -m hail --regenerate-snapshots`,
> then commit the result. Broader integration coverage lives in `hvantk/tests/test_ptm.py`.

The fixture is synthetic: all 6 rows are fabricated (real gene symbols/accessions, invented coordinates and modifications), not a real UniProt excerpt. No `THIRD_PARTY_DATA.md` entry applies.
