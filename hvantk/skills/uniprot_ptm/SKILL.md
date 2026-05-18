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

This skill covers BUILD and UPDATE of the UniProt PTM sites Hail Table. It does NOT cover the upstream PTM coordinate mapping pipeline (`hvantk/ptm/pipeline.py`), variant annotation against PTM regions (`hvantk/ptm/annotate.py`), or the constraint LMM (`hvantk/ptm/analysis.py`, `hvantk/ptm/test.py`) — those modules are PTM analysis tooling that *consumes* the table built here.

## 2. Source identity

UniProt REST API endpoint, query, fields, and batch size live in `hvantk/ptm/constants.py` (`UNIPROT_API_URL`, `UNIPROT_HUMAN_PTM_QUERY`, `UNIPROT_API_FIELDS`, `UNIPROT_BATCH_SIZE`). Read those; do not restate.

Stable provider notes:
- UniProt provides a live JSON search endpoint at `https://rest.uniprot.org/uniprotkb/search`. There is no archival versioning at the query level — `Last-Modified` may be absent. The version date on the downloaded TSV is the date of download.
- The downloader (`UniProtPTMDataset.download`) pages by cursor (HTTP `Link: <...>; rel="next"`) and writes one TSV row per MOD_RES feature. It does NOT produce mapped genomic coordinates.
- The builder consumes the *mapped* TSV produced by `hvantk.ptm.pipeline.map_ptm_sites`, not the raw UniProt download. The mapping step joins the UniProt TSV against an Ensembl GTF (and optionally CPTAC/PeptideAtlas phospho sources) to derive `chrom`, `codon_start`, `codon_end`, `strand`, `residue_pos`, `tissue_type`, etc.

## 3. Backend choice + reasoning

`hail`. The mapped PTM TSV is small (~250k rows across all sources) but every downstream consumer joins it locus-wise against gnomAD-scale variant tables, so the table is materialised as a Hail Table keyed by `locus` to keep those joins cheap.

## 4. Raw format & gotchas

- File: tab-separated, single header row. Columns: `chrom`, `codon_start`, `codon_end`, `strand`, `uniprot_id`, `gene_symbol`, `residue_pos`, `amino_acid`, `ptm_type`, `ptm_category`, `source_db`, `evidence_type`, `n_observations`, `tissue_type` (see `PTM_OUTPUT_COLUMNS` in `hvantk/ptm/constants.py`).
- The UniProt raw TSV (`uniprot-ptm-human-<YYYY-MM-DD>.tsv`) and the mapped TSV are different files. The raw TSV has columns `accession`, `gene_symbol`, `position`, `description`, `amino_acid`, `ensembl_xrefs`, `sequence_length` (`_TSV_COLUMNS` in `shared/datasets.py`) and feeds the mapping step.
- Import: `hl.import_table(impute=False, min_partitions=16, types={"codon_start": tstr, "codon_end": tstr, "residue_pos": tstr, "n_observations": tstr})`. Numeric fields stay as strings on import and are cast inside the builder so that empty cells do not abort `hl.import_table`.
- Contig remap: `MT → chrM` then prepend `chr` to all non-`chr`-prefixed values. Rows whose remapped contig is not in the reference genome contig list are filtered out.
- `tissue_type` is provenance-bearing: CPTAC contributes `"normal"`/`"tumor"`, curated UniProt contributes `""`. The field is kept as-is by the builder so downstream stratification works.
- Keying: `key_by="locus"`. A `flanking_interval` field (`codon_start - flanking_codons*3` to `codon_end + flanking_codons*3`, clamped to chromosome bounds) is added for proximity joins.
- `flanking_codons` defaults to 5; the builder raises `ValueError` if a negative value is passed.

## 5. Output contract

Hail Table at `<output_path>.ht`. Keyed by `locus`.

Row schema includes: `locus`, `chrom`, `codon_start`, `codon_end`, `strand`, `uniprot_id`, `gene_symbol`, `residue_pos`, `amino_acid`, `ptm_type`, `ptm_category`, `source_db`, `evidence_type`, `n_observations`, `tissue_type`, `flanking_interval`.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/uniprot_ptm/plugin.yaml` (drives loader registration; compound dataset key `uniprot-ptm:sites`).
- Builder: `create_ptm_sites_tb` in `hvantk/skills/uniprot_ptm/builder.py` (uses `_create_table_base()` per `_conventions` § 4).
- Downloader CLI: `download_cmd` (Click `uniprot-ptm-download`) in `hvantk/skills/uniprot_ptm/cli.py`; lifecycle entry-point `download_dataset(raw_dir=...)`. Wired into the umbrella `hvantk download uniprot-ptm` group in `hvantk/tools/plugins/download_cli.py`.
- Dataset class: `UniProtPTMDataset` in `hvantk/skills/uniprot_ptm/shared/datasets.py`.
- Build CLI: `hvantk mktable ptm-sites` in `hvantk/tools/build/make_table_cli.py` (user-facing name preserved).
- PTM mapping pipeline (out-of-plugin, intentionally): `hvantk/ptm/pipeline.py` (`map_ptm_sites`, `download_uniprot_ptm`) — the analysis tooling that produces the *mapped* TSV consumed by the builder.
- Constants: `UNIPROT_API_URL`, `UNIPROT_HUMAN_PTM_QUERY`, `UNIPROT_API_FIELDS`, `UNIPROT_BATCH_SIZE`, `PTM_OUTPUT_COLUMNS` in `hvantk/ptm/constants.py`.
- Tests: `hvantk/skills/uniprot_ptm/tests/test_drift_probe.py`. The end-to-end builder is exercised indirectly via `hvantk/tests/test_ptm.py` (mapping + builder integration).

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Download the raw UniProt TSV: `hvantk download uniprot-ptm --output-dir <raw_dir>` (or call `UniProtPTMDataset.from_latest().download(raw_dir)`).
3. Run the mapping pipeline to attach genomic coordinates: see `hvantk.ptm.pipeline.map_ptm_sites`. Output is the mapped TSV consumed by step 4.
4. Build via Python (`create_ptm_sites_tb(input_path, output_path, reference_genome=..., flanking_codons=..., overwrite=True)`) or CLI (`hvantk mktable ptm-sites --raw-input ... --output-ht ...`).
5. Sanity-check the output: row count plausible (hundreds of thousands across sources); `locus` populated; `flanking_interval` length is roughly `(codon_end - codon_start) + 2 * flanking_codons * 3`.
6. Run validation: `pytest hvantk/skills/uniprot_ptm/tests -m hail`.

## 8. Update playbook

When UniProt releases a new monthly snapshot (typically late in the month):

1. Re-download: `hvantk download uniprot-ptm --output-dir <raw_dir> --overwrite`. Note the date label.
2. Diff the new raw TSV header against the previous fixture. New `_TSV_COLUMNS` entries imply a downloader schema change and require updating `shared/datasets.py`.
3. Regenerate the drift fingerprint: invoke `fetch_fingerprint()` against the live endpoint and overwrite `tests/drift_fingerprint.json`. Inspect the `uniprot_entry_keys` field for additions or removals.
4. Re-run the mapping pipeline (`hvantk.ptm.pipeline.map_ptm_sites`) to refresh the mapped TSV input.
5. Re-run `pytest hvantk/skills/uniprot_ptm/tests -m hail` and the broader PTM test suite (`pytest hvantk/tests/test_ptm.py -m hail`).
6. Open PR; reviewer checks the drift fingerprint diff and the row-count delta in the PR description.

## 9. Validation contract

- `fixture`: `hvantk/skills/uniprot_ptm/tests/testdata/raw/uniprot-ptm/` (placeholder; integration coverage lives in `hvantk/tests/test_ptm.py`)
- `drift_fingerprint`: `hvantk/skills/uniprot_ptm/tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/uniprot_ptm/tests -m hail`

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run in a hail-enabled environment, use
> `pytest hvantk/skills/uniprot_ptm/tests/test_builder.py --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
