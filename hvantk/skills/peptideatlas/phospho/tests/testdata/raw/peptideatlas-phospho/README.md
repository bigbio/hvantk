# Synthetic PeptideAtlas raw-build fixture

SYNTHETIC. Contains no PeptideAtlas rows. Table headers and modification
notation checked against human phospho build 202512/606 on 2026-10-06.

`atlas_build_606-synthetic.tsv.zip` is a miniature, fully fabricated stand-in for
a PeptideAtlas `atlas_build_<id>.tsv.zip` raw TSV dump. It lets
`parse_raw_dir` / `parse_peptideatlas_zip`
(`hvantk/skills/peptideatlas/phospho/shared/datasets.py`) run against
something real-shaped — real table names, real column headers in real order,
real modification notation — end to end, rather than only the
already-parsed intermediate TSV the builder itself consumes.

## What's real, what's fabricated

- **Real:** the five table names, every column name, and column order in
  each, captured verbatim via `unzip -p atlas_build_606.tsv.zip <table> |
  head -3` on the cluster login node. The phospho modification notation
  (`<residue>[Phospho]`) and a non-phospho example
  (`<residue>[Carbamidomethyl]`) are the real notations seen in
  `modified_peptide_sequence` in that build — confirmed as the *only* bracket
  shapes across a 1.55M-row / 300MB sample (zero numeric-mass notation such
  as `[167]`/`[181]`/`[243]`). The `\N` NULL sentinel and the trailing extra
  tab at the end of every header and data row are also real. See the
  plugin's `SKILL.md` s 4 for the full list of format facts.
- **Fabricated:** every accession (except the reused placeholder `P04637`,
  as the pre-existing `test_phospho.py` mocks already did), every sequence,
  every count, every ID. No value in this fixture was copied from a real
  PeptideAtlas row.

## Row-purpose table

| Rows | Purpose |
|---|---|
| `P04637` (fake seq), `peptide_instance` 201 + 202 -> `modified_peptide_instance` 301 + 302, mapped onto `P04637` at offsets that land on the same residue | two peptides covering the same site; counts must sum (42 + 17 = 59) |
| `P04637`, `peptide_instance` 203 -> `modified_peptide_instance` 303 (`AAS[Phospho]AAY[Phospho]AA`) | one peptide with two phospho sites (S and Y) — multi-site offset extraction |
| `P04637`, `peptide_instance` 204 -> `modified_peptide_instance` 304 (`AAC[Carbamidomethyl]AAAA`) | a non-phospho modification; must contribute zero output rows |
| `P04637` + `P99999`, `peptide_instance` 205 -> `modified_peptide_instance` 305, two `peptide_mapping` rows (`matched_biosequence_id` 100 and 101) | one peptide mapping to two proteins — multi-mapping |
| `Q99999` (`biosequence_id` 102), `protein_identification.presence_level_id` = `3`, `peptide_instance` 206 -> `modified_peptide_instance` 306 | non-canonical isoform — excluded by the canonical filter |
| `DECOY_FAKE1` (`biosequence_id` 103), `peptide_instance` 207 -> `modified_peptide_instance` 307 | dropped by the `DECOY_` accession-prefix filter |
| `CONTAM_FAKE1` (`biosequence_id` 104), `peptide_instance` 208 -> `modified_peptide_instance` 308 | dropped by the `CONTAM_` accession-prefix filter |

Together `P04637` ends up with phospho sites on S, T, *and* Y residues,
covering all three phosphorylatable amino acids on one canonical protein.

## Expected `parse_peptideatlas_zip` output (5 rows)

| accession | position | amino_acid | n_observations |
|---|---|---|---|
| P04637 | 14 | T | 59 |
| P04637 | 32 | S | 30 |
| P04637 | 35 | Y | 30 |
| P04637 | 53 | S | 8 |
| P99999 | 23 | S | 8 |

`DECOY_FAKE1`, `CONTAM_FAKE1`, `Q99999` and the `Carbamidomethyl`-only
peptide contribute zero rows.

## Recipe

Built by a one-off local script (not committed — the row design above and
in `test_builder.py`'s module docstring is the source of truth for
regenerating it): write the five TSVs in memory using the real headers from
this file's "What's real" section, with every row (header and data) ending
in the extra trailing tab the real dump also carries, then
`zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED).writestr(table_name,
content)` once per table, named `atlas_build_606-synthetic.tsv.zip`.
