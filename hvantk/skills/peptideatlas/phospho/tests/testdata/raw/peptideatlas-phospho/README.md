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
  `modified_peptide_sequence` in that build. A full-table scan of
  `modified_peptide_instance.tsv` across all 3,097,865 rows of the real
  build 606 archive (run as an HPC batch job, not on a login node) found
  exactly three `[STY]\[...\]` bracket shapes — `S[Phospho]` (2,420,951),
  `T[Phospho]` (558,273), `Y[Phospho]` (111,329) — and zero non-phospho
  modification or numeric-mass notation (such as `[167]`/`[181]`/`[243]`) on
  any S, T, or Y residue anywhere in the real build. The `\N` NULL sentinel
  and the trailing extra tab at the end of every header and data row are
  also real. See the plugin's `SKILL.md` s 4 for the full list of format
  facts.
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
| `DECOY_FAKE1` (`biosequence_id` 103, `protein_identification.presence_level_id` = `1` — canonical), `peptide_instance` 207 -> `modified_peptide_instance` 307 | dropped by the `DECOY_` accession-prefix filter ONLY. The canonical `presence_level_id` means the canonical-ID check alone would NOT remove this row, so a mutation that deletes the prefix filter now surfaces this site and fails the test |
| `CONTAM_FAKE1` (`biosequence_id` 104, `protein_identification.presence_level_id` = `1` — canonical), `peptide_instance` 208 -> `modified_peptide_instance` 308 | dropped by the `CONTAM_` accession-prefix filter ONLY, for the same reason as `DECOY_FAKE1` above |
| `P04637`, `peptide_instance` 209 (`n_observations` 100) -> `modified_peptide_instance` 309 + 310, the same `AAS[Phospho]AAAA` at charge 2 (12 observations) and charge 3 (8), mapped at `start_in_biosequence` 40 | one peptide in two phospho forms; the site counts each form's own observations, 12 + 8 = 20, not the peptide's 100 once per form (#425) |

Together `P04637` ends up with phospho sites on S, T, *and* Y residues,
covering all three phosphorylatable amino acids on one canonical protein.

## Expected `parse_peptideatlas_zip` output (6 rows)

| accession | position | amino_acid | n_observations |
|---|---|---|---|
| P04637 | 14 | T | 59 |
| P04637 | 32 | S | 30 |
| P04637 | 35 | Y | 30 |
| P04637 | 42 | S | 20 |
| P04637 | 53 | S | 8 |
| P99999 | 23 | S | 8 |

`DECOY_FAKE1`, `CONTAM_FAKE1`, `Q99999` and the `Carbamidomethyl`-only
peptide contribute zero rows.

## Not covered: non-phospho modifications on S, T or Y

A bug that counted any bracketed modification on S, T or Y as phospho would pass these
tests: the only non-phospho modification here is `C[Carbamidomethyl]`. The full scan
above is why: real build 606 carries no non-phospho modification on any S, T or Y
residue, so a row testing it would have to invent a notation, and data from this build
cannot trigger the bug.

## Recipe

Built by a one-off local script (not committed — the row design above and
in `test_builder.py`'s module docstring is the source of truth for
regenerating it): write the five TSVs in memory using the real headers from
this file's "What's real" section, with every row (header and data) ending
in the extra trailing tab the real dump also carries, then
`zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED).writestr(table_name,
content)` once per table, named `atlas_build_606-synthetic.tsv.zip`. Every
`biosequence_id` referenced by a peptide mapping — including the `DECOY_`/
`CONTAM_` rows — must also get a `protein_identification` row (see the
row-purpose table above): leaving one out lets the canonical-ID check drop
it, masking whichever other filter the row is meant to test.

The `peptide_instance` 209 rows (#425) were appended to the committed zip by copying an
existing row of each table and changing only the fields the row-purpose table names;
member order, compression and timestamps were kept.
