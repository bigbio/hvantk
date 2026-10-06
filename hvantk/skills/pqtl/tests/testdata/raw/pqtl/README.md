# Synthetic pQTL fixture

SYNTHETIC. Contains no rows from Fang et al. (medRxiv 10.1101/2025.01.10.25320181,
licence `cc_no`). Header, delimiter and value conventions checked against
`/dss/work/heto4575/hvantk-datasets/raw/qtls/gtex/raw/Liver.allpairs_nobsGE72.txt.gz`
(whole file, 34,480,083 data rows) and
`.../Heart.allpairs_nobsGE72.txt.gz` (first 5,000,000 lines) on 2026-10-06.

## Format facts verified against the real files (both clean, zero anomalies)

- Header is exactly `gene_name SNP CHR BP A1 NMISS BETA STAT P`.
- Single-space delimited; every row has exactly 9 fields (no padding, no leading
  spaces, no double-spaces) in both files checked.
- `CHR` has no `chr` prefix (`17`, `13`, `X`, ...); the SNP id's contig does, and the
  two never disagreed in either file (0 mismatches across ~39.5M rows checked).
- `A1` always equals the ALT allele embedded in the `SNP` id (4th underscore field).
- `NMISS` is an integer.
- `BETA`, `STAT`, `P` are always numeric in the real files -- zero `NA` or other
  non-numeric values found in either file, so **no `NA` row is included here** (would
  misrepresent the format).
- `STAT == 0` never occurs in either real file, but the builder still filters it
  (`ht_part.filter(ht_part.STAT != 0.0)`), so this fixture fabricates one such row to
  exercise that path.
- `gene_name` is never empty and never looks like an Ensembl id (`^ENSG...`) in either
  file checked (0 of ~39.5M rows) -- so **no Ensembl-id `gene_name` row is included
  here** either, for the same reason.
- The only `SNP` id suffix observed is `_b38` (100% of rows in both files).

## Rows

| Row | gene_name | SNP | Purpose |
|---|---|---|---|
| 1 | `BRCA1` | `chr17_43094687_A_G_b38` | mapped symbol, positive BETA, valid GRCh38 chr17 position |
| 2 | `BRCA2` | `chr13_32340073_G_C_b38` | mapped symbol, negative BETA/STAT -- SE must still come out positive |
| 3 | `TP53` | `chr17_43094687_A_G_b38` | second gene at the *same* locus/alleles as row 1, different `gene_id` |
| 4 | `BRCA1` | `chr17_43095211_A_G_b38` | `STAT == 0` -- must be dropped by the builder |
| 5 | `SYNTHGENE1` | `chr5_100000_A_T_b38` | unmapped symbol -- `gene_id` falls back to the raw symbol |
| 6 | `BRCA2` | `chrX_71130000_C_T_b38` | chrX variant -- non-numeric contig handling |
| 7 | `TP53` | `chr1_1000000_AT_A_b38` | indel id (multi-base REF `AT` / ALT `A`) |
| 8 | `BRCA1` | `chr17_43110000_T_C_b38` | `P = 1e-300` -- extreme float, underflow-adjacent |

The two conditional rows described in the fixture-authoring task ("an `NA` row" and "a
`gene_name` that is an Ensembl id") were deliberately **omitted** because the Step 0
verification above found zero instances of either pattern in ~39.5M real rows across
two tissues -- fabricating them would assert a format behavior the real source does not
exhibit.

All associations (gene/variant pairings, BETA/STAT/P values) are fabricated. Genomic
positions are valid GRCh38 coordinates for their stated chromosome but were not checked
against real gene boundaries beyond rows 1 and 3 (chosen inside the BRCA1 locus,
chr17:43,044,295-43,125,364) and row 2 (inside BRCA2, chr13:32,315,474-32,400,266).

## Generation

1. Verified the real-format facts above by streaming `zcat <file> | awk ...` over the
   Liver file (whole file) and the first 5,000,000 lines of the Heart file on the
   cluster (`/dss/work/heto4575/agent-runs/pqtl-step0/`), via `sbatch` on
   `rosa_express.p` (job 20318211, ~1m43s elapsed).
2. Hand-authored 8 rows (above) as plain text with the verified header and single-space
   delimiter.
3. `gzip -c Liver.allpairs_nobsGE72-synthetic.txt > Liver.allpairs_nobsGE72-synthetic.txt.gz`
   (plain gzip, not bgzip -- matches the real files and what
   `hvantk.core.utils.qtl_helpers.scan_tissue_files` / `hl.import_table(force=True)`
   expect).

`scan_tissue_files` only globs `*.txt.gz` and `*.tsv.gz`, so this README is never picked
up as a data file.
