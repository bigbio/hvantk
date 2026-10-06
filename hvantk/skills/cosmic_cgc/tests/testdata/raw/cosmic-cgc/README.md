# cosmic-cgc synthetic fixture

SYNTHETIC. Contains no COSMIC data. Header and value conventions checked against
Cosmic_CancerGeneCensus_v103_GRCh38 (licensed copy) on 2026-10-06.

The COSMIC Cancer Gene Census licence (amended 30 Oct 2025, cl. 2.1.2, 4.7, C.h)
forbids redistributing COSMIC data, including individual rows or gene-level
annotations. `cosmic-cgc-synthetic.tsv.gz` carries the real v103 column header
and the real delimiter/separator/flag conventions, but every value — gene
symbols, names, IDs, coordinates, syndromes — is fabricated. No row in this
file corresponds to a real COSMIC Cancer Gene Census entry.

## Format

- 21 tab-separated columns, exact v103 header order:
  `GENE_SYMBOL NAME COSMIC_GENE_ID CHROMOSOME GENOME_START GENOME_STOP
  CHR_BAND SOMATIC GERMLINE TUMOUR_TYPES_SOMATIC TUMOUR_TYPES_GERMLINE
  CANCER_SYNDROME TISSUE_TYPE MOLECULAR_GENETICS ROLE_IN_CANCER
  MUTATION_TYPES TRANSLOCATION_PARTNER OTHER_GERMLINE_MUT OTHER_SYNDROME
  TIER SYNONYMS`
- Compressed with plain gzip (Python's stdlib `gzip` module), matching the
  real export (`file` on the licensed `.tsv.gz` reports standard gzip, not
  BGZF — confirmed from the magic bytes: `FLG=0x08` (FNAME set), no BGZF
  `"BC"` extra-field subfield). This exercises the same `resolve_compression()`
  branch (`force_bgz=False`) as the production file.
- No quoting; no literal `NA`/`N/A`/`-`/`null` missing-value sentinels —
  missing values are the empty string.
- `TIER` is the bare digit `"1"` or `"2"`. `SOMATIC`/`GERMLINE` are `"y"`/`"n"`
  (never `"yes"`; both columns are always populated in the real v103 export).
- Multi-value fields (`TUMOUR_TYPES_SOMATIC`, `TUMOUR_TYPES_GERMLINE`,
  `ROLE_IN_CANCER`, `MUTATION_TYPES`) are comma-separated, usually `", "`
  (comma-space) but not always — the real export has rows with a bare `,`
  and no following space, which is why the builder splits on `,` and then
  strips each token rather than splitting on the literal `", "`.
- `ROLE_IN_CANCER` real vocabulary: `oncogene`, `TSG`, `fusion`.
  `TISSUE_TYPE` real vocabulary: `E`, `L`, `M`, `O`.
  `COSMIC_GENE_ID` real shape: `COSG` + 5-6 digits.
  `CHROMOSOME` real form: bare `1`-`22` or `X` (no `chr` prefix).
  `GENOME_START`/`GENOME_STOP` are numeric or empty (6/763 rows are empty
  in the licensed v103 export).

## Rows

| Row | Purpose |
|---|---|
| SYNTHA | `TIER` 1; `SOMATIC`/`GERMLINE` in the real `"y"` form; `TUMOUR_TYPES_SOMATIC` has a stray-space empty entry (`"A, , B"`) — exercises strip + empty-entry dropping |
| SYNTHB | `TIER` 2; `SOMATIC` empty; `GERMLINE` in the real `"y"` form — exercises tier normalisation and a false boolean from an empty cell |
| SYNTHC | every multi-value field empty — exercises `[]` arrays (not missing) |
| SYNTHD | `ROLE_IN_CANCER` with two roles, `MUTATION_TYPES` with two codes, joined with the real `", "` separator — exercises multi-value parsing with 2+ elements |
| SYNTHE | a gene symbol the stub `gene_catalog` in `test_builder.py` does not map — exercises dropping an unresolved row when keyed by `hgnc_id`; also carries empty `GENOME_START`/`GENOME_STOP` to exercise the missing-tolerant `hl.parse_int32` cast before the row is filtered out |

## Recipe

Generated with a small script (not committed) that wrote the header and the
five rows above as tab-separated text, then compressed with Python's
`gzip.open(path, "wt", encoding="utf-8", newline="\n")` — no COSMIC input
file was read to produce it. Regenerate by recreating that script from this
table if the fixture ever needs to change; do not derive it from a real
export.
