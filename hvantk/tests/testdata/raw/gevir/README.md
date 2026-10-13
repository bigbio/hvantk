# Synthetic GeVIR fixture

SYNTHETIC. Contains no rows from the GeVIR paper (Abramovs, Brass & Tassabehji,
*Nature Genetics* 52(1):35-39, 2020, PMID 31873297, DOI 10.1038/s41588-019-0560-2).

## Why synthetic, not a real excerpt (issue #386)

The journal is subscription-access (`journal_is_oa: false` per Unpaywall); this
article's only open-access route is a green, self-archived copy in the University
of Manchester's institutional repository under an unspecified ("other-oa") licence
-- not a gold/CC-BY licence from the publisher. The Supplementary Table 2 workbook
is served from the publisher's own CDN (static-content.springer.com) with no CC
licence statement found, so redistributing real rows (even a handful) is not
clearly permitted. Checked 2026-10-13 via Unpaywall (`api.unpaywall.org`) and
PubMed's copyright-status lookup (no PMC deposit, no open-access licence on
record).

## Format facts verified against the real workbook (checked 2026-10-13)

The real workbook was downloaded directly from
`static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0560-2/MediaObjects/
41588_2019_560_MOESM3_ESM.xlsx` (10,270,511 bytes, matching the committed drift
fingerprint) and parsed with `hvantk/skills/gevir/cli.py`'s own extractor:

- Sheet `table_2`, header on row 2, data on rows 3-19363 (19,361 genes) -- confirmed
  by reading the live workbook end to end.
- Output header is exactly `gnomad_gene_name	gene_id	canonical_transcript
  gevir_pct	virlof_pct` (tab-separated), matching the pre-existing (pre-#386)
  derived TSV's header.
- `gevir_pct`/`virlof_pct` are plain decimal strings (no scientific notation seen
  in columns F/H of the live file), full IEEE-754 double precision as Excel stored
  them (e.g. `6.7455193429999998`) -- this fixture mirrors that precision
  convention rather than the pre-#386 11-character-rounded one.
- 29 of 19,361 `gnomad_gene_name` cells are Excel date serials, not strings (see
  `hvantk/skills/gevir/cli.py`'s module docstring for the five distinct gene
  families this covers, confirmed against the live workbook plus HGNC's and
  Ensembl's REST APIs). Row 26 below (`MARCH6`) exercises the restored-name shape
  for builder-input format fidelity; the restoration logic itself (decoding the
  Excel serial) is covered by `hvantk/skills/gevir/tests/test_downloader.py`
  against its own tiny synthetic *workbook*, not this TSV.

## Rows

26 rows total:

- Rows 1-6: the exact `gene_id` values `hvantk/skills/gevir/tests/test_builder.py`'s
  `SAMPLE_KEYS` pins (`ENSG00000000003`, `ENSG00000000005`, `ENSG00000000419`,
  `ENSG00000000457`, `ENSG00000000460`, `ENSG00000000938`) -- real, public Ensembl
  gene identifiers (not GeVIR's data), each paired with a fabricated symbol
  (`SYNG1`-`SYNG6`), fabricated transcript id, and fabricated percentiles. None of
  these values match the real, published GeVIR numbers for these genes.
- Rows 7-25: 19 more fabricated rows under an out-of-range, obviously-synthetic
  gene-id namespace (`ENSG00000900XXX`), giving the fixture a believable
  multi-row table shape without approaching the real 19,361-row size.
- Row 26: `gnomad_gene_name = MARCH6`, `gene_id = ENSG00000900099` -- an ordinary
  (already-restored) string value exercising the same shape a restored,
  date-mangled cell produces downstream of the real downloader, with fabricated
  percentiles.

All symbols, transcript ids, and percentile values are fabricated.

## Generation

1. Verified the real-format facts above against a live download of the workbook
   on 2026-10-13 (see above).
2. Generated 26 rows with a fixed random seed (20261013): the 6 required sample
   keys, 19 filler rows, and 1 format-fidelity row (script not committed -- a
   one-off, reproducible from the facts in this README).
3. BGZF-compressed with `hvantk.core.utils.bgzf.BgzfWriter`, matching the real
   downloader's output compression.
