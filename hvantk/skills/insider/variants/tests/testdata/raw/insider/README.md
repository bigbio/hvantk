# Synthetic INSIDER (variants) fixture

SYNTHETIC. Contains no rows from the real INSIDER BED
(`Whole_Human_Interactome_Interface_hg38.bed`, ~1.17 GB, 208,448 PPI tracks, 18.6M data
rows; see `hvantk/skills/insider/catalog/datasets.json`, accession `INSIDER_v1.0`).

## Why synthetic

The previous fixture (`insider_sample.bed`, removed here) was a real excerpt, sliced
from the live BED by a track-aware sub-sampler. #426 asked for a licence audit of every
real fixture without a `THIRD_PARTY_DATA.md` entry. For INSIDER (Yu lab / Interactome
Insider, Wei et al. 2017), no published licence or terms-of-use statement could be
found, checked on 2026-10-13:

- `http://interactomeinsider.yulab.org/` carries no About/Downloads/Terms page with a
  licence statement (the nav's "About" and "Downloads" items are in-page JS tabs, not
  separate URLs; no `/download.html` or similar page exists — confirmed 404).
- The underlying paper (Meyer MJ et al., "Interactome INSIDER: a structural
  interactome browser for genomic studies", *Nat Methods* 15, 290-298 (2018), PMID
  29355848, DOI 10.1038/nmeth.4540) carries no open licence tag in Europe PMC
  (`license: None`); the project's own `SKILL.md` citation of PMID 29036289 for this
  paper is itself wrong (that PMID belongs to an unrelated Bioinformatics paper,
  "VICTOR: genome-based phylogeny and classification of prokaryotic viruses" — a
  separate, pre-existing documentation defect, not fixed here).
- The closest formal licence signal found is the bioRxiv preprint of the same work
  (DOI 10.1101/126862), tagged **CC BY-NC-ND** in Europe PMC — non-commercial, no
  derivatives, and in any case a licence on the preprint text, not a stated licence for
  the downloadable BED/TXT data files themselves.

An absent/unclear licence limits what a fixture may contain, not whether one exists
(`_conventions` SKILL.md §9), so this fixture is now a synthetic, format-faithful BED
rather than a real excerpt.

## Format facts preserved from the real file

- `browser hide all` directive line, followed by `track name=<P1>_ppi_<P2>
  description="Interface Residues for <P1>_<P2> (SOURCE: Predicted Interface)"
  visibility=dense itemRgb="On"` headers, one per PPI pair, each followed by its data
  rows (a plain `head -N` would split a track block, which is why the original
  sub-sampler was track-aware).
- Data rows: `<contig> <start> <end> <name> <score> <strand> <thickStart> <thickEnd>
  <itemRgb>` (tab-separated, UCSC-style BED, `itemRgb` fixed at `247,176,91` as in the
  live file).
- The same interval can appear under more than one `track`, which the builder
  aggregates into one row keyed by `interval` with `ppi_ids` as an array — this
  fixture keeps that shape: `chr3:2000301-2000303` (the Hail-interval form of BED row
  `chr3 2000300 2000302`) appears under all four `SYNTH0003_ppi_*` tracks.
- The first track (`SYNTH0001_ppi_SYNTH0002`) opens with a zero-length row
  (`1000000 1000000`), preserving the real file's "first track has a filtered,
  zero-length opening row" shape noted in `test_builder.py`.

## What is fabricated

- All "UniProt accession" strings (`SYNTH0001`-`SYNTH0007`) are placeholders, not real
  UniProt IDs.
- All genomic coordinates are arbitrary round numbers on `chr11`/`chr3`, chosen only to
  stay within real GRCh38 chromosome bounds; they are not real interface positions and
  were not checked against any real protein's sequence.

## Generation

1. Verified the absence of a published licence as described above (2026-10-13).
2. Hand-authored 5 tracks / 20 data rows reproducing the real file's block structure
   (one zero-length-opening track, one multi-row track, and three tracks sharing one
   interval to exercise the multi-PPI-per-interval aggregation), replacing all
   accessions and coordinates with placeholders.
3. Regenerated `../../snapshots/schema.json` and `../../snapshots/sample_rows.json` via
   `pytest hvantk/skills/insider/variants/tests/test_builder.py -m hail
   --regenerate-snapshots`.
