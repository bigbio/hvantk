# NOTICE: AlphaGenome output data

**These AlphaGenome outputs are provided under and subject to the AlphaGenome
Output Terms of Use (<https://deepmind.google.com/science/alphagenome/output-terms>).
Non-commercial use only; they must not be used to train machine-learning models.
These outputs are not covered by this repository's MIT licence.**

## What this is

`clinvar-subset.parquet` is a small, real subset of AlphaGenome prediction
scores, in the AlphaGenome SDK's `variant_scorers.tidy_scores()` long format
(one row per variant x output track x metric).

## Modifications

A subset of 3 variants and 150 rows of an AlphaGenome SDK `score_variant` ->
`variant_scorers.tidy_scores()` run on 2026-06-24 (RECOMMENDED scorers, 1 Mb
interval, human, API backend). For each variant and each of the 19 scorers it
keeps the row with the largest |quantile_score|, the row with the largest
|raw_score| and one random row, plus every heart SPLICE_SITE_USAGE row, and
writes them in reverse source order; `make_fixture.py` records exactly how.
Columns unchanged from the original AlphaGenome output; rows were only selected
and reordered.

## ClinVar identities

ClinVar variant identities are public domain.
