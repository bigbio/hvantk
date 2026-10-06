# Third-Party Data Attribution

This document lists third-party datasets committed to this repository as test fixtures. These are real excerpts from public data sources, committed to support reproducible testing. Synthetic fixtures (fabricated values for testing purposes) are not listed here. Test data is excluded from built packages via the `exclude` list in `pyproject.toml`, so this document concerns the git repository only. A pull request that commits a real excerpt as a fixture adds its entry here.

## Expression Atlas

### Fixture: `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/`

**Source**: Expression Atlas experiment E-MTAB-6798, "Mouse RNA-seq time-series of the development of seven major organs"

**Version**: no upstream release version; derived from upstream files retrieved before 2025-05-12 (a full copy was committed in bf2b310c and later removed as unused).

**Files**:
- `E-MTAB-6798-transcripts-tpms.tsv`: Transcript-level TPM expression matrix
- `E-MTAB-6798.condensed-sdrf.tsv`: Sample and data relationship format metadata

**Licence**: CC BY 4.0

**Attribution**: Expression Atlas, EMBL-EBI — https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6798. Citation: "Expression Atlas in 2026: enabling FAIR and open expression data through community collaboration and integration" (*Nucleic Acids Research*, 10 December 2025).

**Modifications**:
- Expression matrix: first 21 transcript rows (20 distinct genes; one gene's second transcript retained to exercise the fixture's multi-transcript-per-gene shape)
- Samples: first 4 of 317 samples (ERR2588382, ERR2588384, ERR2588383, ERR2588399 in header order)
- SDRF: lines for the 4 retained samples only, preserving the original ragged field layout

## AlphaGenome

### Fixture: `hvantk/skills/alphagenome/tests/testdata/raw/alphagenome/`

**Source**: AlphaGenome (Google DeepMind) variant-effect predictions for 411 ClinVar canonical-splice variants, from an AlphaGenome SDK `score_variant` → `variant_scorers.tidy_scores()` run on 2026-06-24 (RECOMMENDED scorers, 1 Mb interval, organism human, API backend).

**Version**: the AlphaGenome SDK and hosted model as of that run; the exact SDK version was not recorded (the run's manifest records the backend, scorer set, interval and organism, but no SDK version). Source shard `scores.shard-00000-of-00001.parquet`, 16,642,316 rows, sha256 `261cf87273f1b36d20db640305f414b2e9bdc55824cddff9e4cc4bcc4c6745a1`.

**Files**:
- `clinvar-subset.parquet`: AlphaGenome scores in the `tidy_scores()` long format
- `NOTICE.md`: the AlphaGenome Output Terms notice that must stay with the data
- `make_fixture.py`: the script that selected the rows (code, not data)

**Other copies of these outputs** (same source and terms):
- `hvantk/skills/alphagenome/tests/snapshots/sample_rows.json`: the builder's summary of the fixture, with scores and track names verbatim; its notice is `tests/snapshots/NOTICE.md` beside it
- `hvantk/skills/alphagenome/tests/test_builder.py`: the oracle constants `AG_SPLICE_SITES` and `HEART_SPLICE_SITE_USAGE`, marked by a comment

**Licence**: AlphaGenome Output Terms of Use (https://deepmind.google.com/science/alphagenome/output-terms): non-commercial use only, and the outputs must not be used to train machine-learning models. Not covered by this repository's MIT licence; the conspicuous notice the terms require is `NOTICE.md` in the fixture directory, with a second one beside the snapshot.

**Attribution**: Google DeepMind, AlphaGenome. ClinVar variant identities are public domain.

**Modifications**:
- Rows only: 3 of the 411 variants (`chr3:39408741:T>C`, `chr6:112216367:C>A`, `chrX:153694448:T>G`) and 150 of the 16,642,316 rows, selected by `make_fixture.py` with seed 20261006 and written in reverse source order. For each variant and scorer it keeps the row with the largest |quantile_score|, the row with the largest |raw_score| and one random row, plus every heart SPLICE_SITE_USAGE row
- Columns: unchanged

## dbNSFP

### Fixture: `hvantk/tests/testdata/raw/dbnsfp/`

**Source**: dbNSFP, academic branch, release v4.9a: functional predictions and annotations for all potential human non-synonymous and splice-site SNVs, distributed through https://www.dbnsfp.org after academic registration.

**Version**: v4.9a (academic). Committed to this repository in efd6c3bd (2025-08-21).

**Files**:
- `dbNSFP4_v49a_example_variants.bgz`: rows of the dbNSFP variant table in its distributed tab-separated layout, BGZF-compressed
- `NOTICE.md`: the licence notice that stays with the data

**Other copies of these data** (same source and terms):
- `hvantk/skills/dbnsfp/tests/snapshots/sample_rows.json`: six fixture rows as the dbNSFP builder outputs them, with the values parsed into typed fields; its notice is `tests/snapshots/NOTICE.md` beside it. `schema.json` there holds only field names and types.

**Licence**: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/), the licence of the dbNSFP academic branch (https://www.dbnsfp.org/license/): non-commercial use only, with attribution, and no sharing of modified versions. Commercial use requires dbNSFP's paid commercial licence, and the CADD, VEST, M-CAP, MutScore, PolyPhen-2, PrimateAI and RGC Million Exome scores in the academic branch also need commercial licences from their authors. Not covered by this repository's MIT licence; the notice is `NOTICE.md` in the fixture directory, with a second one beside the snapshot.

**Attribution**: dbNSFP, © 2024–2026 Genos Bioinformatics LLC — https://www.dbnsfp.org. Citation: Liu X, Li C, Mou C, Dong Y, Tu Y. "dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs." *Genome Medicine* 12, 103 (2020), https://doi.org/10.1186/s13073-020-00803-9.

**Modifications**:
- None to the content: the header (all 458 columns) and the first 4,999 variants on chromosome 10 of the v4.9a variant table, byte-identical to the release once decompressed, stored BGZF-compressed
- To regenerate it, extract whole lines unedited: trimming columns or editing values would make a modified version, which the licence does not allow sharing
