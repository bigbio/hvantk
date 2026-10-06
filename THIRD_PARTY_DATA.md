# Third-Party Data Attribution

This document lists third-party datasets committed to this repository as test fixtures. These are real excerpts from public data sources, committed to support reproducible testing. Synthetic fixtures (fabricated values for testing purposes) are not listed here. Test data is excluded from built packages via the `exclude` list in `pyproject.toml`, so this document concerns the git repository only. A pull request that commits a real excerpt as a fixture adds its entry here.

## Expression Atlas

### Fixture: `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/`

**Source**: Expression Atlas experiment E-MTAB-6798, "Mouse RNA-seq time-series of the development of seven major organs"

**Version**: no upstream release version; retrieved before 2025-05-12 (the files were first committed in bf2b310c).

**Files**:
- `E-MTAB-6798-transcripts-tpms.tsv`: Transcript-level TPM expression matrix
- `E-MTAB-6798.condensed-sdrf.tsv`: Sample and data relationship format metadata

**Licence**: CC BY 4.0

**Attribution**: Expression Atlas, EMBL-EBI — https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6798. Citation: "Expression Atlas in 2026: enabling FAIR and open expression data through community collaboration and integration" (*Nucleic Acids Research*, 10 December 2025).

**Modifications**:
- Expression matrix: first 21 transcript rows (20 distinct genes; one gene's second transcript retained to exercise the fixture's multi-transcript-per-gene shape)
- Samples: first 4 of 317 samples (ERR2588382, ERR2588384, ERR2588383, ERR2588399 in header order)
- SDRF: lines for the 4 retained samples only, preserving the original ragged field layout

### Upstream Files: `hvantk/tests/testdata/raw/expression_atlas/`

**Source**: Expression Atlas experiment E-MTAB-6798

**Version**: no upstream release version; retrieved before 2025-05-12 (the files were first committed in bf2b310c).

**Files**:
- `E-MTAB-6798-transcripts-tpms.tsv.bgz`: Full transcript-level TPM matrix (~116k transcripts × 317 samples), compressed
- `E-MTAB-6798.condensed-sdrf.tsv`: Full sample and data relationship format metadata

**Licence**: CC BY 4.0

**Attribution**: Expression Atlas, EMBL-EBI — https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6798. Citation: "Expression Atlas in 2026: enabling FAIR and open expression data through community collaboration and integration" (*Nucleic Acids Research*, 10 December 2025).

**Modifications**: None. These are the unmodified upstream files.
