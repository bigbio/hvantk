# pqtl

Protein quantitative trait loci (pQTL) link genetic variants to protein
abundance levels. The primary source currently supported is Fang et al. (2025),
which provides pQTL allpairs from 5 tissues profiled by TMT mass spectrometry
in the GTEx cohort. Gene symbols are mapped to Ensembl gene IDs via HGNC for
compatibility with the eQTL cascade join. The table is keyed by
(locus, alleles, gene_id).

## Dataset

- `pqtl:metrics` — per-variant-gene pQTL association metrics, keyed by (locus, alleles, gene_id)

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The builder uses a delegation-stub to the legacy `create_pqtl_tb` because
the pQTL builder requires an external HGNC table resource and complex
gene-symbol → Ensembl ID mapping that shares helpers with other builders.
The drift probe is a stub. The downloader is not implemented.

## Schema

Field documentation TBD. Key fields: `locus`, `alleles`, `gene_id`,
`beta`, `se`, `p_value`, `tissue`, `source`, `is_cis`.

## Tests

No fixture is available for the pqtl source. The test file contains a
registration-only test and a skipped round-trip test.
