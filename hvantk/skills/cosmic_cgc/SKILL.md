# cosmic-cgc

The COSMIC Cancer Gene Census (CGC) is a curated catalogue of genes with
causal roles in human cancer. Each entry includes a tier classification
(Tier 1 — strong experimental evidence; Tier 2 — literature-supported),
mutation types, tumour types (somatic/germline), and role in cancer
(oncogene/TSG/fusion). COSMIC is maintained by the Wellcome Sanger Institute.

Upstream: https://cancer.sanger.ac.uk/census

## Dataset

- `cosmic-cgc:submissions` — gene-level cancer gene census table, keyed by gene_symbol (or hgnc_id if hgnc_path is provided)

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The drift probe is a stub; a real probe should be implemented in a
follow-up. The downloader is not implemented; upstream files are
expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_symbol` by default,
or by `hgnc_id` if the `hgnc_path` parameter is provided. Multi-value
fields (`tumour_types_somatic`, `tumour_types_germline`, `role_in_cancer`,
`mutation_types`) are parsed into arrays.

## Tests

No fixture is available for the cosmic-cgc source (COSMIC requires account
login). The test file contains a registration-only test and a skipped
round-trip test.
