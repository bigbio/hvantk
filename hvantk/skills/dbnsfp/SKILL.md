# dbnsfp

dbNSFP (database of Non-synonymous functional predictions) is a comprehensive
database of functional predictions and annotations for non-synonymous single
nucleotide variants (SNVs) in the human genome. It aggregates scores from
many tools (SIFT, PolyPhen-2, CADD, REVEL, etc.) and population frequencies
from gnomAD, ExAC, and 1000 Genomes. Version 4.x covers all possible
non-synonymous SNVs on the reference genome.

Upstream: https://sites.google.com/site/jpopgen/dbNSFP

## Dataset

- `dbnsfp:variants` — per-variant functional annotation table, keyed by (locus, alleles)

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The drift probe is a stub; a real probe should be implemented in a
follow-up. The downloader is not implemented; upstream files are
expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `(locus, alleles)`.
Transcript-specific score fields (ending in `_score` or `CADD_phred`) are
parsed into dicts keyed by Ensembl transcript ID. Population frequency
fields from gnomAD, ExAC, 1000Gp3, and ESP6500 are grouped into structs.

## Tests

A conformance test in `tests/test_dbnsfp.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz`.
