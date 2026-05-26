# ensembl-gene

Ensembl BioMart gene annotation provides authoritative gene coordinates,
transcript–protein mappings, gene types, and synonyms for human genes.
The table is exported from Ensembl BioMart with canonical transcripts
filtered and grouped by Ensembl gene ID. It is a primary source for
locus-to-gene mapping in the hvantk annotation pipeline.

Upstream: https://www.ensembl.org/biomart/

## Dataset

- `ensembl-gene:genes` — gene-level annotation table, keyed by gene_id (Ensembl gene ID)

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The drift probe is a stub; a real probe should be implemented in a
follow-up. The downloader is not implemented; upstream files are
expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_id`. Aggregated
fields per gene include: `transcript_id` (set), `protein_id` (set),
`gene_synonym` (set), `gene_name`, `chromosome`, `gene_start`, `gene_end`,
`gene_type`.

## Tests

A conformance test in `tests/test_ensembl_gene.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/ensembl/ensembl_gene_biomart.tsv.bgz`.
