# ensembl-gene

Ensembl BioMart gene annotation provides authoritative gene coordinates,
transcript–protein mappings, gene types, and synonyms for human genes.
The table is exported from Ensembl BioMart with canonical transcripts
filtered and grouped by Ensembl gene ID. It is a primary source for
locus-to-gene mapping in the hvantk annotation pipeline.

Upstream: https://www.ensembl.org/biomart/

## Dataset

- `ensembl-gene:genes` — gene-level annotation table, keyed by gene_id (Ensembl gene ID)

## Build

```bash
hvantk reprocess ensembl-gene:genes \
  --raw-dir <dir-containing-biomart-tsv> \
  --output <out.ht> \
  [--plugin-arg canonical=true] \
  [--plugin-arg fields=gene_id,gene_name]
```

The build is resolved by the plugin loader
(`hvantk/core/plugin/loader.py`) from `plugin.yaml` via
`get_registry().get_dataset("ensembl-gene:genes")` and executed through
`run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`).

The builder is `build_ensembl_gene_genes` in
`hvantk/skills/ensembl_gene/builder.py`, with signature
`build_ensembl_gene_genes(parsed_input, ctx, **params) -> AnnotationTable`.
It imports the BioMart TSV inline via `hl.import_table`, renames fields
using `ENSEMBL_BIOMART_FIELDS` from
`hvantk/skills/ensembl_gene/shared/constants.py`, optionally filters to
canonical transcripts (`canonical` param, default `True`), groups by
`gene_id`, optionally selects `fields`, and wraps the result with
`AnnotationTable.from_hail`.

## Notes

The drift probe is a stub (`drift_probe.py`); a real probe should be
implemented in a follow-up. No downloader is implemented; upstream
BioMart exports are expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_id`. Aggregated
fields per gene include: `transcript_id` (set), `protein_id` (set),
`gene_synonym` (set), `gene_name`, `chromosome`, `gene_start`, `gene_end`,
`gene_type`.

## Tests

A conformance test in `tests/test_ensembl_gene.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/ensembl/ensembl_gene_biomart.tsv.bgz`.
