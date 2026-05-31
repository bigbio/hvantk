# cosmic-cgc

The COSMIC Cancer Gene Census (CGC) is a curated catalogue of genes with
causal roles in human cancer. Each entry includes a tier classification
(Tier 1 — strong experimental evidence; Tier 2 — literature-supported),
mutation types, tumour types (somatic/germline), and role in cancer
(oncogene/TSG/fusion). COSMIC is maintained by the Wellcome Sanger Institute.

Upstream: https://cancer.sanger.ac.uk/census

## Dataset

- `cosmic-cgc:submissions` — gene-level cancer gene census table, keyed by `gene_symbol` by default (or `hgnc_id` if a `gene_catalog` is provided)

## Build

```bash
hvantk reprocess cosmic-cgc:submissions --raw-dir <dir> --output <out.ht>
```

Optional plugin args (passed via `--plugin-arg KEY=VALUE`):

- `mutation_context` — `both` (default), `somatic`, or `germline`
- `min_classification` — filter to tiers at or above this level
- `fields` — restrict to a subset of output fields

The builder function is `build_cosmic_cgc_submissions` in
`hvantk/skills/cosmic_cgc/builder.py`, with signature
`build_cosmic_cgc_submissions(parsed_input, ctx, **params) -> AnnotationTable`.
It imports the COSMIC CGC TSV inline via `hl.import_table`, renames fields,
normalises tier classifications, parses multi-value fields, and wraps the
result with provenance (`schema_id="cosmic-cgc-v1"`). The downloader is not
implemented; upstream files are expected to be externally materialized for now.

## Schema

Field documentation TBD. The table is keyed by `gene_symbol` by default,
or by `hgnc_id` if a `gene_catalog` (an `HGNCGeneCatalogStreamer`) is passed
as a build parameter to resolve gene symbols to HGNC IDs. Multi-value
fields (`tumour_types_somatic`, `tumour_types_germline`, `role_in_cancer`,
`mutation_types`) are parsed into arrays.

## Tests

```bash
pytest hvantk/skills/cosmic_cgc/tests
```

Test fixtures and snapshots are plugin-relative (declared in `plugin.yaml`):

- fixture: `tests/testdata/raw/cosmic-cgc`
- schema snapshot: `tests/snapshots/schema.json`
- row snapshot: `tests/snapshots/sample_rows.json`
- drift fingerprint: `tests/drift_fingerprint.json`

Note: COSMIC requires account login, so a full raw fixture may not be
available; tests may be limited to registration/round-trip coverage.
