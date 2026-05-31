# pqtl

Protein quantitative trait loci (pQTL) link genetic variants to protein
abundance levels. The primary source currently supported is Fang et al. (2025),
which provides pQTL allpairs from 5 tissues profiled by TMT mass spectrometry
in the GTEx cohort. Gene symbols are mapped to Ensembl gene IDs via HGNC for
compatibility with the eQTL cascade join. The table is keyed by
(locus, alleles, gene_id).

## Dataset

- `pqtl:metrics` — per-variant-gene pQTL association metrics, keyed by (locus, alleles, gene_id)

## Build

```bash
hvantk reprocess pqtl:metrics \
  --raw-dir <dir-with-fang-allpairs> \
  --output <out.ht> \
  --plugin-arg hgnc_ht=<hgnc-lookup.ht>
```

The HGNC Hail Table (built via `hvantk reprocess hgnc:lookup`) is required so
gene symbols can be mapped to Ensembl gene IDs for the eQTL cascade join. To
opt out and produce a symbol-keyed table for non-cascade use, pass
`--plugin-arg no_gene_map=true`.

## Builder

The builder `build_pqtl_metrics` lives in `hvantk/skills/pqtl/builder.py`. Its
signature is `(parsed_input, ctx, **params) -> AnnotationTable`. It imports the
Fang allpairs inline via `hl.import_table`, parses GTEx variant IDs with
`parse_gtex_variant_id` from `hvantk/core/utils/qtl_helpers.py`, derives SE as
`|BETA / STAT|` (Fang files lack an SE column), and maps gene symbols to
Ensembl gene IDs through a `GeneCatalogStreamer` (base class in
`hvantk/core/streamers/gene_catalog.py`). The drift probe is a stub. The
downloader is not implemented.

## Schema

Field documentation TBD. Key fields: `locus`, `alleles`, `gene_id`,
`beta`, `se`, `p_value`, `tissue`, `source`, `is_cis`.

## Tests

```bash
pytest hvantk/skills/pqtl/tests
```

No fixture is available for the pqtl source. `hvantk/skills/pqtl/tests/test_pqtl.py`
contains a registration-only test (`test_pqtl_metrics_registered`) and a skipped
round-trip test (`test_pqtl_metrics_round_trip`). The `tests:` block in
`plugin.yaml` declares plugin-relative fixture/snapshot paths
(`tests/testdata/raw/pqtl`, `tests/snapshots/schema.json`,
`tests/snapshots/sample_rows.json`, `tests/drift_fingerprint.json`) that are
not yet populated.
