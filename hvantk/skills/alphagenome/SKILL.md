# alphagenome

AlphaGenome is a deep learning model from Google DeepMind that predicts the
functional effects of genetic variants on gene expression, chromatin accessibility,
and other molecular phenotypes at nucleotide resolution. The hvantk builder runs
the AlphaGenome API via the AlphaGenomePipeline for a given set of input variants
and produces a Hail Table keyed by (locus, alleles).

## Dataset

- `alphagenome:predictions` — per-variant AlphaGenome effect predictions, keyed by (locus, alleles)

## Builder

The builder is `build_alphagenome_predictions` in
`hvantk/skills/alphagenome/builder.py`, with the standard plugin signature
`(parsed_input, ctx, **params) -> AnnotationTable`.
Internally it drives `AlphaGenomePipeline` (in
`hvantk/skills/alphagenome/pipelines.py`) to call the external API for each
variant, then builds the Hail Table inline from the input variants
(`hl.read_table` for `.ht` input, or `hl.import_table` + `hl.locus` for a TSV
with `chrom`/`pos`/`ref`/`alt` columns) and keys it by `(locus, alleles)`.

A `config_path` param pointing to an AlphaGenome YAML config is required (see
`tests/testdata/alphagenome_config.yaml`); `no_resume` (bool, default False) is
optional. There is no built-in downloader (no `lifecycle.download` in
`plugin.yaml`). The drift probe (`drift_probe.py`, `fetch_fingerprint`) is a
placeholder stub.

## Build invocation

```bash
hvantk reprocess alphagenome:predictions \
  --raw-dir <dir-with-input-variants> \
  --output <out.ht> \
  --plugin-arg config_path=<path/to/alphagenome_config.yaml>
```

## Schema

Field documentation TBD. The table is keyed by `(locus, alleles)`.
Prediction outputs depend on the AlphaGenome config and model version.

## Tests

```bash
pytest hvantk/skills/alphagenome/tests
```

No raw-data fixture is available for the alphagenome source (requires AlphaGenome
API access). `tests/test_alphagenome.py` contains a registration-only test
(`test_alphagenome_predictions_registered`) and a skipped round-trip test. The
only checked-in fixture is `tests/testdata/alphagenome_config.yaml`; the schema/row
snapshot and drift fingerprint paths declared in `plugin.yaml` are not yet
populated.
