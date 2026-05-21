# alphagenome

AlphaGenome is a deep learning model from Google DeepMind that predicts the
functional effects of genetic variants on gene expression, chromatin accessibility,
and other molecular phenotypes at nucleotide resolution. The hvantk builder runs
the AlphaGenome API via the AlphaGenomeStreamer for a given set of input variants
and produces a Hail Table keyed by (locus, alleles).

## Dataset

- `alphagenome:predictions` — per-variant AlphaGenome effect predictions, keyed by (locus, alleles)

## Phase K notes

This plugin was promoted from the hardcoded `_TABLE_BUILDERS` entry of
the same name as part of Phase K of the data-model platform refactor.
The builder uses a delegation-stub to the legacy `create_alphagenome_tb`
because the AlphaGenome builder wraps the AlphaGenomeStreamer with complex
setup/teardown. The drift probe is a stub. The downloader is not implemented;
a config_path parameter pointing to an AlphaGenome YAML config is required.

## Schema

Field documentation TBD. The table is keyed by `(locus, alleles)`.
Prediction outputs depend on the AlphaGenome config and model version.

## Tests

No fixture is available for the alphagenome source (requires AlphaGenome
API access). The test file contains a registration-only test and a skipped
round-trip test.
