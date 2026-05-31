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

## Build

The builder is `build_dbnsfp_variants` in `hvantk/skills/dbnsfp/builder.py`
(signature `(parsed_input, ctx, **params) -> AnnotationTable`). The plugin
loader resolves it from `plugin.yaml` via `get_registry().get_dataset("dbnsfp:variants")`;
top-level builds run through `run_builder_for_spec`.

Invoke via the CLI:

```bash
hvantk reprocess dbnsfp:variants --raw-dir <dir> --output <out.ht> \
  [--plugin-arg reference_genome=GRCh38] [--plugin-arg parse_transcript_scores=true]
```

Supported `params`: `reference_genome` (default `GRCh38`),
`min_partitions` (default 200), `force_bgz` (default true),
`parse_transcript_scores` (default true), `group_prefixes` (list of str),
`auto_convert_bgz` (default false).

## Notes

The drift probe (`drift_probe.fetch_fingerprint`) is a stub; a real probe
should be implemented in a follow-up. No downloader is wired in
`plugin.yaml` lifecycle yet; upstream files are expected to be externally
materialized for now.

## Schema

Field documentation TBD. The table is keyed by `(locus, alleles)`.
Transcript-specific score fields (ending in `_score` or `CADD_phred`) are
parsed into dicts keyed by Ensembl transcript ID. Population frequency
fields from gnomAD, ExAC, 1000Gp3, and ESP6500 are grouped into structs.

## Tests

A conformance test in `hvantk/skills/dbnsfp/tests/test_dbnsfp.py` exercises
the build via `run_builder_for_spec`. Per `plugin.yaml`, test artifacts are
plugin-relative: the raw fixture lives under `tests/testdata/raw/dbnsfp/`
(`dbNSFP4_v49a_example_variants.bgz`), with the drift fingerprint at
`tests/drift_fingerprint.json` and snapshots at `tests/snapshots/`. Run with
`pytest hvantk/skills/dbnsfp/tests`.
