# ClinVar Streaming Example

Example for querying ClinVar variant annotations with hvantk's `ClinVarVariantTableStreamer`.

For full documentation, see the [ClinVar examples guide](https://bigbio.github.io/hvantk/examples/clinvar/).

## Contents

- **`clinvar_streamer_example.py`** - loads a built ClinVar table and filters, aggregates,
  and labels it

## Prerequisites

Build a ClinVar Hail Table first. The ClinVar plugin has a built-in downloader, so this
single command downloads the latest VCF into `data/clinvar/` and builds the table in one go:

```bash
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar \
  --output clinvar.ht
```

If you already have a VCF under `data/clinvar/`, add `--skip-download` to reuse it.

## Quick Start

```bash
python examples/clinvar/clinvar_streamer_example.py
```

Reads `clinvar.ht` in the current directory by default; pass a different path as the first
argument.

## What the Script Demonstrates

- Loading a built ClinVar `AnnotationTable` with `ClinVarVariantTableStreamer.from_path`
- Filtering by gene set (`filter_to_genes`) and by clinical significance
  (`filter_by_pathogenicity`)
- Aggregating clinical-significance counts over a filtered table
- Deriving a TP/TN training-label column with `label_training_set`
