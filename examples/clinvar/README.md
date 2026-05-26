# ClinVar Streaming Example

Example for working with ClinVar variant annotations using hvantk's data streaming API.

For full documentation, see the [ClinVar examples guide](https://bigbio.github.io/hvantk/examples/clinvar/).

## Contents

- **`clinvar_streamer_example.py`** - ClinVar data streaming and chunk processing

## Prerequisites

Build a ClinVar Hail Table first:

```bash
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar \
  --output clinvar.ht
```

## Quick Start

```bash
python examples/clinvar/clinvar_streamer_example.py
```

## What the Script Demonstrates

- Loading ClinVar annotation tables
- Streaming and processing variant data in chunks
- Filtering by gene set or disease terms
- Aggregating and analyzing results
