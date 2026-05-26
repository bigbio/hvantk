# ClinVar Streaming Examples

Stream and filter ClinVar variant annotations using hvantk's `ClinvarDataStreamer`.

## Prerequisites

Build a ClinVar Hail Table first. The ClinVar plugin has a built-in downloader, so this single command downloads the latest VCF into `data/clinvar/` and builds the table in one go:

```bash
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar/ \
  --output clinvar.ht \
  --plugin-arg reference_genome=GRCh38
```

If you already have `clinvar.vcf.bgz` (or `.vcf.gz`), place it under `data/clinvar/` and add `--skip-download`.

## Basic Streaming

```python
from hvantk.data.clinvar_streamer import ClinvarDataStreamer

streamer = ClinvarDataStreamer(
    clinvar_path="clinvar.ht",
    chunk_size=5000,
)

streamer.setup()
try:
    for chunk in streamer.stream():
        print(f"Chunk rows: {chunk.count()}")
finally:
    streamer.teardown()
```

## Filtering by Gene Set

```python
streamer = ClinvarDataStreamer(
    clinvar_path="clinvar.ht",
    gene_set={"BRCA1", "BRCA2", "TP53"},
)
streamer.setup()
for chunk in streamer.stream():
    # Process chunk
    pass
streamer.teardown()
```

## Filtering by Disease Terms

```python
streamer = ClinvarDataStreamer(
    clinvar_path="clinvar.ht",
    disease_terms={"breast_cancer", "ovarian_cancer"},
)
streamer.setup()
for chunk in streamer.stream():
    pass
streamer.teardown()
```

## Aggregating Results

```python
from collections import Counter

gene_counts = Counter()
streamer.setup()
for chunk in streamer.stream():
    rows = chunk.select("gene").collect()
    for row in rows:
        gene_counts[row.gene] += 1
streamer.teardown()
print(gene_counts)
```

## Training Set Generation

```python
from hvantk.data.clinvar_streamer import create_clinvar_training_set_streamer

processor = create_clinvar_training_set_streamer(
    clinvar_path="clinvar.ht",
    output_dir="./data/training_set",
    gene_set={"BRCA1", "BRCA2"},
)
result = processor.process()
if result:
    print(f"Generated training set with {result.count()} variants")
```

## Tips

- **Chunk size**: Adjust for memory/performance tradeoffs (default: 10000)
- **Data freshness**: Update ClinVar tables regularly (monthly recommended)
- **Gene names**: ClinVar uses HGNC symbols; verify your gene names match

## Runnable Scripts

See the [`examples/clinvar/`](https://github.com/bigbio/hvantk/tree/main/examples/clinvar/) directory for a complete runnable example.
