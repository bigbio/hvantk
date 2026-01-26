# ClinVar Data Streaming Examples

This directory contains examples for working with ClinVar variant annotations using hvantk's data streaming capabilities.

## Contents

**`clinvar_streamer_example.py`** - ClinVar data streaming and chunk processing

Demonstrates:
- Loading ClinVar annotation tables
- Streaming and processing variant data in chunks
- Filtering by gene set or disease terms
- Aggregating and analyzing results

## Quick Start

```bash
# Activate environment
poetry shell

# Run the example
python examples/clinvar/clinvar_streamer_example.py
```

## Example Workflow

The example demonstrates common ClinVar streaming operations:

```python
from hvantk.utils.clinvar_streamer import ClinvarDataStreamer

# Path to ClinVar VCF or Hail Table
clinvar_path = "path/to/clinvar.vcf.bgz"  # or .ht

# Optional: define a gene set to filter
gene_set = {"BRCA1", "BRCA2"}

# Create the streamer
streamer = ClinvarDataStreamer(
    clinvar_path=clinvar_path,
    gene_set=gene_set,  # Optional
    chunk_size=5000,    # Optional, default 10000
)

# Setup the streamer (loads data)
streamer.setup()

try:
    for chunk in streamer.stream():
        # chunk is a Hail Table with a subset of variants
        print(f"Chunk rows: {chunk.count()}")
        # You can process, filter, or export each chunk here
finally:
    streamer.teardown()
```

## Common Operations

### Filter by Gene Set

```python
# Only variants in specified genes
streamer = ClinvarDataStreamer(
    clinvar_path=clinvar_path,
    gene_set={"BRCA1", "BRCA2", "TP53"},
)
streamer.setup()
for chunk in streamer.stream():
    # Process chunk
    pass
streamer.teardown()
```

### Filter by Disease Terms

```python
# Only variants with matching disease terms (case-insensitive, underscores allowed)
streamer = ClinvarDataStreamer(
    clinvar_path=clinvar_path,
    disease_terms={"breast_cancer", "ovarian_cancer"},
)
streamer.setup()
for chunk in streamer.stream():
    # Process chunk
    pass
streamer.teardown()
```

### Aggregating Results

```python
# Aggregate variant counts by gene across all chunks
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

## Building ClinVar Tables

Before using the streamer, you need a ClinVar Hail Table or VCF:

```bash
# Option 1: CLI
hvantk mktable clinvar \
  --raw-input clinvar.vcf.bgz \
  --output-ht clinvar.ht \
  --reference-genome GRCh38

# Option 2: Python API
from hvantk.tables.table_builders import create_clinvar_tb

create_clinvar_tb(
    input_path="clinvar.vcf.bgz",
    output_path="clinvar.ht",
    reference_genome="GRCh38"
)
```

## Advanced Usage

### Custom Chunk Processing

```python
# Example: count pathogenic variants in each chunk
def process_chunk(chunk):
    return chunk.filter(chunk.info.CLNSIG.contains("Pathogenic")).count()

streamer.setup()
for chunk in streamer.stream():
    n_pathogenic = process_chunk(chunk)
    print(f"Pathogenic variants in chunk: {n_pathogenic}")
streamer.teardown()
```

### Using with create_clinvar_training_set_streamer

```python
from hvantk.utils.clinvar_streamer import create_clinvar_training_set_streamer

processor = create_clinvar_training_set_streamer(
    clinvar_path=clinvar_path,
    output_dir="./data/training_set",
    gene_set={"BRCA1", "BRCA2"},
)
result = processor.process()
if result:
    print(f"Generated training set with {result.count()} variants")
```

## Documentation

- [ClinVar Table Builder](../../docs/library/usage.md)
- [Table Builders Guide](../../docs/library/usage.md)
- [Architecture](../../docs/ARCHITECTURE.md)

## Resources

- [ClinVar Database](https://www.ncbi.nlm.nih.gov/clinvar/)
- [ClinVar VCF Format](https://www.ncbi.nlm.nih.gov/clinvar/docs/vcf/)
- [ACMG Variant Interpretation Guidelines](https://www.acmg.net/)

## Tips

1. **Gene set and disease term filtering**: Use the gene_set and disease_terms arguments for targeted streaming
2. **Chunk size**: Adjust chunk_size for memory/performance tradeoffs
3. **Data freshness**: Update ClinVar tables regularly (monthly recommended)
4. **Memory**: Large filters may require sufficient Spark memory

## Troubleshooting

### Table not found

Ensure you've built the ClinVar table first:
```bash
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht
```

### Gene names not matching

ClinVar uses HGNC gene symbols. Check gene naming:
```python
# View available genes
gene_list = streamer.clinvar_ht.aggregate(hl.agg.collect_as_set(streamer.clinvar_ht.gene))
print(gene_list)
```

### Out of memory

For large filters, increase Spark memory or filter in stages:
```python
# Filter incrementally
streamer = ClinvarDataStreamer(clinvar_path=clinvar_path, gene_set=gene_list)
streamer.setup()
for chunk in streamer.stream():
    # Process chunk
    pass
streamer.teardown()
```
