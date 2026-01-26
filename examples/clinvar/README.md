# ClinVar Data Streaming Examples

This directory contains examples for working with ClinVar variant annotations using hvantk's data streaming capabilities.

## Contents

**`clinvar_streamer_example.py`** - ClinVar data filtering and export

Demonstrates:
- Loading ClinVar annotation tables
- Filtering variants by clinical significance
- Extracting pathogenic variants
- Exporting results to various formats

## Quick Start

```bash
# Activate environment
poetry shell

# Run the example
python examples/clinvar/clinvar_streamer_example.py
```

## Example Workflow

The example demonstrates common ClinVar operations:

```python
from hvantk.data.data_streamer import ClinVarStreamer

# Initialize streamer
streamer = ClinVarStreamer("path/to/clinvar.ht")

# Filter for pathogenic variants
pathogenic = streamer.filter_by_significance(["Pathogenic", "Likely_pathogenic"])

# Filter by genes
brca_variants = streamer.filter_by_genes(["BRCA1", "BRCA2"])

# Filter by review status
high_confidence = streamer.filter_by_stars(min_stars=2)

# Export results
streamer.export_tsv("pathogenic_variants.tsv")
streamer.export_vcf("pathogenic_variants.vcf.bgz")
```

## Common Operations

### Filter by Clinical Significance

```python
# Pathogenic variants only
pathogenic = streamer.filter_by_significance(["Pathogenic"])

# Pathogenic or Likely pathogenic
likely_pathogenic = streamer.filter_by_significance([
    "Pathogenic",
    "Likely_pathogenic"
])

# Exclude VUS
no_vus = streamer.exclude_vus()

# Benign variants
benign = streamer.filter_by_significance(["Benign", "Likely_benign"])
```

### Filter by Genes

```python
# Single gene
brca1 = streamer.filter_by_genes(["BRCA1"])

# Multiple genes
cancer_genes = streamer.filter_by_genes(["BRCA1", "BRCA2", "TP53", "PTEN"])

# From file
streamer.filter_by_genes_file("gene_list.txt")
```

### Filter by Review Status

```python
# At least 2 stars (multiple submitters)
high_quality = streamer.filter_by_stars(min_stars=2)

# At least 3 stars (expert panel)
expert_reviewed = streamer.filter_by_stars(min_stars=3)
```

### Combine Filters

```python
# Pathogenic variants in BRCA1/2 with high confidence
result = (streamer
    .filter_by_genes(["BRCA1", "BRCA2"])
    .filter_by_significance(["Pathogenic"])
    .filter_by_stars(min_stars=2))

result.export_tsv("brca_pathogenic_high_conf.tsv")
```

## ClinVar Significance Levels

| CLNSIG | Description | Common Use |
|--------|-------------|------------|
| Pathogenic | Disease-causing | Disease studies |
| Likely_pathogenic | Probably disease-causing | Disease studies |
| Uncertain_significance | Unknown impact (VUS) | Usually excluded |
| Likely_benign | Probably not harmful | Control sets |
| Benign | Not harmful | Control sets |

## Review Status (Stars)

| Stars | Status | Description |
|-------|--------|-------------|
| 0 | No assertion | Low confidence |
| 1 | Single submitter | Basic confidence |
| 2 | Multiple submitters | Good confidence |
| 3 | Expert panel | High confidence |
| 4 | Practice guideline | Highest confidence |

## Building ClinVar Tables

Before using the streamer, you need a ClinVar Hail Table:

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

## Export Formats

### TSV Export

```python
streamer.export_tsv("output.tsv")
```

Includes columns:
- `locus` - Genomic position (chr:pos)
- `alleles` - Reference and alternate alleles
- `CLNSIG` - Clinical significance
- `CLNREVSTAT` - Review status
- `GENEINFO` - Associated gene(s)
- Additional INFO fields

### VCF Export

```python
streamer.export_vcf("output.vcf.bgz")
```

Standard VCF format with ClinVar annotations in INFO field.

### Hail Table

```python
filtered_ht = streamer.get_table()
filtered_ht.write("filtered_clinvar.ht")
```

Preserves full Hail Table structure for downstream analysis.

## Advanced Usage

### Custom Filters

```python
# Filter by allele frequency (if annotated)
rare_variants = streamer.filter_expr("ht.info.AF < 0.01")

# Filter by variant type
snvs = streamer.filter_expr("hl.len(ht.alleles[0]) == 1 & hl.len(ht.alleles[1]) == 1")
```

### Statistics

```python
# Count variants
total = streamer.count()

# Get significance distribution
sig_counts = streamer.get_significance_counts()

# Get per-gene counts
gene_counts = streamer.get_gene_counts()
```

## Use Cases

### Case 1: Disease Variant Discovery

```python
# Find pathogenic variants in candidate genes
disease_variants = (ClinVarStreamer("clinvar.ht")
    .filter_by_genes_file("candidate_genes.txt")
    .filter_by_significance(["Pathogenic", "Likely_pathogenic"])
    .filter_by_stars(min_stars=1)
    .export_tsv("disease_variants.tsv"))
```

### Case 2: Benign Control Set

```python
# Extract high-confidence benign variants
benign_controls = (ClinVarStreamer("clinvar.ht")
    .filter_by_significance(["Benign", "Likely_benign"])
    .filter_by_stars(min_stars=2)
    .export_vcf("benign_controls.vcf.bgz"))
```

### Case 3: Expert-Curated Variants

```python
# Get only expert-reviewed pathogenic variants
expert_path = (ClinVarStreamer("clinvar.ht")
    .filter_by_stars(min_stars=3)
    .filter_by_significance(["Pathogenic"])
    .get_table())
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

1. **Review status matters**: Use `min_stars=2` for reliable results
2. **VUS variants**: Usually excluded from analysis (uncertain significance)
3. **Multiple genes**: ClinVar variants can affect multiple genes
4. **Data freshness**: Update ClinVar tables regularly (monthly recommended)
5. **Memory**: Large filters may require sufficient Spark memory

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
gene_list = streamer.get_unique_genes()
print(gene_list)
```

### Out of memory

For large filters, increase Spark memory or filter in stages:
```python
# Filter incrementally
filtered = (streamer
    .filter_by_significance(["Pathogenic"])
    .checkpoint("temp.ht")  # Checkpoint intermediate result
    .filter_by_genes(gene_list))
```
