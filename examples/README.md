# hvantk Examples

This directory contains practical examples demonstrating how to use hvantk for various genomic data analysis tasks.

## Quick Links

- [Recipe Templates](#recipe-templates)
- [Python Examples](#python-examples)
- [Output Examples](#output-examples)

## Recipe Templates

The `recipes/` directory contains ready-to-use recipe templates for batch processing:

### Table Recipes

**[tables.example.json](recipes/tables.example.json)** - Batch-create annotation tables (Hail Tables)

```bash
# Edit the recipe with your paths
# Then run:
hvantk mktable-batch --recipe recipes/tables.example.json
```

Demonstrates:
- ClinVar annotation table creation
- INSIDER interactome table creation
- Reference genome specification
- Multiple output formats

### MatrixTable Recipes

**[matrices.example.json](recipes/matrices.example.json)** - Batch-create expression matrices

```bash
hvantk mkmatrix-batch --recipe recipes/matrices.example.json
```

Demonstrates:
- UCSC Cell Browser data conversion
- Expression Atlas data processing
- Multiple matrix creation in one run

**[cptac.example.json](recipes/cptac.example.json)** - CPTAC protein expression data

```bash
hvantk mkmatrix-batch --recipe recipes/cptac.example.json
```

Demonstrates:
- Protein expression matrix creation
- Categorical metadata handling
- Custom column mapping

### Analysis Recipes

**[heart_rare_variants.yaml](recipes/heart_rare_variants.yaml)** - Example analysis workflow

Demonstrates:
- YAML recipe format
- Multi-step analysis pipeline
- Combining multiple data sources

## Python Examples

### ClinVar Streamer Examples

**[clinvar_streamer_examples.py](clinvar_streamer_examples.py)** - Working with ClinVar data

```bash
python examples/clinvar_streamer_examples.py
```

Demonstrates:
- Loading ClinVar annotations
- Filtering variants by clinical significance
- Extracting pathogenic variants
- Exporting results

### HGC Quality Control

**[hgc_qc_example.py](hgc_qc_example.py)** - Joint genotyping QC workflow

```bash
python examples/hgc_qc_example.py
```

Demonstrates:
- Computing QC metrics for combined cohorts
- Generating QC visualizations
- Creating HTML reports
- Interactive plots
- Quality-based filtering

Features:
- Sample-level QC (call rate, heterozygosity)
- Variant-level QC (call rate, Hardy-Weinberg)
- Ti/Tv ratio analysis
- Professional HTML reports

## Output Examples

### QC Output

**[qc_output_20251127_204630/](qc_output_20251127_204630/)** - Example QC report output

Contains:
- `qc_report_20251127_204630.html` - Interactive HTML report
- `qc_dashboard_20251127_204630.png` - QC dashboard visualization

Open the HTML file in a browser to see:
- Summary statistics
- Interactive plots
- Quality recommendations
- Filtering suggestions

## Usage Patterns

### Running Recipes

```bash
# JSON recipe (built-in support)
hvantk mktable-batch --recipe path/to/recipe.json

# YAML recipe (requires PyYAML)
hvantk mkmatrix-batch --recipe path/to/recipe.yaml
```

### Running Python Examples

```bash
# Activate environment first
poetry shell

# Run example
python examples/clinvar_streamer_examples.py

# Or with custom parameters
python examples/hgc_qc_example.py --input my_data.mt
```

## Customizing Examples

### Modifying Recipes

1. Copy an example recipe:
   ```bash
   cp recipes/tables.example.json my_recipe.json
   ```

2. Edit paths and parameters:
   ```json
   {
     "tables": [
       {
         "name": "clinvar",
         "input": "/path/to/your/clinvar.vcf.bgz",
         "output": "/path/to/your/output.ht",
         "params": {"reference_genome": "GRCh38"}
       }
     ]
   }
   ```

3. Run your recipe:
   ```bash
   hvantk mktable-batch --recipe my_recipe.json
   ```

### Adapting Python Examples

Most examples accept command-line arguments:

```python
# Check available options
python examples/hgc_qc_example.py --help

# Run with custom input
python examples/hgc_qc_example.py \
  --input /path/to/cohort.mt \
  --output /path/to/qc_report.html
```

## Example Datasets

For testing, use the small test datasets in `hvantk/tests/testdata/`:

```bash
# Example with test data
hvantk mkmatrix ucsc \
  -e hvantk/tests/testdata/raw/ucsc/exprMatrix.test.tsv.bgz \
  -m hvantk/tests/testdata/raw/ucsc/meta.test.tsv \
  -o /tmp/test_output.mt
```

## Common Workflows

### Workflow 1: Build Annotation Database

```bash
# 1. Create recipe with all your annotation sources
cat > annotations.json << 'EOF'
{
  "tables": [
    {"name": "clinvar", "input": "clinvar.vcf.bgz", "output": "clinvar.ht"},
    {"name": "gnomad-metrics", "input": "gnomad.tsv.bgz", "output": "gnomad.ht"}
  ]
}
EOF

# 2. Build all tables at once
hvantk mktable-batch --recipe annotations.json
```

### Workflow 2: Expression Analysis

```bash
# 1. Build expression matrix
hvantk mkmatrix ucsc \
  --expression-matrix expr.tsv.bgz \
  --metadata meta.tsv \
  --output-mt heart_expr.mt

# 2. Use in downstream analysis
python your_analysis.py --expression heart_expr.mt
```

### Workflow 3: Joint Genotyping with QC

```bash
# 1. Combine GVCFs
hvantk hgc gvcf-combine -g gvcfs/ -o cohort.vds

# 2. Convert to MatrixTable
hvantk hgc vds2mt -i cohort.vds -o cohort.mt --adjust-genotypes

# 3. Generate QC report
hvantk hgc qc-report -i cohort.mt -o qc_report.html

# 4. Filter based on QC
hvantk hgc filter-qc -i cohort.mt -o filtered.mt --min-sample-call-rate 0.95
```

## Tips and Best Practices

1. **Start small**: Test with small datasets first
2. **Use recipes**: Automate repetitive tasks with recipe files
3. **Check outputs**: Verify outputs before running large-scale analyses
4. **Save recipes**: Keep recipe files for reproducibility
5. **Document custom recipes**: Add comments to complex recipes

## Need Help?

- **Documentation**: See [docs/](../docs/)
- **Usage Guide**: [docs/library/usage.md](../docs/library/usage.md)
- **HGC Documentation**: [docs/tools/hgc.md](../docs/tools/hgc.md)
- **Architecture**: [docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md)
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)

## Contributing Examples

Have a useful example? Please contribute!

1. Add your example to this directory
2. Document it in this README
3. Include sample data or reference test data
4. Submit a pull request

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.
