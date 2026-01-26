# hvantk Examples

This directory contains end-to-end workflow examples demonstrating how to use hvantk for genomic variant analysis. Each subdirectory represents a complete analysis pipeline with scripts, results, and documentation.

## Quick Links

- [HGC (Joint Genotyping)](#hgc-joint-genotyping)
- [PSROC (Prediction Score ROC Analysis)](#psroc-prediction-score-roc-analysis)
- [ClinVar Data Streaming](#clinvar-data-streaming)
- [Recipe Templates](#recipe-templates)

## Workflow Examples

### HGC (Joint Genotyping)

**Directory:** [`hgc/`](hgc/)

Complete pipeline for joint genotyping of GVCF cohorts with quality control and benchmarking.

**Key scripts:**
- `hgc_qc_example.py` - QC workflow with visualization
- `hgc_cpu_scaling_benchmark.py` - CPU scalability testing
- `hgc_scalability_benchmark.py` - Sample size scalability

**Quick start:**
```bash
python examples/hgc/hgc_qc_example.py
```

**Outputs:** QC reports (HTML), dashboards (PNG), metrics (JSON)

**Documentation:** [HGC README](hgc/README.md) | [HGC Docs](../docs/tools/hgc.md)

---

### PSROC (Prediction Score ROC Analysis)

**Directory:** [`psroc/`](psroc/)

ROC curve analysis for evaluating variant pathogenicity prediction scores using ClinVar labels.

**Key scripts:**
- `run_psroc_example.py` - Complete PSROC workflow with synthetic data

**Quick start:**
```bash
python examples/psroc/run_psroc_example.py
```

**Outputs:** ROC curves (PNG), AUC metrics (JSON), annotated variants (TSV)

**Documentation:** [PSROC README](psroc/README.md) | [PSROC Docs](../docs/tools/psroc.md)

---

### ClinVar Data Streaming

**Directory:** [`clinvar/`](clinvar/)

Examples for filtering and processing ClinVar variant annotations.

**Key scripts:**
- `clinvar_streamer_example.py` - ClinVar filtering and export

**Quick start:**
```bash
python examples/clinvar/clinvar_streamer_example.py
```

**Outputs:** Filtered variant tables (TSV/VCF)

**Documentation:** [ClinVar README](clinvar/README.md)

---

## Recipe Templates

**Directory:** [`recipes/`](recipes/)

Ready-to-use recipe templates for batch processing annotation tables and expression matrices.

### Table Recipes

**`tables.example.json`** - Batch-create annotation tables (Hail Tables)

```bash
hvantk mktable-batch --recipe recipes/tables.example.json
```

Demonstrates:
- ClinVar annotation table creation
- INSIDER interactome table creation
- Reference genome specification
- Multiple output formats

### MatrixTable Recipes

**`matrices.example.json`** - Batch-create expression matrices

```bash
hvantk mkmatrix-batch --recipe recipes/matrices.example.json
```

Demonstrates:
- UCSC Cell Browser data conversion
- Expression Atlas data processing
- Multiple matrix creation in one run

**`cptac.example.json`** - CPTAC protein expression data

```bash
hvantk mkmatrix-batch --recipe recipes/cptac.example.json
```

Demonstrates:
- Protein expression matrix creation
- Categorical metadata handling
- Custom column mapping

### Analysis Recipes

**`heart_rare_variants.yaml`** - Example analysis workflow (YAML format)

Demonstrates:
- Multi-step analysis pipeline
- Combining multiple data sources
- Downstream analysis patterns

---

## Running Examples

### Prerequisites

```bash
# Install dependencies
poetry install

# Activate environment
poetry shell
```

### Basic Usage

```bash
# Run a workflow example
python examples/<workflow>/script.py

# Run with custom parameters (if supported)
python examples/hgc/hgc_qc_example.py --input my_data.mt --output my_results/
```

### Using Recipes

```bash
# JSON recipe (built-in support)
hvantk mktable-batch --recipe path/to/recipe.json

# YAML recipe (requires PyYAML)
hvantk mkmatrix-batch --recipe path/to/recipe.yaml
```

---

## Example Structure

Each workflow directory follows this structure:

```
<workflow>/
├── README.md              # Workflow-specific documentation
├── *.py                   # Main example scripts
├── scripts/               # Supporting scripts (optional)
│   ├── plot_*.py         # Visualization scripts
│   └── run_*.sh          # Shell automation scripts
└── results/              # Example outputs
    ├── *.json            # Metrics and metadata
    ├── *.tsv             # Tabular results
    └── plots/            # Visualizations (PNG/PDF)
```

---

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

### Workflow 2: Joint Genotyping + QC

```bash
# 1. Combine GVCFs
hvantk hgc gvcf-combine -g gvcfs/ -o cohort.vds

# 2. Convert to MatrixTable
hvantk hgc vds2mt -i cohort.vds -o cohort.mt --adjust-genotypes

# 3. Run QC workflow
python examples/hgc/hgc_qc_example.py --input cohort.mt
```

### Workflow 3: PSROC Analysis

```bash
# 1. Build ClinVar table
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht

# 2. Build dbNSFP table
hvantk mktable dbnsfp --raw-input dbnsfp.tsv.bgz --output-ht dbnsfp.ht

# 3. Run PSROC
hvantk psroc \
  --genes-file genes.txt \
  --clinvar-ht clinvar.ht \
  --dbnsfp-ht dbnsfp.ht \
  --scores CADD_phred REVEL_score \
  --output-dir results/
```

---

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

### Adapting Python Scripts

Most examples accept command-line arguments:

```bash
# Check available options
python examples/hgc/hgc_qc_example.py --help

# Run with custom input
python examples/hgc/hgc_qc_example.py \
  --input /path/to/cohort.mt \
  --output /path/to/results/
```

---

## Test Data

For quick testing, use the test datasets in `hvantk/tests/testdata/`:

```bash
# Example with test data
hvantk mkmatrix ucsc \
  -e hvantk/tests/testdata/raw/ucsc/exprMatrix.test.tsv.bgz \
  -m hvantk/tests/testdata/raw/ucsc/meta.test.tsv \
  -o /tmp/test_output.mt
```

**Test datasets:**
- `testdata/psroc/` - Synthetic PSROC data (100 variants, 4 scores)
- `testdata/raw/ucsc/` - UCSC Cell Browser test data
- `testdata/raw/tsv/` - Small TSV datasets
- Additional test data throughout `testdata/`

---

## Tips and Best Practices

1. **Start small**: Test workflows with small datasets first
2. **Check outputs**: Verify results before scaling to production data
3. **Use recipes**: Automate repetitive tasks with recipe files
4. **Save recipes**: Keep recipe files for reproducibility
5. **Monitor resources**: Check memory/CPU usage for large datasets
6. **Read READMEs**: Each workflow has specific tips in its README

---

## Troubleshooting

### Common Issues

**Import errors:**
```bash
# Ensure environment is activated
poetry shell

# Reinstall if needed
poetry install
```

**Out of memory:**
```python
# Increase Spark memory in Hail init
import hail as hl
hl.init(default_reference='GRCh38', driver_memory='16g', executor_memory='8g')
```

**File not found:**
```bash
# Use absolute paths in recipes
{
  "input": "/full/path/to/file.vcf.bgz"
}
```

**Hail not initialized:**
```python
# Initialize Hail before using hvantk functions
from hvantk.core.hail_context import init_hail
init_hail()
```

### Getting Help

- **Documentation**: See [`docs/`](../docs/)
- **Usage Guide**: [`docs/library/usage.md`](../docs/library/usage.md)
- **Architecture**: [`docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md)
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)

---

## Contributing Examples

Have a useful example? Please contribute!

1. Create a new workflow directory: `examples/<workflow>/`
2. Add your example scripts and README
3. Include sample data or reference test data in `hvantk/tests/testdata/`
4. Update this README with a link
5. Submit a pull request

See [`CONTRIBUTING.md`](../CONTRIBUTING.md) for guidelines.

---

## Additional Resources

### Documentation

- [Library Usage Guide](../docs/library/usage.md)
- [Annotation Sources](../docs/library/annotation-sources.md)
- [HGC Documentation](../docs/tools/hgc.md)
- [PSROC Documentation](../docs/tools/psroc.md)
- [Architecture Overview](../docs/ARCHITECTURE.md)

### External Resources

- [Hail Documentation](https://hail.is/docs/0.2/)
- [ClinVar Database](https://www.ncbi.nlm.nih.gov/clinvar/)
- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [UCSC Cell Browser](https://cells.ucsc.edu/)

---

## Directory Contents

```
examples/
├── README.md                    # This file
│
├── hgc/                         # HGC joint genotyping workflow
│   ├── README.md
│   ├── hgc_qc_example.py
│   ├── hgc_cpu_scaling_benchmark.py
│   ├── hgc_scalability_benchmark.py
│   ├── test_benchmark_setup.py
│   ├── scripts/                 # Supporting scripts
│   └── results/                 # Example outputs
│
├── psroc/                       # PSROC ROC analysis workflow
│   ├── README.md
│   ├── run_psroc_example.py
│   └── results/                 # Example outputs
│       ├── *.json              # Metrics
│       ├── *.tsv               # Annotated variants
│       └── plots/              # ROC curves and dashboards
│
├── clinvar/                     # ClinVar data streaming
│   ├── README.md
│   └── clinvar_streamer_example.py
│
└── recipes/                     # Recipe templates
    ├── tables.example.json
    ├── matrices.example.json
    ├── cptac.example.json
    └── heart_rare_variants.yaml
```

---

## License

See [`LICENSE`](../LICENSE) for details.
