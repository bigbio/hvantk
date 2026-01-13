[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)

# hvantk

**Hail-based toolkit for multi-omics variant annotation and analysis.**

`hvantk` is an annotation toolkit that uses [Hail](https://hail.is/) to annotate variants and genes with multiple omics data types. It enables integration of variant prediction scores, gene expression data, protein expression, and clinical annotations to improve genetic variant interpretation.

`hvantk` is a modular toolkit that uses [Apache Hail](https://hail.is/) to annotate and analyze variants, genes, proteins, and expression data from heterogeneous omics sources. The library enables multi-omics integration to improve the interpretation of genetic variants.

**Core Capabilities:**
- Variant annotations (ClinVar, dbNSFP, gnomAD, CCR scores)
- Gene annotations (Ensembl, GeVIR, gene constraints)
- Protein annotations (INSIDER protein-protein interactions)
- Expression data (bulk & single-cell RNA-seq from UCSC, GTEx)
- Joint genotyping workflows (GVCF combining, QC, format conversion)
- Recipe-based batch processing

- **Multiomics Integration**: Combine variant annotations, gene expression, and clinical data
- **Hail-Powered**: Leverage Hail's scalable genomic data processing
- **Modular Design**: Extensible framework for adding new data sources
- **Joint Genotyping**: HGC module for efficient GVCF combination and quality control
- **Recipe-Based Workflows**: JSON/YAML recipes for reproducible analyses
- **Multiple Formats**: Support for VCF, MatrixTable, and Hail Table formats

### Using Poetry (recommended)

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
poetry shell
```

### Using pip

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install -e .
```

**Prerequisites**: Python ≥3.10, Apache Hail

## Main Tools

### HGC: Joint Genotyping Pipeline

High-performance joint genotyping for large cohorts using Hail. Combines thousands of GVCF files with integrated quality control.

**Key features:**
- GVCF combination at scale
- VDS ↔ MatrixTable ↔ VCF format conversion
- Comprehensive QC metrics and visualization
- Professional HTML QC reports

**Quick example:**
```bash
# Combine GVCF files
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds

# Convert to MatrixTable and run QC
hvantk hgc vds2mt -i cohort.vds -o cohort.mt
hvantk hgc compute-qc -i cohort.mt -o cohort_qc.mt

# Generate QC report
hvantk hgc qc-report -i cohort_qc.mt -o qc_report.html
```

📖 **[Full HGC Documentation](docs/tools/hgc.md)**

### Annotation Tables: Build Custom Annotation Resources

Create Hail Tables from public annotation databases (ClinVar, gnomAD, Ensembl, etc.).

**Single table creation:**
```bash
# Build ClinVar annotation table
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht --ref-genome GRCh38

# Build Ensembl gene table
hvantk mktable ensembl-gene --raw-input biomart.tsv.bgz --output-ht ensembl.ht
```

**Batch processing with recipes:**
```bash
hvantk mktable-batch --recipe tables_recipe.json
```

📖 **[Annotation Tables Guide](docs/library/usage.md#1-build-a-single-annotation-table-ht)**

### Expression Matrices: Process Omics Data

Build Hail MatrixTables from bulk and single-cell expression data.

**Example:**
```bash
# Convert UCSC Cell Browser data to MatrixTable
hvantk mkmatrix ucsc -e expr.tsv.bgz -m metadata.tsv -o ucsc.mt

# Batch processing
hvantk mkmatrix-batch --recipe matrices_recipe.json
```

📖 **[Expression Data Guide](docs/library/usage.md#3-build-a-single-matrixtable-mt)**

### Data Downloaders

Download curated datasets directly from public repositories.

**Example:**
```bash
# Download UCSC Cell Browser dataset
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc
```

📖 **[Data Sources](docs/library/annotation-sources.md)**

## Quick Start Example

```bash
# 1. Download a dataset
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc

# 2. Convert to Hail MatrixTable
hvantk mkmatrix ucsc \
  -e data/ucsc/exprMatrix.tsv.bgz \
  -m data/ucsc/meta.tsv \
  -o data/ucsc/adultPancreas.mt

# 3. Build annotation tables via recipe
cat > recipe.json << EOF
{
  "tables": [
    {
      "name": "clinvar",
      "input": "/data/clinvar.vcf.bgz",
      "output": "/out/clinvar.ht",
      "params": {"reference_genome": "GRCh38"}
    }
  ]
}
EOF

hvantk mktable-batch --recipe recipe.json
```

## Documentation

- **[Architecture Overview](docs/ARCHITECTURE.md)** - Module organization and design patterns
- **[Usage Guide](docs/library/usage.md)** - Detailed usage examples and recipes
- **[Data Sources](docs/library/annotation-sources.md)** - Available annotation sources and download instructions
- **[HGC Tool](docs/tools/hgc.md)** - Joint genotyping and quality control
- **[API Reference](docs/ARCHITECTURE.md#extension-points)** - Extending hvantk with custom builders
- **[Full Documentation Index](docs/README.md)** - Complete documentation structure

## Citation

If you use hvantk in your research, please cite:

```bibtex
@software{hvantk2024,
  title = {hvantk: Hail-based toolkit for multi-omics variant annotation and analysis},
  author = {Perez-Riverol, Yasset and Audain, Enrique},
  year = {2024},
  url = {https://github.com/bigbio/hvantk}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for detailed information on:
- Development workflow and setup
- Adding new data sources
- Code style guidelines
- Testing requirements
- Pull request process

**Developer quick start:**
```bash
poetry install
pytest -q
hvantk --help
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Questions**: Open a discussion on GitHub
- **Documentation**: [docs/](docs/)

## Acknowledgments

- Built on [Apache Hail](https://hail.is/) for distributed genomic data processing
- Integrates data from ClinVar, gnomAD, Ensembl, UCSC, and other public resources
