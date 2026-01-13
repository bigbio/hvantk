[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)

# hvantk

**Hail-based multiomics variant annotation toolkit**

`hvantk` is a powerful annotation toolkit that uses [Hail](https://hail.is/) to annotate variants and genes with multiple omics data types. It enables integration of variant prediction scores, gene expression data, protein expression, and clinical annotations to improve genetic variant interpretation.

## ✨ Key Features

- **Multiomics Integration**: Combine variant annotations, gene expression, and clinical data
- **Hail-Powered**: Leverage Hail's scalable genomic data processing
- **Modular Design**: Extensible framework for adding new data sources
- **Joint Genotyping**: HGC module for efficient GVCF combination and quality control
- **Recipe-Based Workflows**: JSON/YAML recipes for reproducible analyses
- **Multiple Formats**: Support for VCF, MatrixTable, and Hail Table formats

## 📦 Installation

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
poetry shell
```

Prerequisites: [Poetry](https://python-poetry.org/) for dependency management.

## 🚀 Quick Start

### Build Annotation Tables

```bash
# ClinVar annotations
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht --ref-genome GRCh38

# Gene constraint metrics
hvantk mktable gnomad-metrics --raw-input gnomad.tsv.bgz --output-ht gnomad.ht
```

### Create Expression MatrixTables

```bash
# UCSC Cell Browser data
hvantk mkmatrix ucsc \
  --expression-matrix expr.tsv.bgz \
  --metadata meta.tsv \
  --output-mt ucsc.mt
```

### Joint Genotyping with HGC

```bash
# Combine GVCF files
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds

# Convert to MatrixTable
hvantk hgc vds2mt -i cohort.vds -o cohort.mt --adjust-genotypes

# Generate QC report
hvantk hgc qc-report -i cohort.mt -o qc_report.html
```

## 📚 Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- **[Usage Guide](docs/library/usage.md)** - Detailed examples and workflows
- **[HGC Tool](docs/tools/hgc.md)** - Joint genotyping and quality control
- **[Annotation Sources](docs/library/annotation-sources.md)** - Available data sources
- **[Developer Guide](docs/planning/DEVELOPING.md)** - Contributing and development
- **[Full Documentation Index](docs/README.md)** - Complete documentation structure

## 🧬 Supported Data Sources

### Variant Annotations
- **ClinVar** - Clinical significance annotations
- **dbNSFP** - Missense variant prediction scores
- **gnomAD** - Population allele frequencies and constraint metrics
- **INSIDER** - Protein-protein interaction sites

### Expression Data
- **UCSC Cell Browser** - Single-cell RNA-seq datasets
- **Expression Atlas** - Bulk RNA-seq across tissues and conditions
- **CPTAC** - Protein expression data

See [Annotation Sources](docs/library/annotation-sources.md) for complete list and download instructions.

## 🔧 Tools

### HGC: Hail-based Genotype Combiner

High-performance joint genotyping workflow for combining GVCF files:

- Scalable GVCF combination (1000s of samples)
- VDS ↔ MatrixTable ↔ VCF format conversion
- Comprehensive quality control and visualization
- Interactive HTML reports

[HGC Documentation](docs/tools/hgc.md) | [Examples](examples/hgc_qc_example.py)

## 💻 For Developers

```bash
# Run tests
pytest -q

# Explore CLI
hvantk --help

# Add a new data source
# 1. Define schema contract
# 2. Write builder function
# 3. Register in catalog
# 4. Add tests
```

See [Developer Guide](docs/planning/DEVELOPING.md) for detailed workflow.

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 🤝 Contributing

Contributions welcome! Please see our [Developer Guide](docs/planning/DEVELOPING.md) for:
- Development workflow
- Code style guidelines
- Testing requirements
- Pull request process

## 📧 Support

- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Documentation**: [docs/](docs/)
- **Examples**: [examples/](examples/)
