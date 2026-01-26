[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)

# hvantk

**Hail-based toolkit for multiomics variant annotation and analysis.**

`hvantk` is a modular toolkit that uses [Hail](https://hail.is/) to annotate and analyze variants, genes, proteins, and expression data from heterogeneous omics sources. The library enables multiomics integration to improve the interpretation of genetic variants.

**Core Capabilities:**
- Variant annotations (ClinVar, dbNSFP, gnomAD, CCR scores)
- Gene annotations (Ensembl, GeVIR, gene constraints)
- Protein annotations (INSIDER protein-protein interactions)
- Expression data (bulk & single-cell RNA-seq from UCSC, GTEx)
- Joint genotyping workflows (GVCF combining, QC, format conversion)
- Recipe-based batch processing

## Installation

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

**Prerequisites**: Python ≥3.10, Hail

## Core Workflows

### HGC: Joint Genotyping Pipeline

High-performance joint genotyping for large cohorts. Combines thousands of GVCF files with integrated QC.

```bash
# End-to-end pipeline
hvantk hgc pipeline -i /data/gvcfs -o /output

# Or run individual steps
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds
hvantk hgc vds2mt -i cohort.vds -o cohort.mt
hvantk hgc qc-report -i cohort.mt -o qc_report.html
```

📖 **[Full HGC Documentation](docs/tools/hgc.md)**

### PSROC: Variant Score Evaluation

Evaluate pathogenicity prediction scores (CADD, REVEL, MetaLR) using ClinVar truth labels. Generate ROC curves and performance metrics.

```bash
# Run ROC analysis
hvantk psroc \
  --genes-file genes.txt \
  --clinvar-ht clinvar.ht \
  --dbnsfp-ht dbnsfp.ht \
  --scores "CADD_phred,REVEL_score" \
  --output-dir results/
```

📖 **[PSROC Documentation](docs/tools/psroc.md)** | **[Example](examples/psroc/)**

### Annotation Tables

Create Hail Tables from public databases (ClinVar, gnomAD, Ensembl).

```bash
# Single table
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht

# Batch processing
hvantk mktable-batch --recipe tables_recipe.json
```

📖 **[Tables Guide](docs/library/usage.md#1-build-a-single-annotation-table-ht)**

### Expression Matrices

Build Hail MatrixTables from bulk and single-cell expression data.

```bash
# UCSC Cell Browser data
hvantk mkmatrix ucsc -e expr.tsv.bgz -m metadata.tsv -o ucsc.mt

# Batch processing
hvantk mkmatrix-batch --recipe matrices_recipe.json
```

📖 **[Expression Guide](docs/library/usage.md#3-build-a-single-matrixtable-mt)**

### Data Downloaders

Download curated datasets from public repositories.

```bash
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc
```

📖 **[Data Sources](docs/library/annotation-sources.md)**

## Quick Start

```bash
# Download and process expression data
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc
hvantk mkmatrix ucsc -e data/ucsc/exprMatrix.tsv.bgz -m data/ucsc/meta.tsv -o data/ucsc/adultPancreas.mt

# Build annotation tables
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht --ref-genome GRCh38

# Or use batch processing with recipes (see examples/recipes/)
hvantk mktable-batch --recipe recipe.json
```

## Documentation

- **[Usage Guide](docs/library/usage.md)** - Examples and recipes
- **[HGC Tool](docs/tools/hgc.md)** - Joint genotyping pipeline
- **[PSROC Tool](docs/tools/psroc.md)** - Variant score evaluation
- **[Data Sources](docs/library/annotation-sources.md)** - Available annotations
- **[Architecture](docs/ARCHITECTURE.md)** - Design and extension points
- **[Full Index](docs/README.md)** - Complete documentation

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

- Built on [Hail](https://hail.is/) for distributed genomic data processing
- Integrates data from ClinVar, gnomAD, Ensembl, UCSC, and other public resources
