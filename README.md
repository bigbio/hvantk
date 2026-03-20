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
- Ancestry inference (PCA + Random Forest classification)
- Recipe-based batch processing

## Installation

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

### Using pip

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install -e .
```

**Prerequisites**: Python ≥3.10, Hail

## Core Workflows

```bash
# Joint genotyping pipeline
hvantk hgc pipeline -i /data/gvcfs -o /output

# Pathogenicity score ROC analysis
hvantk psroc --genes BRCA1,BRCA2 --clinvar-ht clinvar.ht --dbnsfp-ht dbnsfp.ht --scores "CADD_phred,REVEL_score" -o results/

# Gene set enrichment (overlap + burden)
hvantk enrichex overlap -g genes.txt -s gene_sets.json -o overlap.tsv
hvantk enrichex burden -m cohort.mt -p phenotypes.ht -s gene_sets.json -o burden.tsv

# Ancestry inference
hvantk ancestry-inference -q cohort.mt -r 1kg_reference.mt --ancestry-col super_pop -o ancestry.ht

# Annotation tables and expression matrices
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht
hvantk mkmatrix ucsc -e expr.tsv.bgz -m metadata.tsv -o ucsc.mt
hvantk mktable-batch --recipe recipe.json
```

See the [Quick Start Guide](docs_site/getting-started/quickstart.md) for detailed walkthroughs.

## Documentation

**Browse the full documentation site:** [https://bigbio.github.io/hvantk](https://bigbio.github.io/hvantk)

Or read the source markdown directly:

- **[Usage Guide](docs_site/guide/usage.md)** - Examples and recipes
- **[HGC Tool](docs_site/tools/hgc.md)** - Joint genotyping pipeline
- **[PSROC Tool](docs_site/tools/psroc.md)** - Variant score evaluation
- **[EnrichEx Tool](docs_site/tools/enrichex.md)** - Gene set enrichment analysis
- **[Ancestry Tool](docs_site/tools/ancestry.md)** - Genetic ancestry inference
- **[Data Sources](docs_site/guide/data-sources.md)** - Available annotations
- **[Architecture](docs_site/architecture.md)** - Design and extension points

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
- **Documentation**: [docs_site/](docs_site/)

## Acknowledgments

- Built on [Hail](https://hail.is/) for distributed genomic data processing
- Integrates data from ClinVar, gnomAD, Ensembl, UCSC, and other public resources
