[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://bigbio.github.io/hvantk)

# hvantk

**Hail-based toolkit for multiomics variant annotation and analysis.**

`hvantk` is a modular toolkit that uses [Hail](https://hail.is/) to annotate and analyze variants, genes, proteins, and expression data from heterogeneous omics sources. The library enables multiomics integration to improve the interpretation of genetic variants.

## Installation

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

**Prerequisites**: Python >=3.10, Hail

Verify it works:

```bash
hvantk utils check-install
hvantk --help
```

## Toolkit

| Tool | Description | Command | Docs |
|------|-------------|---------|------|
| **Downloads** | Acquire external datasets (ClinVar, ClinGen, HGNC, etc.) | `hvantk download <source>` | [Data Sources](docs_site/guide/data-sources.md) |
| **Annotation builders** | Variant, gene, and protein tables (ClinVar, dbNSFP, gnomAD, Ensembl, HGNC, INSIDER, CCR) | `hvantk mktable <source>` | [Usage Guide](docs_site/guide/usage.md) |
| **Expression builders** | Bulk and single-cell matrices (UCSC, GTEx, Expression Atlas, CPTAC) | `hvantk mkmatrix <source>` | [Usage Guide](docs_site/guide/usage.md) |
| **Batch recipes** | Recipe-based batch processing for tables and matrices | `hvantk mktable-batch` | [Recipes](docs_site/examples/recipes.md) |
| **HGC** | Joint genotyping pipeline (GVCF combining, QC, format conversion) | `hvantk hgc` | [HGC](docs_site/tools/hgc.md) |
| **Ancestry** | Ancestry inference (PCA + Random Forest classification) | `hvantk ancestry-inference` | [Ancestry](docs_site/tools/ancestry.md) |
| **QTL Cascade** | Molecular QTL integration (eQTL + pQTL cascade, colocalization ABF) | `hvantk qtlcascade` | [QTL Cascade](docs_site/tools/qtlcascade.md) |
| **EnrichEx** | Gene set enrichment (overlap testing + rare variant burden) | `hvantk enrichex` | [EnrichEx](docs_site/tools/enrichex.md) |
| **PS-ROC** | Pathogenicity score ROC evaluation against ClinVar labels | `hvantk psroc` | [PS-ROC](docs_site/tools/psroc.md) |
| **PTM** | Post-translational modification variant classification | `hvantk ptm` | [PTM](docs_site/tools/ptm.md) |
| **Expression** | Expression analysis (summarize, marker extraction) | `hvantk expression` | [Usage Guide](docs_site/guide/usage.md) |

## Documentation

**Full docs site:** [https://bigbio.github.io/hvantk](https://bigbio.github.io/hvantk)

- [Data Sources](docs_site/guide/data-sources.md) -- Available annotations and how to acquire them
- [Examples](docs_site/examples/index.md) -- Tutorials and walkthroughs for each tool
- [Architecture](docs_site/architecture.md) -- Design patterns and extension points

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

See [CONTRIBUTING.md](CONTRIBUTING.md) for development workflow, code style, and testing requirements.

```bash
poetry install
pytest -q
hvantk --help
```

## License

MIT License - see [LICENSE](LICENSE).

## Support

- [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- [Documentation](https://bigbio.github.io/hvantk)
