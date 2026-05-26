# hvantk

**Hail-based toolkit for multiomics variant annotation and analysis.**

`hvantk` is a modular toolkit that uses [Hail](https://hail.is/) to annotate and analyze variants, genes, proteins, and expression data from heterogeneous omics sources. The library enables multiomics integration to improve the interpretation of genetic variants.

## Core Capabilities

- **Variant annotations** — ClinVar, dbNSFP, gnomAD, CCR scores
- **Gene annotations** — Ensembl, GeVIR, gene constraints
- **Protein annotations** — INSIDER protein-protein interactions
- **Expression data** — bulk & single-cell RNA-seq from UCSC, GTEx
- **Joint genotyping** — GVCF combining, QC, format conversion ([HGC](tools/hgc.md))
- **Ancestry inference** — PCA + Random Forest classification ([Ancestry](tools/ancestry.md))
- **Enrichment analysis** — overlap and burden testing ([EnrichEx](tools/enrichex.md))
- **Score evaluation** — pathogenicity score ROC analysis ([PS-ROC](tools/psroc.md))

## Get Started

- [Installation](getting-started/installation.md) — set up hvantk with Poetry or pip
- [Quick Start](getting-started/quickstart.md) — core workflows and first commands

## Documentation

- [Building Tables](guide/usage.md) — usage examples
- [Data Sources](guide/data-sources.md) — available data sources and acquisition
- [Architecture](architecture.md) — design and extension points
- [Examples](examples/index.md) — end-to-end workflow examples
- [Contributing](contributing.md) — development workflow

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

## License

This project is licensed under the MIT License — see the [LICENSE](https://github.com/bigbio/hvantk/blob/main/LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Questions**: Open a discussion on GitHub

## Acknowledgments

- Built on [Hail](https://hail.is/) for distributed genomic data processing
- Integrates data from ClinVar, gnomAD, Ensembl, UCSC, and other public resources
