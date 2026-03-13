# hvantk Documentation

Welcome to the hvantk documentation! This directory contains comprehensive documentation for the hvantk library.

## Documentation Structure

The documentation is organized into the following sections:

### 📚 Library Documentation (`library/`)

Core library documentation for end users:

- **[Usage Guide](library/usage.md)** - Complete guide to using hvantk for building Hail Tables and MatrixTables
- **[Annotation Sources](library/annotation-sources.md)** - Catalog of available annotation sources and how to obtain them

### 🛠️ Tools Documentation (`tools/`)

Documentation for specific tools within hvantk:

- **[HGC (Hail-based Genotype Combiner)](tools/hgc.md)** - Joint genotyping workflows, VDS/MatrixTable conversion, and quality control
- **[PSROC (Pathogenicity Score ROC)](tools/psroc.md)** - Variant score evaluation using ClinVar truth labels
- **[EnrichEx](tools/enrichex.md)** - Gene set enrichment analysis (overlap and burden testing)
- **[Ancestry Inference](tools/ancestry.md)** - Genetic ancestry prediction using PCA and Random Forest classification

### 📋 Planning & Development

Documentation for developers and contributors:

- **[Architecture](ARCHITECTURE.md)** - System design, protocols, and extension points
- **[Contributing](../CONTRIBUTING.md)** - Development workflow and contribution guidelines

## Quick Links

### Getting Started
1. [Installation](../README.md#installation) - Install hvantk
2. [Usage Guide](library/usage.md) - Build your first Table/MatrixTable
3. [Examples](../examples/) - Ready-to-use examples and recipes

### Common Tasks
- [Download annotation data](library/annotation-sources.md) - Get source data
- [Build annotation tables](library/usage.md#1-build-a-single-annotation-table-ht) - Create Hail Tables
- [Create expression matrices](library/usage.md#3-build-a-single-matrixtable-mt) - Create MatrixTables
- [Joint genotyping with HGC](tools/hgc.md#quick-start) - Combine GVCF files
- [Ancestry inference](tools/ancestry.md#quick-start) - Predict genetic ancestry

### For Developers
- [Contributing](../CONTRIBUTING.md) - How to contribute
- [Add a new data source](ARCHITECTURE.md#adding-a-new-data-source) - Add new annotations
- [Architecture](ARCHITECTURE.md) - System design and extension points

## Documentation Conventions

- **Code examples**: All code examples are tested and ready to use
- **File paths**: Paths are examples; adjust to your environment
- **Commands**: Shell commands assume you're in the hvantk root directory
- **Formats**: Recipes support both JSON and YAML (YAML requires PyYAML)

## Need Help?

- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Examples**: [Example scripts and recipes](../examples/)
- **Main README**: [Project overview and quick start](../README.md)

## Contributing to Documentation

Documentation contributions are welcome! Please:

1. Keep documentation up-to-date with code changes
2. Include practical examples
3. Use clear, concise language
4. Add links between related sections
5. Test all code examples

See the [Contributing Guide](../CONTRIBUTING.md) and [Architecture](ARCHITECTURE.md) for more information.
