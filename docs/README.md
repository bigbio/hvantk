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

### 📋 Planning & Development (`planning/`)

Documentation for developers and contributors:

- **[Developer Guide](../local/docs/DEVELOPING.md)** - Development workflow and contracts
- **[Streamers and Recipes](../local/docs/STREAMERS_AND_RECIPES.md)** - Streamer interface and JSON/YAML recipe format
- **[Data Catalog](../local/docs/DATA_CATALOG.md)** - Dataset registry, versioning, and hosting strategy

### 📊 Registry (`registry/`)

Dataset registry and web interface for browsing available datasets.

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

### For Developers
- [Contributing](../local/docs/DEVELOPING.md#contributing) - How to contribute
- [Add a new data source](../local/docs/DEVELOPING.md#development-cycle) - Add new annotations
- [Create custom streamers](../local/docs/STREAMERS_AND_RECIPES.md#creating-custom-streamers) - Build transformers
- [Register datasets](../local/docs/DATA_CATALOG.md#registering-a-new-dataset) - Add to catalog

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

See the [Developer Guide](../local/docs/DEVELOPING.md) for more information.
