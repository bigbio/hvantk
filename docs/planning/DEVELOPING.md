# Developer Guide

This guide covers the development workflow and contracts for contributing to hvantk.

## Development Workflow

### Setting Up Development Environment

1. Clone the repository:
   ```bash
   git clone https://github.com/bigbio/hvantk
   cd hvantk
   ```

2. Install dependencies using Poetry:
   ```bash
   poetry install
   poetry shell
   ```

3. Run tests to verify setup:
   ```bash
   pytest -q
   ```

### Development Cycle

1. **Define a data product contract**: Specify the Table/MatrixTable schema and metadata
2. **Write a downloader** (optional): Create a data fetcher if the source requires it
3. **Build a builder**: Create a function that outputs a Hail Table/MatrixTable
4. **Register the dataset**: Add to the data catalog with provenance, versions, and checksums
5. **Create streamers**: Write transformers that process the data
6. **Compose a recipe**: Define how to answer biological questions using your data
7. **Add tests**: Create test cases using fixtures in `hvantk/tests/testdata`

### Code Style

- Follow PEP 8 for Python code
- Use type hints where applicable
- Write docstrings for public functions and classes
- Keep functions focused and modular

### Testing

Run the test suite:
```bash
# All tests
pytest

# Specific module
pytest hvantk/tests/hgc/

# With verbose output
pytest -v
```

### Contributing

1. Create a feature branch from `main`
2. Make your changes
3. Add/update tests
4. Update documentation
5. Submit a pull request

## Data Product Contracts

### Table Schema Contract

When adding a new annotation source:

- Define clear key fields (e.g., locus/alleles for variants, gene_id for genes)
- Document expected data types
- Specify required vs optional fields
- Include provenance metadata

### MatrixTable Contract

For expression or omics data:

- Row keys: genes, variants, or features
- Column keys: samples, cells, or conditions
- Entry fields: measurements (expression, dosage, etc.)
- Column metadata: sample/cell annotations
- Row metadata: feature annotations

## See Also

- [Data Catalog](DATA_CATALOG.md) - Dataset registry and hosting strategy
- [Streamers and Recipes](STREAMERS_AND_RECIPES.md) - Transformer interface and recipe format
- [Usage Guide](../library/usage.md) - End-user documentation
