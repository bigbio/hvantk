# Contributing to hvantk

Thank you for your interest in contributing to hvantk! This guide will help you get started.

## Quick Start

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/hvantk
   cd hvantk
   ```

2. **Set up the development environment**
   ```bash
   poetry install
   poetry shell
   ```

3. **Run tests to verify setup**
   ```bash
   pytest -q
   ```

## Development Workflow

### 1. Create a Feature Branch

```bash
git checkout -b feature/your-feature-name
```

Use descriptive branch names:
- `feature/add-new-datasource` for new features
- `fix/bug-description` for bug fixes
- `docs/update-readme` for documentation changes

### 2. Make Your Changes

Follow these guidelines:

#### Code Style
- Follow PEP 8 for Python code
- Use type hints where applicable
- Write clear docstrings for public functions and classes
- Keep functions focused and modular

#### Testing
- Add tests for new functionality
- Update existing tests if behavior changes
- Ensure all tests pass before submitting
- Use fixtures in `hvantk/tests/testdata` for test data

#### Documentation
- Update relevant documentation in `docs/`
- Add docstrings to new functions and classes
- Include usage examples where appropriate
- Update README.md if adding major features

### 3. Test Your Changes

```bash
# Run all tests
pytest

# Run specific test file
pytest hvantk/tests/test_file.py

# Run with verbose output
pytest -v

# Run with coverage
pytest --cov=hvantk
```

### 4. Commit Your Changes

Write clear, descriptive commit messages:

```bash
git add .
git commit -m "Add support for new annotation source

- Implement builder for XYZ database
- Add tests for XYZ builder
- Update documentation
"
```

Follow conventional commits format:
- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation changes
- `test:` for test additions/changes
- `refactor:` for code refactoring

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub with:
- Clear title describing the change
- Description of what changed and why
- Reference to related issues (if any)
- Screenshots (if UI changes)

## Adding a New Data Source

Follow this workflow when adding a new annotation source:

1. **Define the data product contract**
   - Specify Table/MatrixTable schema
   - Document key fields and data types
   - Include metadata requirements

2. **Write the builder function**
   ```python
   def build_my_datasource(raw_input: str, output_ht: str, **kwargs):
       """Build Hail Table from raw data source.

       Args:
           raw_input: Path to raw input file
           output_ht: Path for output Hail Table
           **kwargs: Additional parameters
       """
       # Implementation
   ```

3. **Register in the catalog**
   - Add to dataset registry
   - Include provenance and version info
   - Add checksums for data integrity

4. **Create tests**
   - Add test data to `hvantk/tests/testdata`
   - Write unit tests for the builder
   - Test with various input scenarios

5. **Update documentation**
   - Add to [Annotation Sources](guide/annotation-sources.md)
   - Add usage example to [Usage Guide](guide/usage.md)

See [Architecture](architecture.md) for detailed information on system design and extension points.

## Code Review Process

All submissions require review:

1. **Automated checks** must pass:
   - Tests must pass
   - Code must follow style guidelines
   - No merge conflicts

2. **Manual review** will check:
   - Code quality and design
   - Test coverage
   - Documentation completeness
   - Alignment with project goals

3. **Feedback and iteration**:
   - Address reviewer comments
   - Push additional commits to same branch
   - Request re-review when ready

## Community Guidelines

- Be respectful and inclusive
- Provide constructive feedback
- Ask questions if something is unclear
- Help others when you can

## Questions?

- **Documentation**: See the [guide](guide/) section
- **Architecture**: See [architecture.md](architecture.md)
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Discussions**: [GitHub Discussions](https://github.com/bigbio/hvantk/discussions)

## License

By contributing to hvantk, you agree that your contributions will be licensed under the MIT License.
