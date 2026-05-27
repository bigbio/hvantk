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
   eval "$(poetry env activate)"
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
- Update relevant documentation in `docs_site/`
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

New data sources are added as plugins under `hvantk/skills/<provider>/`. The
canonical reference is `hvantk/skills/_conventions/SKILL.md` — start there
and copy from an existing plugin (clinvar, hgnc, msigdb are good models).

1. **Scaffold the plugin folder**
   - `hvantk/skills/<provider>/plugin.yaml` — manifest declaring the
     `Builder` / `DriftProbe` / optional `DownloadFn` + `ParseFn`,
     `artifact_type`, `schema_id`, and test fixture paths.
   - `hvantk/skills/<provider>/builder.py` — implement the Phase B
     contract: `build_<provider>_<dataset>(parsed_input, ctx, **params) -> Artifact`.
   - `hvantk/skills/<provider>/drift_probe.py` — return a dict the platform
     hashes into a `source_fingerprint`.
   - `hvantk/skills/<provider>/SKILL.md` — author-facing operational guide.

2. **Write the Phase B builder**
   ```python
   from hvantk.core.models import AnnotationTable

   def build_myprovider_dataset(parsed_input, ctx, *, **params):
       """Phase B builder — returns an AnnotationTable."""
       # ... import + transform ...
       return AnnotationTable.from_hail(
           ht, provenance=ctx.provenance(schema_id="myprovider-v1")
       )
   ```
   `BuildContext` (`ctx`) supplies plugin name, version, and source
   fingerprint; the platform stamps Provenance and validates the artifact
   type / schema_id against `plugin.yaml`. **The CLI `hvantk reprocess` is
   the only public build path** — there is no separate programmatic API.

3. **Create tests**
   - `hvantk/skills/<provider>/tests/test_builder.py` — snapshot
     round-trip test (use `phase_b_snapshot_adapter` from
     `hvantk/tests/_snapshot_utils.py`).
   - `hvantk/skills/<provider>/tests/testdata/` — minimal fixture.
   - `hvantk/skills/<provider>/tests/snapshots/` — generated via
     `pytest --regenerate-snapshots`.

4. **Update documentation**
   - Add to [Data Sources](guide/data-sources.md).
   - Add usage example to [Usage Guide](guide/usage.md) showing
     `hvantk reprocess <provider>:<dataset>`.

See [Architecture](architecture.md) for detailed information on the plugin
system and extension points.

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
