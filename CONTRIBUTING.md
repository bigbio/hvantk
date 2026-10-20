# Contributing to hvantk

The full guide — setup, branch and PR flow, code style, tests, documentation
rules, and the plugin checklist — lives at
[docs_site/contributing.md](docs_site/contributing.md). This page is a quick
entry point.

## Setup

```bash
poetry install
eval "$(poetry env activate)"
pytest -q
```

## Where PRs go

Branch from `dev` and open the PR against `dev`, not `main`
([details](docs_site/contributing.md#branch-and-pr-flow)).

## Before opening a PR

```bash
ruff format .
ruff check . --select=E9,F63,F7,F82
pytest -q
```

For a documentation change, also run:

```bash
mkdocs build --strict --site-dir /tmp/site-check
```

## More

- [Setup](docs_site/contributing.md#setup)
- [Branch and PR flow](docs_site/contributing.md#branch-and-pr-flow)
- [Code style](docs_site/contributing.md#code-style)
- [Tests](docs_site/contributing.md#tests)
- [Documentation](docs_site/contributing.md#documentation)
- [Adding a new data provider plugin](docs_site/contributing.md#adding-a-new-data-provider-plugin)
- [Plugin SKILL.md maintenance](docs_site/contributing.md#plugin-skillmd-maintenance)
- [Snapshot regeneration](docs_site/contributing.md#snapshot-regeneration)
- [Drift-fingerprint maintenance](docs_site/contributing.md#drift-fingerprint-maintenance)

By contributing to hvantk, you agree that your contributions will be
licensed under the MIT License.
