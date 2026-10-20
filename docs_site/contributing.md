# Contributing to hvantk

Thank you for your interest in contributing to hvantk. This is the full
development guide. [`CONTRIBUTING.md`](https://github.com/bigbio/hvantk/blob/main/CONTRIBUTING.md)
at the repo root is a short entry point that links back into the sections below.

## Setup

```bash
git clone https://github.com/YOUR_USERNAME/hvantk
cd hvantk
poetry install
```

Activate the environment (Poetry 2.0+):

```bash
eval "$(poetry env activate)"
```

Verify the setup:

```bash
pytest -q
```

## Branch and PR flow

This project integrates on `dev`, not `main`. Branch from `dev`, and open your
pull request against `dev` — PRs are squash-merged there, and `dev` is
promoted to `main` separately at release time.

```bash
git checkout -b feature/your-feature-name dev
```

Use a prefix that matches the change: `feature/` for new functionality,
`fix/` for bug fixes, `docs/` for documentation-only changes.

Commit messages follow the conventional-commits style used throughout the
history: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`.

```bash
git push origin feature/your-feature-name
```

Open the PR on GitHub with a clear title, a description of what changed and
why, and a reference to any related issue. All submissions are reviewed
before merge; automated checks (tests, lint, format, docs build) must pass.

## Code style

- Follow PEP 8. Use type hints and docstrings on public functions and classes.
- Format with `ruff format .`. CI checks formatting with `ruff format --check .`.
- Lint with `ruff check . --select=E9,F63,F7,F82` — this is CI's blocking gate.
  CI also runs an advisory, non-blocking pass with the preview rule set:
  `ruff check . --exit-zero --statistics --preview`.
- Ruff configuration lives in `pyproject.toml` (`[tool.ruff]`, `[tool.ruff.lint]`,
  `[tool.ruff.format]`).
- One-time: `git config blame.ignoreRevsFile .git-blame-ignore-revs`, so
  `git blame` skips the bulk ruff-format reformat commit.

## Tests

```bash
pytest -q
```

`pytest.ini` sets `testpaths` to both `hvantk/tests` and `hvantk/skills`, so a
plain `pytest` run already collects every plugin's own tests alongside the
top-level suite. The default run deselects slower/optional categories via
`addopts`.

| Marker | Meaning | In default run? |
|---|---|---|
| `slow` | Long-running tests | No |
| `network` | Requires network I/O | No |
| `hail` | Requires Hail/Spark (auto-applies the session-scoped `hail_session` fixture) | No |
| `integration` | Integration tests | No |
| `round_trip(dataset)` | A plugin dataset's round-trip test, tagged with `"<provider>:<dataset>"` | Yes |
| `round_trip_contract` | The slow per-manifest round-trip collection check | No |

Common invocations:

```bash
pytest                          # everything, including Hail/network/slow tests
pytest -m hail                  # Hail-dependent tests only
pytest -m "not slow"            # exclude slow tests
pytest -m network                # network-dependent tests only
pytest hvantk/tests/test_file_utils.py   # one file
pytest -v                        # verbose output
```

CI additionally runs `pytest --require-round-trips`, which fails the run if a
`@pytest.mark.round_trip(...)` test ends up skipped or `xfail`ed for real (a
failure mode plain collection cannot see). Run it locally with the same flag
when you touch a plugin's round-trip test.

There is no coverage tooling wired up: `pytest-cov` is not a project
dependency, so `pytest --cov=hvantk` fails unless you install it yourself.

## Documentation

- Source lives under `docs_site/`; `mkdocs.yml` owns the navigation.
- Wrap every `hvantk ...` invocation you document in a ` ```bash ` fence.
  `hvantk/tests/test_docs_claims.py` checks that the flags you document on a
  fenced `hvantk` command actually exist on that command — it does not check
  defaults, semantics, Python snippets, or commands in other fence languages.
  Avoid `<placeholder>` tokens in the command path itself (one in the command
  path disables the check for that invocation).
- Before opening a docs PR, build with the same strict settings CI uses —
  it catches a broken internal link or anchor that a normal build only
  warns about, and keeps the output out of the repo root (`site/` is not
  gitignored):
  ```bash
  mkdocs build --strict --site-dir /tmp/site-check
  ```
- Renaming a heading changes its anchor. Before renaming one, check who
  links to it: `git grep -n "<page>.md#<anchor>"`.
- Update `README.md` too when a change is user-facing enough to belong there.

## Adding a new data provider plugin

Each external data source ships as a self-contained plugin under
`hvantk/skills/<provider>/`. The HGNC plugin (`hvantk/skills/hgnc/`) is the
canonical reference. Step by step:

1. Create the folder layout. Single-dataset providers:
   `hvantk/skills/<provider>/{plugin.yaml,builder.py,cli.py,drift_probe.py,SKILL.md,tests/}`.
   Multi-dataset providers: add one subfolder per dataset and a `shared/`
   folder for code reused across them.
2. Write `plugin.yaml` with `api_version: 2`, a `source.catalog_ref` (a bare
   provider identifier, e.g. `catalog_ref: hgnc`), and one `datasets:` entry
   per dataset declaring `builder`, `drift_probe`, optional
   `lifecycle.{download,parse}`, and the `tests:` block (fixture,
   schema/row snapshots, drift fingerprint, command). `command` must be
   exactly `pytest hvantk/skills/<this provider>/...`, with `-m hail`
   appended only when `backend: hail` and no other flags otherwise — that is
   the only shape CI's `plugin-contract` job and the default `pytest`
   selection actually run, and `hvantk/tests/test_plugin_contract_artifacts.py`
   checks it. If the builder can stamp more than one output schema, declare
   the default as `schema_id` and the others under `schema_ids`; a build
   that stamps any other ID is rejected.
3. Implement the plugin builder `build_<source>(parsed_input, ctx, **params)`
   in `builder.py`, returning a dataset (`AnnotationTable` / `ExpressionMatrix`
   / `VariantMatrix` / `GeneSet`) stamped via `ctx.provenance(schema_id=…)`.
   To record options that change the dataset's contents (a filter, a
   threshold), also pass `build_parameters={…}`. It is optional and must be
   JSON-serialisable; a non-JSON value raises when stamped. Generic Hail
   helpers (`create_table_base`, `cleanup_temp_file`) live in
   `hvantk/core/utils/hail_helpers.py`.
4. Implement the downloader in `cli.py` (Click command + a `download_dataset`
   function the loader can wire to `lifecycle.download`). Not every source
   needs one — see
   [Downloader or documented acquisition](guide/data-sources.md#downloader-or-documented-acquisition)
   for when to build a downloader versus document manual acquisition.
5. Implement the drift probe — a zero-arg function returning the fingerprint
   dict described in `hvantk/skills/_conventions/SKILL.md` § 12. Capture the
   expected fingerprint at `tests/drift_fingerprint.json` with
   `hvantk drift --regenerate <provider:dataset>` (not by hand) and commit
   it. A first capture records no rebuild-ledger row; a later hand-run
   regenerate takes `--accepted-in "PR #N"` and records one when upstream
   changed (`_conventions` § 12). See
   [Drift-fingerprint maintenance](#drift-fingerprint-maintenance) below.
6. Add the round-trip test under `tests/`, tag it
   `@pytest.mark.round_trip("<provider>:<dataset>")` (for a round trip
   parametrized over several datasets, put the mark on each case's
   `pytest.param(..., marks=pytest.mark.round_trip(...))` instead of stacking
   every id on the function — otherwise every case carries every dataset's
   id), and mark it `@pytest.mark.hail` too when the dataset's `backend` is
   `hail`. Then snapshot the schema and a small set of rows.
   `hvantk/tests/test_plugin_contract_artifacts.py` replays the manifest's
   `tests.command` in a subprocess and asserts it selects exactly that one
   marked item, clean of a skip/xfail mark, and hail-marked when required;
   run it locally with `-m round_trip_contract` (it is slow, so `pytest.ini`
   deselects it by default — CI's `plugin-contract` job runs it on every PR).
   `--require-round-trips` additionally fails a run if the round trip ends
   skipped or xfailed for real, which is the one failure mode (e.g. a
   `pytest.skip()` call inside the test body) collection cannot see.
7. Write `SKILL.md` following the nine-section template in
   `_conventions/SKILL.md` § 2.
8. Verify: `hvantk plugins validate hvantk/skills/<provider>/plugin.yaml`
   should accept the manifest; `hvantk plugins describe <provider>` should
   list the new dataset without a `FAILED TO LOAD` mark, and
   `hvantk plugins errors` should report `(no load errors)` (it exits 1 when
   anything failed to load, while `describe` still exits 0 when only one of
   the provider's datasets failed); the round-trip test must pass;
   `hvantk drift <provider:dataset>` must match the committed fingerprint.

Read `hvantk/skills/_conventions/SKILL.md` for the full contract (registry
keys, drift-probe shape, lifecycle stages, validation paths, hard
guardrails). See [Architecture](architecture.md#3-plugin-contract-adding-a-data-source)
for how the pieces fit together at runtime.

## Plugin SKILL.md maintenance

`hvantk/skills/` contains agent-readable methodology for maintaining stable
resources. Each per-resource `SKILL.md` is a design contract (source
identity, output schema, builder + CLI + registry wiring, validation
contract) that the round-trip test pins down. Shared conventions live in
`hvantk/skills/_conventions/SKILL.md` — read that first before any
per-resource skill.

**When you change builder code**: if a PR changes a function or path named
in any `SKILL.md` — including the `_conventions/` skill — update the
affected skill in the same PR. A PR with obvious skill drift will be
rejected in review.

**Skill status**: each `SKILL.md` carries a `status` frontmatter field:

- `provisional` — work in progress, not yet validated against the
  round-trip + update tests.
- `stable` — actively maintained, agent can rely on it.
- `deprecated` — superseded by a code change; do not use.

Do not invoke an agent against a `deprecated` skill.

**Skill-driven changes**: skill-assisted code changes (where an agent
generates the builder) are reviewed identically to hand-written code. Never
auto-merge a skill-generated PR. Always run the full test suite, especially
the affected plugin's `tests/test_builder.py`.

## Snapshot regeneration

When a builder change legitimately changes output:

```bash
poetry run pytest hvantk/skills/<provider>/tests/test_builder.py --regenerate-snapshots
# For Hail-marked builders only:
poetry run pytest hvantk/skills/<provider>/tests/test_builder.py -m hail --regenerate-snapshots
```

Review the diff in `hvantk/skills/<provider>/tests/snapshots/` before
committing. The regenerated snapshot is the new ground truth. (A handful of
plugins — `dbnsfp`, `gevir`, `gnomad_metrics` — share a fixture under
`hvantk/tests/testdata/raw/` instead of their own plugin folder; see that
plugin's `SKILL.md` § 9 for its exact snapshot path.)

## Drift-fingerprint maintenance

A plugin's `tests/drift_fingerprint.json` is the committed baseline that
`hvantk drift <provider:dataset>` compares the live source against. Update
it with the CLI, never by hand:

```bash
hvantk drift --regenerate <provider:dataset>
```

A first capture records no rebuild-ledger row. A later regenerate, run
because the upstream source materially changed, should name the PR that
accepts the new fingerprint, which also records a rebuild-ledger row in
`hvantk/resources/drift_ledger.json`:

```bash
hvantk drift --regenerate <provider:dataset> --accepted-in "PR #123"
```

See `hvantk/skills/_conventions/SKILL.md` § 12 for the fingerprint's exact
shape and what each field means.

## Code review process

All submissions require review:

1. **Automated checks** must pass: tests, lint, format, and (for docs
   changes) the strict mkdocs build.
2. **Manual review** checks code quality and design, test coverage,
   documentation completeness, and alignment with project goals.
3. Address reviewer comments by pushing additional commits to the same
   branch, then request re-review.

## Questions and links

- **Usage**: see [Building datasets](guide/usage.md) and [Data sources](guide/data-sources.md)
- **Architecture**: see [Architecture](architecture.md)
- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Discussions**: [GitHub Discussions](https://github.com/bigbio/hvantk/discussions)

By contributing to hvantk, you agree that your contributions will be
licensed under the MIT License.
