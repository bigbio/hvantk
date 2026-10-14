# Contributing to hvantk

See the full contributing guide at [docs_site/contributing.md](docs_site/contributing.md).

## Adding a new data provider plugin

Each external data source ships as a self-contained plugin under `hvantk/skills/<provider>/`. The HGNC plugin (`hvantk/skills/hgnc/`) is the canonical reference. Step by step:

1. Create the folder layout. Single-dataset providers: `hvantk/skills/<provider>/{plugin.yaml,builder.py,cli.py,drift_probe.py,SKILL.md,tests/}`. Multi-dataset providers: add one subfolder per dataset and a `shared/` folder for code reused across them.
2. Write `plugin.yaml` with `api_version: 2`, a `source.catalog_ref` (a bare provider identifier, e.g. `catalog_ref: hgnc`), and one `datasets:` entry per dataset declaring `builder`, `drift_probe`, optional `lifecycle.{download,parse}`, and the `tests:` block (fixture, schema/row snapshots, drift fingerprint, command). `command` must be exactly `pytest hvantk/skills/<this provider>/...`, with `-m hail` appended only when `backend: hail` and no other flags otherwise — that is the only shape CI's `plugin-contract` job and the default `pytest` selection actually run, and `hvantk/tests/test_plugin_contract_artifacts.py` checks it. If the builder can stamp more than one output schema, declare the default as `schema_id` and the others under `schema_ids`; a build that stamps any other ID is rejected.
3. Implement the plugin builder `build_<source>(parsed_input, ctx, **params)` in `builder.py`, returning a typed artifact (`AnnotationTable` / `ExpressionMatrix` / `VariantMatrix` / `GeneSet`) stamped via `ctx.provenance(schema_id=…)`. To record options that change the artifact's contents (a filter, a threshold), also pass `build_parameters={…}`. It is optional and must be JSON-serialisable; a non-JSON value raises when stamped. Generic Hail helpers (`create_table_base`, `cleanup_temp_file`) live in `hvantk/core/utils/hail_helpers.py`.
4. Implement the downloader in `cli.py` (Click command + a `download_dataset` function the loader can wire to `lifecycle.download`). Not every source needs one — see [Downloader or documented acquisition](docs_site/guide/data-sources.md#downloader-or-documented-acquisition) for when to build a downloader versus document manual acquisition.
5. Implement the drift probe — a zero-arg function returning the fingerprint dict described in `hvantk/skills/_conventions/SKILL.md` § 12. Commit the expected fingerprint at `tests/drift_fingerprint.json`.
6. Add the round-trip test under `tests/`, tag it `@pytest.mark.round_trip("<provider>:<dataset>")` (for a round trip parametrized over several datasets, put the mark on each case's `pytest.param(..., marks=pytest.mark.round_trip(...))` instead of stacking every id on the function — otherwise every case carries every dataset's id), and mark it `@pytest.mark.hail` too when the dataset's `backend` is `hail`. Then snapshot the schema and a small set of rows. `hvantk/tests/test_plugin_contract_artifacts.py` replays the manifest's `tests.command` in a subprocess and asserts it selects exactly that one marked item, clean of a skip/xfail mark, and hail-marked when required; run it locally with `-m round_trip_contract` (it is slow, so pytest.ini deselects it by default — CI's `plugin-contract` job runs it on every PR). `--require-round-trips` additionally fails a run if the round trip ends skipped or xfailed for real, which is the one failure mode (e.g. a `pytest.skip()` call inside the test body) collection cannot see (#432).
7. Write `SKILL.md` following the nine-section template in `_conventions/SKILL.md` § 2.
8. Verify: `hvantk plugins validate` should accept the manifest; `hvantk plugins describe <provider>` should list the new dataset without a `FAILED TO LOAD` mark, and `hvantk plugins errors` should report `(no load errors)` (it exits 1 when anything failed to load, while `describe` still exits 0 when only one of the provider's datasets failed); the round-trip test must pass; `hvantk drift <provider:dataset>` must match the committed fingerprint.

Read `hvantk/skills/_conventions/SKILL.md` for the full contract (registry keys, drift-probe shape, lifecycle stages, validation paths, hard guardrails).

## Skill maintenance (resource-centric skills)

`hvantk/skills/` contains agent-readable methodology for maintaining stable resources. Each per-resource `SKILL.md` is a design contract (source identity, output schema, builder + CLI + registry wiring, validation contract) that the round-trip test pins down. Shared conventions live in `hvantk/skills/_conventions/SKILL.md` — read that first before any per-resource skill.

### When you change builder code

If a PR changes a function or path named in any `SKILL.md` — including the `_conventions/` skill — update the affected skill in the same PR. Reviewer will reject PRs where skill drift is obvious.

### Skill status

Each `SKILL.md` carries a `status` frontmatter field:

- `provisional` — work in progress, not yet validated against the round-trip + update tests.
- `stable` — actively maintained, agent can rely on it.
- `deprecated` — superseded by code change; do not use.

Do not invoke an agent against a `deprecated` skill.

### Snapshot regeneration

When a builder change legitimately changes output:

```bash
poetry run pytest hvantk/tests/test_<source>_builder.py --regenerate-snapshots
# For Hail-marked builders only:
poetry run pytest hvantk/tests/test_<source>_builder.py -m hail --regenerate-snapshots
```

Review the diff in `hvantk/tests/snapshots/<source>/` before committing. The regenerated snapshot is the new ground truth.

### Skill-driven changes

Skill-assisted code changes (where an agent generates the builder) are reviewed identically to hand-written code. NEVER auto-merge a skill-generated PR. Always run the full test suite, especially the affected `test_<source>_builder.py`.
