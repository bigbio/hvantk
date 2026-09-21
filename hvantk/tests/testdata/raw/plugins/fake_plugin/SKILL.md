---
name: hvantk:resource-fake
description: Fake plugin used by hvantk loader and CLI tests. Not a real provider.
status: provisional
backend: hail
domain: genomics
---

# Fake test plugin

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

Used by hvantk loader tests. **Not a real provider** — there is no upstream source
behind it and nothing here should be copied as fact.

It does, however, carry the full nine-section contract from
`hvantk/skills/_conventions/SKILL.md` § 2, because `hvantk plugins validate` checks
every spec a manifest declares and this fixture is validated by
`hvantk/tests/test_plugins_cli.py::test_validate_command_accepts_valid_manifest`.
That makes this file the minimal conforming example: if the contract gains a section,
this fixture fails until it is added here too.

## 1. Status & scope

Test fixture. Exercises manifest loading, dataset registration under the compound key
`fake:default`, and the CLI validation path. No data is built.

## 2. Source identity

None. There is no upstream source.

## 3. Backend choice + reasoning

`hail`, declared in `plugin.yaml` only so the manifest exercises the same code path a
Hail-backed provider does.

## 4. Raw format & gotchas

None — the builder in `builder.py` is a stub and reads no input.

## 5. Output contract

Whatever the stub `build` returns. No schema is asserted.

## 6. hvantk integration points

- Manifest: `plugin.yaml` (`api_version: 1`, deliberately older than the `2` real
  providers use, so the loader's version handling stays covered).
- Builder: `build` in `builder.py`. Drift probe: `fetch_fingerprint` in `drift_probe.py`.
- Consumers: `hvantk/tests/test_plugins_cli.py`, `hvantk/tests/test_plugin_loader.py`.

## 7. Workflow steps

Not applicable — never built. Load it with `reg.load_from_directory(...)`.

## 8. Update playbook

Change this fixture only when a loader or CLI behaviour under test changes.

## 9. Validation contract

- `fixture`: `tests/testdata/raw/fake` — **deliberately absent**.
- `schema_snapshot`: `tests/snapshots/schema.json` — **deliberately absent**.
- `row_snapshot`: `tests/snapshots/sample_rows.json` — **deliberately absent**.
- `drift_fingerprint`: `tests/drift_fingerprint.json` — the one artifact that IS present.
- `command`: `pytest hvantk/tests/test_plugins_cli.py`.

The absent four are the point: they keep `hvantk plugins validate` exercising its
declared-but-absent-artifact warning. The paths are still named here, because § 9's job
is to say what the manifest declares — whether or not it is on disk yet.
