# Plugin System — Phase 0 Scaffolding + HGNC Stub Plugin

> **For agentic workers:** Execute this plan task-by-task — write each failing test, run it red, implement, run it green, commit. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the loader/registry/drift-runner/CLI scaffolding for the plugin system, prove it end-to-end by exposing HGNC as a discovered plugin that still points back to the existing builder in `table_builders.py`.

**Architecture:** A new `hvantk/core/plugin_loader.py` scans `hvantk/skills/*/plugin.yaml` and the `hvantk.providers` entry-point group, validates manifests against `plugin_manifest.schema.json`, constructs `Provider` dataclasses, and populates the existing `TABLE_BUILDERS` / `MATRIX_BUILDERS` dicts alongside the legacy `create_table_adapter` registrations. No legacy registrations are removed in this plan — coexistence is the design until follow-up plans migrate each provider's code.

**Tech Stack:** Python 3.11+, Poetry build, pytest, Click CLI, `jsonschema` for manifest validation, `PyYAML` for parsing.

**Spec reference:** `docs/specs/2026-05-15-data-handlers-plugin-system-design.md` (relocated from `docs/<sp>/specs/`). This plan implements §5–§11 of the spec; per-provider code relocation (§12 Phase 1) is out of scope for this plan and lives in follow-up plans.

---

## File Structure

**Created in this plan:**

| Path | Responsibility |
|---|---|
| `hvantk/core/plugin_api.py` | `DatasetSpec`, `TestPaths`, `Provider` dataclasses; `PluginLoadError`, `DriftProbeError` exceptions |
| `hvantk/core/plugin_manifest.schema.json` | JSON Schema for `plugin.yaml` (api_version 1) |
| `hvantk/core/plugin_loader.py` | `PluginRegistry` class, filesystem + entry-point discovery, lazy module-level `REGISTRY` singleton |
| `hvantk/core/drift_runner.py` | `DriftResult` dataclass, `run_drift_check` function, fingerprint comparator |
| `hvantk/commands/plugins_cli.py` | `hvantk plugins {list,describe,errors,validate}` Click commands |
| `hvantk/commands/drift_cli.py` | `hvantk drift` Click command (`--all`, `--domain`, `--json`, `--regenerate`) |
| `hvantk/skills/hgnc/plugin.yaml` | HGNC stub manifest — `builder.function` points at the existing `create_hgnc_gene_tb` in `table_builders.py` |
| `hvantk/skills/hgnc/__init__.py` | Empty package marker |
| `hvantk/skills/hgnc/drift_probe.py` | HTTP HEAD + first-line fetch against `HGNC_DOWNLOAD_URL`; produces fingerprint |
| `hvantk/skills/hgnc/tests/__init__.py` | Empty test package marker |
| `hvantk/skills/hgnc/tests/drift_fingerprint.json` | Expected fingerprint, seeded by running the probe once |
| `hvantk/tests/test_plugin_api.py` | Unit tests for the dataclasses + error types |
| `hvantk/tests/test_plugin_manifest_schema.py` | JSON schema validation tests (valid + invalid manifests) |
| `hvantk/tests/test_plugin_loader.py` | Loader unit + integration tests using a fake plugin under `hvantk/tests/testdata/raw/plugins/` |
| `hvantk/tests/test_drift_runner.py` | Drift runner unit tests with stubbed probes |
| `hvantk/tests/test_plugins_cli.py` | Click command tests for `hvantk plugins …` |
| `hvantk/tests/test_drift_cli.py` | Click command tests for `hvantk drift …` |
| `hvantk/tests/test_plugin_packaging.py` | Asserts plugin.yaml files are importable as package resources |
| `hvantk/tests/testdata/raw/plugins/fake-plugin/plugin.yaml` | Test fixture: a manifest the loader should accept |
| `hvantk/tests/testdata/raw/plugins/fake-plugin/__init__.py` | Empty |
| `hvantk/tests/testdata/raw/plugins/fake-plugin/builder.py` | Stub builder returning a sentinel |
| `hvantk/tests/testdata/raw/plugins/fake-plugin/drift_probe.py` | Stub probe returning fixed fingerprint |
| `hvantk/tests/testdata/raw/plugins/broken-manifest/plugin.yaml` | Test fixture: invalid YAML for negative-path tests |

**Modified in this plan:**

| Path | Change |
|---|---|
| `hvantk/tables/registry.py` | Add `_apply_plugin_registrations()` call after the existing `create_table_adapter` block; module-level call once per process |
| `pyproject.toml` | Add `[tool.poetry.plugins."hvantk.providers"]` group with `hgnc = "hvantk.skills.hgnc"`; add `jsonschema`, `pyyaml` dependencies if missing; add `[tool.poetry] include` for `hvantk/skills/**/plugin.yaml`; add `exclude` for `hvantk/skills/**/tests/**`; extend pytest `testpaths` |
| `hvantk/hvantk.py` | Attach `plugins` and `drift` subgroups to the top-level `cli` Click group |
| `hvantk/skills/hgnc/SKILL.md` | One-line note added under §6 pointing at the new `plugin.yaml` |

**Untouched in this plan (deferred to follow-up plans):**

- `hvantk/tables/table_builders.py` (HGNC builder stays here — stub points at it)
- `hvantk/commands/hgnc_downloader.py` (downloader stays in place)
- `hvantk/tests/test_hgnc_table_hail.py` / `test_hgnc_downloader.py` (tests stay where they are)
- All other provider code

---

## Task 1: Define plugin API dataclasses

**Files:**
- Create: `hvantk/core/plugin_api.py`
- Test: `hvantk/tests/test_plugin_api.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_plugin_api.py
"""Tests for hvantk.core.plugin_api dataclasses and error types."""

import pytest

from hvantk.core.plugin_api import (
    DatasetSpec,
    PluginLoadError,
    Provider,
    TestPaths,
    DriftProbeError,
)


def _make_test_paths() -> TestPaths:
    return TestPaths(
        command="pytest fake",
        fixture="fake/fixture",
        schema_snapshot="fake/schema.json",
        row_snapshot="fake/rows.json",
        drift_fingerprint="fake/fp.json",
    )


def test_dataset_spec_is_frozen():
    spec = DatasetSpec(
        name="hgnc:lookup",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: {"probe_version": 1},
        skill_path="/abs/SKILL.md",
        test_paths=_make_test_paths(),
    )
    with pytest.raises(Exception):
        spec.name = "other"  # type: ignore[misc]


def test_provider_holds_datasets_tuple():
    spec = DatasetSpec(
        name="hgnc:lookup",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: {"probe_version": 1},
        skill_path="/abs/SKILL.md",
        test_paths=_make_test_paths(),
    )
    provider = Provider(name="hgnc", version="0.1.0", datasets=(spec,))
    assert provider.datasets == (spec,)
    assert isinstance(provider.datasets, tuple)


def test_plugin_load_error_is_exception():
    with pytest.raises(PluginLoadError, match="boom"):
        raise PluginLoadError("boom")


def test_drift_probe_error_is_exception():
    with pytest.raises(DriftProbeError, match="network"):
        raise DriftProbeError("network")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_plugin_api.py -v`
Expected: `ImportError: cannot import name 'DatasetSpec' from 'hvantk.core.plugin_api'`

- [ ] **Step 3: Write minimal implementation**

```python
# hvantk/core/plugin_api.py
"""Stable contracts for hvantk provider plugins.

In-tree and out-of-tree plugin authors import their types from this module.
The Provider dataclass is constructed by the loader (`hvantk.core.plugin_loader`)
from a plugin.yaml manifest plus resolved callables — it is not subclassed or
instantiated by plugin authors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping


class PluginLoadError(Exception):
    """Raised when a plugin fails the protocol or its runtime requirements."""


class DriftProbeError(Exception):
    """Raised by a drift probe on transient failure (network, timeout, parse).

    Distinct from drift itself: the runner classifies a DriftProbeError as
    ``status=probe_failed``, separate from a ``status=drifted`` diff.
    """


@dataclass(frozen=True)
class TestPaths:
    """Validation artifacts for one dataset, as paths resolved from plugin.yaml.

    All paths are absolute (resolved by the loader) so callers do not need to
    know where the plugin folder lives.
    """

    command: str
    fixture: str
    schema_snapshot: str
    row_snapshot: str
    drift_fingerprint: str


@dataclass(frozen=True)
class DatasetSpec:
    """One dataset shipped by a provider plugin.

    `name` is the compound key (e.g. ``"hgnc:lookup"``), composed by the loader
    from ``plugin.name + ":" + dataset.name``. Plugin authors only write the
    bare dataset name in the manifest.
    """

    name: str
    domain: str  # genomics | transcriptomics | proteomics | epigenomics | mapping
    backend: str  # hail | anndata | pandas
    builder: Callable[..., Any]
    drift_probe: Callable[[], Mapping[str, Any]]
    skill_path: str
    test_paths: TestPaths


@dataclass(frozen=True)
class Provider:
    """Registry record for one provider plugin.

    Constructed by the loader; not subclassed or instantiated by plugin
    authors. Equality is structural (frozen dataclass), so two Provider
    objects with the same fields compare equal — useful in tests.
    """

    name: str
    version: str
    datasets: tuple[DatasetSpec, ...]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_plugin_api.py -v`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/core/plugin_api.py hvantk/tests/test_plugin_api.py
git commit -m "feat(plugin-api): add Provider/DatasetSpec/TestPaths dataclasses + error types"
```

---

## Task 2: Define the manifest JSON schema

**Files:**
- Create: `hvantk/core/plugin_manifest.schema.json`
- Test: `hvantk/tests/test_plugin_manifest_schema.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_plugin_manifest_schema.py
"""Tests that plugin.yaml content validates against the JSON schema."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

SCHEMA_PATH = (
    Path(__file__).resolve().parents[1] / "core" / "plugin_manifest.schema.json"
)


@pytest.fixture(scope="module")
def schema() -> dict:
    return json.loads(SCHEMA_PATH.read_text())


def _minimal_manifest() -> dict:
    return {
        "api_version": 1,
        "name": "fake",
        "version": "0.1.0",
        "datasets": [
            {
                "name": "default",
                "domain": "genomics",
                "backend": "hail",
                "builder": {
                    "module": "hvantk.tests.testdata.raw.plugins.fake_plugin.builder",
                    "function": "build",
                },
                "drift_probe": {
                    "module": "hvantk.tests.testdata.raw.plugins.fake_plugin.drift_probe",
                    "function": "fetch_fingerprint",
                },
                "skill": "SKILL.md",
                "tests": {
                    "command": "pytest -q",
                    "fixture": "tests/testdata/raw/fake",
                    "schema_snapshot": "tests/snapshots/schema.json",
                    "row_snapshot": "tests/snapshots/sample_rows.json",
                    "drift_fingerprint": "tests/drift_fingerprint.json",
                },
            }
        ],
    }


def test_minimal_manifest_validates(schema: dict):
    jsonschema.validate(_minimal_manifest(), schema)


def test_unknown_backend_rejected(schema: dict):
    m = _minimal_manifest()
    m["datasets"][0]["backend"] = "spark"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_unknown_domain_rejected(schema: dict):
    m = _minimal_manifest()
    m["datasets"][0]["domain"] = "weather"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_missing_required_field_rejected(schema: dict):
    m = _minimal_manifest()
    del m["datasets"][0]["builder"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_api_version_mismatch_rejected(schema: dict):
    m = _minimal_manifest()
    m["api_version"] = 999
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_invalid_cli_command_name_rejected(schema: dict):
    m = _minimal_manifest()
    m["cli"] = [{"command": "no-prefix", "module": "x", "function": "y"}]
    # Command names must start with the provider name; loader enforces this,
    # but the schema enforces general shape: command must be kebab-case.
    m["cli"][0]["command"] = "BAD UPPER"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_full_manifest_validates(schema: dict):
    m = _minimal_manifest()
    m["status"] = "stable"
    m["maintainers"] = ["alice@example.com"]
    m["description"] = "fake provider"
    m["source"] = {"catalog_ref": "fake", "registry_refs": ["genomics:fake"]}
    m["cli"] = [
        {"command": "fake-download", "module": "x.y.z", "function": "download_cmd"}
    ]
    m["runtime_requirements"] = {"python": ["requests>=2"], "system": []}
    jsonschema.validate(m, schema)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_plugin_manifest_schema.py -v`
Expected: `FileNotFoundError` (schema doesn't exist yet).

- [ ] **Step 3: Write the schema**

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://hvantk.dev/schemas/plugin_manifest.schema.json",
  "title": "hvantk plugin manifest",
  "type": "object",
  "required": ["api_version", "name", "version", "datasets"],
  "additionalProperties": false,
  "properties": {
    "api_version": { "const": 1 },
    "name": {
      "type": "string",
      "pattern": "^[a-z][a-z0-9-]*$"
    },
    "version": {
      "type": "string",
      "pattern": "^[0-9]+\\.[0-9]+\\.[0-9]+(?:-[A-Za-z0-9.-]+)?$"
    },
    "status": {
      "type": "string",
      "enum": ["provisional", "stable", "deprecated"]
    },
    "maintainers": {
      "type": "array",
      "items": { "type": "string" }
    },
    "description": { "type": "string" },
    "source": {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "catalog_ref": { "type": "string" },
        "registry_refs": {
          "type": "array",
          "items": { "type": "string" }
        }
      }
    },
    "datasets": {
      "type": "array",
      "minItems": 1,
      "items": {
        "type": "object",
        "required": ["name", "domain", "backend", "builder", "drift_probe", "skill", "tests"],
        "additionalProperties": false,
        "properties": {
          "name": {
            "type": "string",
            "pattern": "^[a-z][a-z0-9-]*$"
          },
          "domain": {
            "type": "string",
            "enum": ["genomics", "transcriptomics", "proteomics", "epigenomics", "mapping"]
          },
          "backend": {
            "type": "string",
            "enum": ["hail", "anndata", "pandas"]
          },
          "builder": { "$ref": "#/$defs/callable_ref" },
          "drift_probe": { "$ref": "#/$defs/callable_ref" },
          "skill": { "type": "string" },
          "tests": {
            "type": "object",
            "required": [
              "command",
              "fixture",
              "schema_snapshot",
              "row_snapshot",
              "drift_fingerprint"
            ],
            "additionalProperties": false,
            "properties": {
              "command": { "type": "string" },
              "fixture": { "type": "string" },
              "schema_snapshot": { "type": "string" },
              "row_snapshot": { "type": "string" },
              "drift_fingerprint": { "type": "string" }
            }
          }
        }
      }
    },
    "cli": {
      "type": "array",
      "items": {
        "type": "object",
        "required": ["command", "module", "function"],
        "additionalProperties": false,
        "properties": {
          "command": {
            "type": "string",
            "pattern": "^[a-z][a-z0-9-]*$"
          },
          "module": { "type": "string" },
          "function": { "type": "string" }
        }
      }
    },
    "runtime_requirements": {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "python": { "type": "array", "items": { "type": "string" } },
        "system": { "type": "array", "items": { "type": "string" } }
      }
    }
  },
  "$defs": {
    "callable_ref": {
      "type": "object",
      "required": ["module", "function"],
      "additionalProperties": false,
      "properties": {
        "module": { "type": "string" },
        "function": { "type": "string" }
      }
    }
  }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_plugin_manifest_schema.py -v`
Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/core/plugin_manifest.schema.json hvantk/tests/test_plugin_manifest_schema.py
git commit -m "feat(plugin-api): add JSON schema for plugin.yaml (api_version 1)"
```

---

## Task 3: Set up fake plugin fixture (used by loader tests)

**Files:**
- Create: `hvantk/tests/testdata/raw/plugins/__init__.py`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/plugin.yaml`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/__init__.py`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/builder.py`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/drift_probe.py`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/SKILL.md`
- Create: `hvantk/tests/testdata/raw/plugins/fake_plugin/tests/drift_fingerprint.json`
- Create: `hvantk/tests/testdata/raw/plugins/broken-manifest/plugin.yaml`

Folder name uses an underscore (`fake_plugin`, not `fake-plugin`) so the module path is directly importable via `importlib.import_module`. The manifest's `name:` field can still be hyphenated — it's not used as a Python identifier.

- [ ] **Step 1: Create the fake plugin package structure**

```bash
mkdir -p hvantk/tests/testdata/raw/plugins/fake_plugin/tests
mkdir -p hvantk/tests/testdata/raw/plugins/broken-manifest
touch hvantk/tests/testdata/raw/plugins/__init__.py
touch hvantk/tests/testdata/raw/plugins/fake_plugin/__init__.py
```

- [ ] **Step 2: Write the fake plugin manifest**

```yaml
# hvantk/tests/testdata/raw/plugins/fake_plugin/plugin.yaml
api_version: 1
name: fake
version: 0.1.0
status: provisional
description: Fake plugin used by loader tests

datasets:
  - name: default
    domain: genomics
    backend: hail
    builder:
      module: hvantk.tests.testdata.raw.plugins.fake_plugin.builder
      function: build
    drift_probe:
      module: hvantk.tests.testdata.raw.plugins.fake_plugin.drift_probe
      function: fetch_fingerprint
    skill: SKILL.md
    tests:
      command: pytest -q
      fixture: tests/testdata/raw/fake
      schema_snapshot: tests/snapshots/schema.json
      row_snapshot: tests/snapshots/sample_rows.json
      drift_fingerprint: tests/drift_fingerprint.json
```

- [ ] **Step 3: Write the fake builder and probe**

```python
# hvantk/tests/testdata/raw/plugins/fake_plugin/builder.py
"""Stub builder for loader tests."""

SENTINEL_OUTPUT = {"called": True}


def build(input_path: str, output_path: str, **kwargs):
    return {"input_path": input_path, "output_path": output_path, **kwargs, **SENTINEL_OUTPUT}
```

```python
# hvantk/tests/testdata/raw/plugins/fake_plugin/drift_probe.py
"""Stub drift probe returning a deterministic fingerprint."""

FIXED_FINGERPRINT = {
    "probe_version": 1,
    "source_version": "fake-v1",
    "headers": {"a.tsv": ["col1", "col2"]},
    "checksums": {"a.tsv": "deadbeef"},
    "fetched_at": "2026-01-01T00:00:00Z",
}


def fetch_fingerprint() -> dict:
    return dict(FIXED_FINGERPRINT)
```

```markdown
<!-- hvantk/tests/testdata/raw/plugins/fake_plugin/SKILL.md -->
# Fake test plugin
Used by hvantk loader tests. Not a real provider.
```

```json
// hvantk/tests/testdata/raw/plugins/fake_plugin/tests/drift_fingerprint.json
{
  "probe_version": 1,
  "source_version": "fake-v1",
  "headers": {"a.tsv": ["col1", "col2"]},
  "checksums": {"a.tsv": "deadbeef"},
  "fetched_at": "2026-01-01T00:00:00Z"
}
```

- [ ] **Step 4: Write a broken-manifest fixture for negative-path tests**

```yaml
# hvantk/tests/testdata/raw/plugins/broken-manifest/plugin.yaml
api_version: 1
name: BROKEN UPPERCASE
version: not-a-semver
datasets: []
```

- [ ] **Step 5: Commit**

```bash
git add hvantk/tests/testdata/raw/plugins/
git commit -m "test: add fake plugin fixtures for loader test suite"
```

---

## Task 4: Implement the PluginRegistry skeleton

**Files:**
- Create: `hvantk/core/plugin_loader.py`
- Test: `hvantk/tests/test_plugin_loader.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_plugin_loader.py
"""Loader unit + integration tests using fake plugin fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.plugin_api import DatasetSpec, PluginLoadError, Provider
from hvantk.core.plugin_loader import PluginRegistry


FIXTURE_ROOT = Path(__file__).parent / "testdata" / "raw" / "plugins"


def test_empty_registry_lists_nothing():
    reg = PluginRegistry()
    assert reg.list_providers() == []
    assert reg.list_datasets() == []
    assert reg.load_errors() == []


def test_load_fake_plugin_from_filesystem():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    providers = reg.list_providers()
    assert len(providers) == 1
    p = providers[0]
    assert isinstance(p, Provider)
    assert p.name == "fake"
    assert p.version == "0.1.0"
    assert len(p.datasets) == 1
    ds = p.datasets[0]
    assert isinstance(ds, DatasetSpec)
    assert ds.name == "fake:default"
    assert callable(ds.builder)
    assert callable(ds.drift_probe)


def test_dataset_lookup_by_compound_key():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    ds = reg.get_dataset("fake:default")
    assert ds.name == "fake:default"
    with pytest.raises(KeyError):
        reg.get_dataset("does:not:exist")


def test_provider_lookup_by_name():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    p = reg.get_provider("fake")
    assert p.name == "fake"
    with pytest.raises(KeyError):
        reg.get_provider("nope")


def test_broken_manifest_records_error_no_crash():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "broken-manifest")
    assert reg.list_providers() == []
    errors = reg.load_errors()
    assert len(errors) == 1
    plugin_id, err = errors[0]
    assert "broken-manifest" in plugin_id
    assert isinstance(err, PluginLoadError)


def test_collision_raises_hard_error(tmp_path: Path):
    # Copy fake plugin to two locations, both claiming name "fake"
    import shutil
    a = tmp_path / "plugin_a"
    b = tmp_path / "plugin_b"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", a)
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", b)
    reg = PluginRegistry()
    reg.load_from_directory(a)
    with pytest.raises(PluginLoadError, match="collision"):
        reg.load_from_directory(b)


def test_builder_is_invokable_through_dataset_spec():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    ds = reg.get_dataset("fake:default")
    result = ds.builder(input_path="/in", output_path="/out", foo="bar")
    assert result["called"] is True
    assert result["input_path"] == "/in"
    assert result["foo"] == "bar"


def test_drift_probe_is_invokable():
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    ds = reg.get_dataset("fake:default")
    fp = ds.drift_probe()
    assert fp["probe_version"] == 1
    assert fp["headers"]["a.tsv"] == ["col1", "col2"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_plugin_loader.py -v`
Expected: `ImportError: cannot import name 'PluginRegistry'`.

- [ ] **Step 3: Write the loader**

```python
# hvantk/core/plugin_loader.py
"""Plugin discovery and registration for hvantk providers.

Scans `hvantk/skills/*/plugin.yaml` and the `hvantk.providers` entry-point
group, validates each manifest against the schema, resolves the declared
builder/probe callables, and constructs Provider records.

The module-level `REGISTRY` is built lazily on first access.
"""

from __future__ import annotations

import importlib
import json
import logging
from dataclasses import dataclass
from importlib.metadata import entry_points
from pathlib import Path
from typing import Any, Callable, Iterable

import jsonschema
import yaml

from .plugin_api import (
    DatasetSpec,
    PluginLoadError,
    Provider,
    TestPaths,
)

logger = logging.getLogger(__name__)

_SCHEMA_PATH = Path(__file__).parent / "plugin_manifest.schema.json"
_SKILLS_ROOT = Path(__file__).resolve().parent.parent / "skills"
_ENTRY_POINT_GROUP = "hvantk.providers"


def _load_schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


class PluginRegistry:
    """Registry of all loaded provider plugins.

    Construction is empty; populate via `load_from_directory`,
    `load_from_skills_root`, and/or `load_from_entry_points`.
    """

    def __init__(self) -> None:
        self._providers: dict[str, Provider] = {}
        self._datasets: dict[str, DatasetSpec] = {}
        self._load_errors: list[tuple[str, Exception]] = []
        self._schema = _load_schema()

    # --- Public lookup API ---

    def get_provider(self, name: str) -> Provider:
        return self._providers[name]

    def get_dataset(self, name: str) -> DatasetSpec:
        return self._datasets[name]

    def list_providers(self) -> list[Provider]:
        return list(self._providers.values())

    def list_datasets(
        self, *, domain: str | None = None, backend: str | None = None
    ) -> list[DatasetSpec]:
        out = list(self._datasets.values())
        if domain is not None:
            out = [d for d in out if d.domain == domain]
        if backend is not None:
            out = [d for d in out if d.backend == backend]
        return out

    def load_errors(self) -> list[tuple[str, Exception]]:
        return list(self._load_errors)

    # --- Loading entry points ---

    def load_from_skills_root(self, root: Path | None = None) -> None:
        """Scan `hvantk/skills/*/plugin.yaml` (or the given root)."""
        root = root or _SKILLS_ROOT
        if not root.is_dir():
            return
        for child in sorted(root.iterdir()):
            if child.is_dir() and not child.name.startswith("_"):
                if (child / "plugin.yaml").is_file():
                    self.load_from_directory(child)

    def load_from_directory(self, plugin_dir: Path) -> None:
        """Load a single plugin from its directory."""
        plugin_id = str(plugin_dir)
        try:
            manifest = self._read_and_validate_manifest(plugin_dir / "plugin.yaml")
            provider = self._build_provider(manifest, plugin_dir)
            self._register(provider, plugin_id)
        except PluginLoadError as exc:
            # Collision is a hard error: re-raise after recording.
            if "collision" in str(exc):
                raise
            self._load_errors.append((plugin_id, exc))
        except Exception as exc:  # noqa: BLE001 — convert into PluginLoadError for telemetry
            self._load_errors.append((plugin_id, PluginLoadError(str(exc))))

    def load_from_entry_points(self) -> None:
        """Iterate `hvantk.providers` entry points and load each."""
        try:
            eps: Iterable = entry_points(group=_ENTRY_POINT_GROUP)
        except TypeError:
            # Older importlib.metadata: entry_points() returns dict.
            eps = entry_points().get(_ENTRY_POINT_GROUP, [])
        for ep in eps:
            try:
                module = ep.load()
                module_path = Path(module.__file__).resolve().parent
                self.load_from_directory(module_path)
            except PluginLoadError:
                raise
            except Exception as exc:  # noqa: BLE001
                self._load_errors.append(
                    (f"entry-point:{ep.name}", PluginLoadError(str(exc)))
                )

    # --- Internal helpers ---

    def _read_and_validate_manifest(self, manifest_path: Path) -> dict:
        if not manifest_path.is_file():
            raise PluginLoadError(f"missing plugin.yaml at {manifest_path}")
        try:
            content = yaml.safe_load(manifest_path.read_text())
        except yaml.YAMLError as exc:
            raise PluginLoadError(f"invalid YAML: {exc}") from exc
        try:
            jsonschema.validate(content, self._schema)
        except jsonschema.ValidationError as exc:
            raise PluginLoadError(f"schema validation failed: {exc.message}") from exc
        return content

    def _build_provider(self, manifest: dict, plugin_dir: Path) -> Provider:
        datasets: list[DatasetSpec] = []
        for ds_manifest in manifest["datasets"]:
            try:
                ds = self._build_dataset_spec(
                    provider_name=manifest["name"],
                    ds_manifest=ds_manifest,
                    plugin_dir=plugin_dir,
                )
                datasets.append(ds)
            except PluginLoadError as exc:
                # Per-dataset failure: record and continue with the rest.
                self._load_errors.append(
                    (f"{manifest['name']}:{ds_manifest.get('name', '?')}", exc)
                )
        return Provider(
            name=manifest["name"],
            version=manifest["version"],
            datasets=tuple(datasets),
        )

    def _build_dataset_spec(
        self, *, provider_name: str, ds_manifest: dict, plugin_dir: Path
    ) -> DatasetSpec:
        compound = f"{provider_name}:{ds_manifest['name']}"
        builder = self._resolve_callable(
            ds_manifest["builder"]["module"], ds_manifest["builder"]["function"]
        )
        probe = self._resolve_callable(
            ds_manifest["drift_probe"]["module"],
            ds_manifest["drift_probe"]["function"],
        )
        tests = ds_manifest["tests"]
        test_paths = TestPaths(
            command=tests["command"],
            fixture=str((plugin_dir / tests["fixture"]).resolve()),
            schema_snapshot=str((plugin_dir / tests["schema_snapshot"]).resolve()),
            row_snapshot=str((plugin_dir / tests["row_snapshot"]).resolve()),
            drift_fingerprint=str((plugin_dir / tests["drift_fingerprint"]).resolve()),
        )
        skill_path = str((plugin_dir / ds_manifest["skill"]).resolve())
        return DatasetSpec(
            name=compound,
            domain=ds_manifest["domain"],
            backend=ds_manifest["backend"],
            builder=builder,
            drift_probe=probe,
            skill_path=skill_path,
            test_paths=test_paths,
        )

    def _resolve_callable(self, module_path: str, func_name: str) -> Callable[..., Any]:
        try:
            module = importlib.import_module(module_path)
        except ImportError as exc:
            raise PluginLoadError(
                f"cannot import module '{module_path}': {exc}"
            ) from exc
        try:
            func = getattr(module, func_name)
        except AttributeError as exc:
            raise PluginLoadError(
                f"module '{module_path}' has no function '{func_name}'"
            ) from exc
        if not callable(func):
            raise PluginLoadError(
                f"'{module_path}.{func_name}' is not callable"
            )
        return func

    def _register(self, provider: Provider, plugin_id: str) -> None:
        if provider.name in self._providers:
            raise PluginLoadError(
                f"provider name collision: '{provider.name}' already registered "
                f"(new attempt from {plugin_id})"
            )
        self._providers[provider.name] = provider
        for ds in provider.datasets:
            if ds.name in self._datasets:
                raise PluginLoadError(
                    f"dataset name collision: '{ds.name}'"
                )
            self._datasets[ds.name] = ds


# --- Module-level singleton (lazy) ---

_REGISTRY: PluginRegistry | None = None


def get_registry() -> PluginRegistry:
    """Return the module-level registry, building it on first access."""
    global _REGISTRY
    if _REGISTRY is None:
        reg = PluginRegistry()
        reg.load_from_skills_root()
        reg.load_from_entry_points()
        _REGISTRY = reg
    return _REGISTRY


def reset_registry_for_tests() -> None:
    """Test-only: drop the cached registry so the next get_registry() rebuilds it."""
    global _REGISTRY
    _REGISTRY = None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_plugin_loader.py -v`
Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/core/plugin_loader.py hvantk/tests/test_plugin_loader.py
git commit -m "feat(plugin-loader): add PluginRegistry with filesystem + entry-point discovery"
```

---

## Task 5: Wire loader into `tables/registry.py` (alongside legacy registrations)

**Files:**
- Modify: `hvantk/tables/registry.py`
- Test: `hvantk/tests/test_plugin_loader_registry_integration.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_plugin_loader_registry_integration.py
"""Verify that the plugin loader populates TABLE_BUILDERS / MATRIX_BUILDERS
alongside the legacy create_table_adapter registrations."""

from __future__ import annotations

from pathlib import Path

from hvantk.core import plugin_loader


def test_legacy_registrations_still_present():
    from hvantk.tables.registry import TABLE_BUILDERS
    # Sanity: an existing legacy builder is still registered.
    assert "clinvar" in TABLE_BUILDERS


def test_plugin_loader_populates_table_builders(monkeypatch, tmp_path):
    # Point the loader at the fake plugin fixture instead of the real skills dir.
    from hvantk.tables import registry as tables_registry

    fixture = (
        Path(__file__).parent
        / "testdata"
        / "raw"
        / "plugins"
        / "fake_plugin"
    )
    plugin_loader.reset_registry_for_tests()
    # Build a registry seeded only with the fake plugin.
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(fixture)
    # Apply that registry to the legacy dict (function added in this task).
    tables_registry._apply_plugin_registrations(reg)
    assert "fake:default" in tables_registry.TABLE_BUILDERS
    # Legacy entry still intact.
    assert "clinvar" in tables_registry.TABLE_BUILDERS
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_plugin_loader_registry_integration.py -v`
Expected: `AttributeError: module 'hvantk.tables.registry' has no attribute '_apply_plugin_registrations'`.

- [ ] **Step 3: Add the plugin-driven registration to `tables/registry.py`**

Append, at the bottom of `hvantk/tables/registry.py` (after the `MATRIX_BUILDERS` dict definition):

```python
# --- Plugin-driven registrations (added by feat/data-handlers-refactoring) ---


def _apply_plugin_registrations(reg: "PluginRegistry") -> None:
    """Add plugin-discovered builders to the legacy dicts.

    Coexistence: this runs ALONGSIDE the create_table_adapter() block above.
    As each provider migrates to the plugin layout, its create_table_adapter
    line is removed in the same migration commit; this function continues
    to populate the dict from the plugin.
    """
    from hvantk.core.plugin_api import DatasetSpec  # local import to avoid cycle

    def _wrap_builder(spec: "DatasetSpec"):
        if spec.backend == "hail":
            def adapter(input_path, output_path, params=None):
                spec.builder(input_path, output_path, **(params or {}))
            return adapter
        if spec.backend == "anndata":
            def adapter(inputs, output_mt, params=None):
                # AnnData builders take multi-input dicts — pass through unchanged.
                spec.builder(**inputs, output_path=output_mt, **(params or {}))
            return adapter
        # pandas backend: same signature as hail for now.
        def adapter(input_path, output_path, params=None):
            spec.builder(input_path, output_path, **(params or {}))
        return adapter

    for ds in reg.list_datasets(backend="hail"):
        TABLE_BUILDERS[ds.name] = _wrap_builder(ds)
    for ds in reg.list_datasets(backend="anndata"):
        MATRIX_BUILDERS[ds.name] = _wrap_builder(ds)


def _initialize_plugin_registrations() -> None:
    """Wire the module-level PluginRegistry into TABLE_BUILDERS / MATRIX_BUILDERS.

    Called once at module import time. Safe to call repeatedly; subsequent
    calls are no-ops because the registry is a module-level singleton.
    """
    from hvantk.core import plugin_loader

    try:
        reg = plugin_loader.get_registry()
    except Exception as exc:  # noqa: BLE001 — never let plugin failure break hvantk import
        logger.warning("plugin loader failed to initialize: %s", exc)
        return
    _apply_plugin_registrations(reg)


_initialize_plugin_registrations()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_plugin_loader_registry_integration.py -v`
Expected: 2 passed.

Also run the full existing suite to confirm no legacy regression:

Run: `pytest hvantk/tests/test_table_builders_cleanup.py -v`
Expected: existing tests still pass.

- [ ] **Step 5: Commit**

```bash
git add hvantk/tables/registry.py hvantk/tests/test_plugin_loader_registry_integration.py
git commit -m "feat(plugin-loader): populate TABLE_BUILDERS/MATRIX_BUILDERS from plugin registry"
```

---

## Task 6: Implement the drift runner

**Files:**
- Create: `hvantk/core/drift_runner.py`
- Test: `hvantk/tests/test_drift_runner.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_drift_runner.py
"""Drift runner tests with stubbed probes."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from hvantk.core.drift_runner import DriftResult, run_drift_check
from hvantk.core.plugin_api import DatasetSpec, DriftProbeError, TestPaths


def _make_spec(
    *,
    probe_return,
    fingerprint_path: Path,
) -> DatasetSpec:
    return DatasetSpec(
        name="fake:default",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: probe_return() if callable(probe_return) else probe_return,
        skill_path="/abs/SKILL.md",
        test_paths=TestPaths(
            command="pytest -q",
            fixture="/abs/fixture",
            schema_snapshot="/abs/schema.json",
            row_snapshot="/abs/rows.json",
            drift_fingerprint=str(fingerprint_path),
        ),
    )


def _write_fingerprint(path: Path, fp: dict) -> None:
    path.write_text(json.dumps(fp))


def test_clean_when_observed_matches_expected(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    expected = {
        "probe_version": 1,
        "source_version": "v1",
        "headers": {"a.tsv": ["col1"]},
        "checksums": {"a.tsv": "deadbeef"},
        "fetched_at": "2026-01-01T00:00:00Z",
    }
    _write_fingerprint(fp_path, expected)
    observed = dict(expected)
    observed["fetched_at"] = "2099-12-31T00:00:00Z"  # should be ignored
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = run_drift_check_with_spec(spec)
    assert result.status == "clean"
    assert result.diff is None


def test_drifted_when_headers_change(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path,
        {"probe_version": 1, "headers": {"a.tsv": ["col1"]}, "checksums": {"a.tsv": "x"}},
    )
    observed = {
        "probe_version": 1,
        "headers": {"a.tsv": ["col1", "col2"]},
        "checksums": {"a.tsv": "x"},
    }
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = run_drift_check_with_spec(spec)
    assert result.status == "drifted"
    assert result.diff is not None


def test_probe_failed_when_probe_raises(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, {"probe_version": 1, "headers": {}, "checksums": {}})

    def boom():
        raise DriftProbeError("network down")

    spec = _make_spec(probe_return=boom, fingerprint_path=fp_path)
    result = run_drift_check_with_spec(spec)
    assert result.status == "probe_failed"
    assert "network down" in str(result.probe_error)


def test_fetched_at_is_excluded_from_comparison(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path, {"probe_version": 1, "fetched_at": "A", "headers": {}, "checksums": {}}
    )
    observed = {"probe_version": 1, "fetched_at": "B", "headers": {}, "checksums": {}}
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = run_drift_check_with_spec(spec)
    assert result.status == "clean"


# Helper that bypasses the registry for unit testing (real run_drift_check
# looks up the dataset by name; we want to inject the spec directly).
def run_drift_check_with_spec(spec: DatasetSpec) -> DriftResult:
    from hvantk.core.drift_runner import _run_drift_check_with_spec
    return _run_drift_check_with_spec(spec)


def test_run_drift_check_resolves_dataset_from_registry(monkeypatch, tmp_path: Path):
    from hvantk.core import drift_runner, plugin_loader

    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path, {"probe_version": 1, "headers": {}, "checksums": {}}
    )
    spec = _make_spec(
        probe_return={"probe_version": 1, "headers": {}, "checksums": {}},
        fingerprint_path=fp_path,
    )

    class FakeReg:
        def get_dataset(self, name: str) -> DatasetSpec:
            assert name == "fake:default"
            return spec

    monkeypatch.setattr(plugin_loader, "get_registry", lambda: FakeReg())
    result = drift_runner.run_drift_check("fake:default")
    assert result.status == "clean"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_drift_runner.py -v`
Expected: `ImportError: cannot import name 'DriftResult'`.

- [ ] **Step 3: Write the drift runner**

```python
# hvantk/core/drift_runner.py
"""Compare a plugin's live drift-probe fingerprint against its committed expected
fingerprint. Returns a structured DriftResult that the CLI / CI can consume.
"""

from __future__ import annotations

import json
import signal
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .plugin_api import DatasetSpec, DriftProbeError

# Keys excluded from the fingerprint comparison.
_IGNORED_KEYS = frozenset({"fetched_at"})


@dataclass
class DriftResult:
    dataset_name: str
    status: str  # clean | drifted | probe_failed
    observed: dict[str, Any] | None = None
    expected: dict[str, Any] | None = None
    diff: dict[str, Any] | None = None
    probe_error: BaseException | None = None


def run_drift_check(dataset_name: str, *, timeout: int = 60) -> DriftResult:
    """Resolve dataset from the module registry, invoke probe, diff."""
    from . import plugin_loader

    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(dataset_name)
    return _run_drift_check_with_spec(spec, timeout=timeout)


def _run_drift_check_with_spec(
    spec: DatasetSpec, *, timeout: int = 60
) -> DriftResult:
    fp_path = Path(spec.test_paths.drift_fingerprint)
    try:
        expected = json.loads(fp_path.read_text())
    except FileNotFoundError:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            probe_error=DriftProbeError(
                f"missing expected fingerprint at {fp_path}"
            ),
        )

    try:
        observed = _invoke_with_timeout(spec.drift_probe, timeout=timeout)
    except DriftProbeError as exc:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            expected=expected,
            probe_error=exc,
        )
    except Exception as exc:  # noqa: BLE001
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            expected=expected,
            probe_error=DriftProbeError(f"probe raised {type(exc).__name__}: {exc}"),
        )

    diff = _compare_fingerprints(expected, observed)
    if diff is None:
        return DriftResult(
            dataset_name=spec.name,
            status="clean",
            observed=observed,
            expected=expected,
        )
    return DriftResult(
        dataset_name=spec.name,
        status="drifted",
        observed=observed,
        expected=expected,
        diff=diff,
    )


def _invoke_with_timeout(fn, *, timeout: int) -> dict:
    """Run fn() under a signal-based timeout (POSIX). Falls back to no
    timeout on platforms where SIGALRM is unavailable."""
    if not hasattr(signal, "SIGALRM"):
        return dict(fn())

    def _handler(signum, frame):
        raise DriftProbeError(f"probe timed out after {timeout}s")

    prev = signal.signal(signal.SIGALRM, _handler)
    try:
        signal.alarm(timeout)
        result = fn()
        return dict(result)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


def _compare_fingerprints(
    expected: dict[str, Any], observed: dict[str, Any]
) -> dict[str, Any] | None:
    """Return None if equal (ignoring _IGNORED_KEYS), else a structured diff."""
    def strip(d: dict) -> dict:
        return {k: v for k, v in d.items() if k not in _IGNORED_KEYS}

    e = strip(expected)
    o = strip(observed)
    if e == o:
        return None
    added = {k: o[k] for k in o.keys() - e.keys()}
    removed = {k: e[k] for k in e.keys() - o.keys()}
    changed = {
        k: {"expected": e[k], "observed": o[k]}
        for k in e.keys() & o.keys()
        if e[k] != o[k]
    }
    return {"added": added, "removed": removed, "changed": changed}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_drift_runner.py -v`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/core/drift_runner.py hvantk/tests/test_drift_runner.py
git commit -m "feat(drift): add drift runner with fingerprint comparator and timeout"
```

---

## Task 7: Implement `hvantk plugins` Click commands

**Files:**
- Create: `hvantk/commands/plugins_cli.py`
- Test: `hvantk/tests/test_plugins_cli.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_plugins_cli.py
"""Tests for `hvantk plugins` Click commands."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.core import plugin_loader
from hvantk.commands.plugins_cli import plugins_group


FIXTURE_ROOT = Path(__file__).parent / "testdata" / "raw" / "plugins"


@pytest.fixture(autouse=True)
def reset_registry(monkeypatch):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)


def test_list_command_shows_loaded_provider():
    runner = CliRunner()
    result = runner.invoke(plugins_group, ["list"])
    assert result.exit_code == 0
    assert "fake" in result.output
    assert "0.1.0" in result.output


def test_describe_command_dumps_manifest():
    runner = CliRunner()
    result = runner.invoke(plugins_group, ["describe", "fake"])
    assert result.exit_code == 0
    assert "fake:default" in result.output


def test_describe_unknown_provider_errors():
    runner = CliRunner()
    result = runner.invoke(plugins_group, ["describe", "nope"])
    assert result.exit_code != 0
    assert "unknown" in result.output.lower()


def test_errors_command_lists_load_errors(monkeypatch):
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "broken-manifest")
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    runner = CliRunner()
    result = runner.invoke(plugins_group, ["errors"])
    assert result.exit_code == 0
    assert "broken-manifest" in result.output


def test_validate_command_accepts_valid_manifest():
    runner = CliRunner()
    manifest = FIXTURE_ROOT / "fake_plugin" / "plugin.yaml"
    result = runner.invoke(plugins_group, ["validate", str(manifest)])
    assert result.exit_code == 0
    assert "ok" in result.output.lower()


def test_validate_command_rejects_invalid_manifest():
    runner = CliRunner()
    manifest = FIXTURE_ROOT / "broken-manifest" / "plugin.yaml"
    result = runner.invoke(plugins_group, ["validate", str(manifest)])
    assert result.exit_code != 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_plugins_cli.py -v`
Expected: `ImportError`.

- [ ] **Step 3: Write the CLI**

```python
# hvantk/commands/plugins_cli.py
"""`hvantk plugins …` commands for introspecting the plugin registry."""

from __future__ import annotations

import json
from pathlib import Path

import click

from hvantk.core import plugin_loader


@click.group(name="plugins")
def plugins_group():
    """Inspect the hvantk plugin registry."""


@plugins_group.command(name="list")
def list_cmd():
    """List loaded provider plugins."""
    reg = plugin_loader.get_registry()
    providers = reg.list_providers()
    if not providers:
        click.echo("(no plugins loaded)")
        return
    click.echo(f"{'NAME':<24} {'VERSION':<12} {'DATASETS':<8}")
    for p in providers:
        click.echo(f"{p.name:<24} {p.version:<12} {len(p.datasets):<8}")


@plugins_group.command(name="describe")
@click.argument("provider")
def describe_cmd(provider: str):
    """Show full details for one provider."""
    reg = plugin_loader.get_registry()
    try:
        p = reg.get_provider(provider)
    except KeyError:
        raise click.ClickException(f"unknown provider: {provider}")
    click.echo(f"name:    {p.name}")
    click.echo(f"version: {p.version}")
    click.echo("datasets:")
    for ds in p.datasets:
        click.echo(f"  - {ds.name}  ({ds.domain}, {ds.backend})")
        click.echo(f"      skill: {ds.skill_path}")
        click.echo(f"      drift_fingerprint: {ds.test_paths.drift_fingerprint}")


@plugins_group.command(name="errors")
def errors_cmd():
    """List plugins that failed to load and why."""
    reg = plugin_loader.get_registry()
    errs = reg.load_errors()
    if not errs:
        click.echo("(no load errors)")
        return
    for plugin_id, exc in errs:
        click.echo(f"{plugin_id}: {exc}")


@plugins_group.command(name="validate")
@click.argument("manifest_path", type=click.Path(exists=True, dir_okay=False))
def validate_cmd(manifest_path: str):
    """Validate a plugin.yaml file against the schema (offline)."""
    import yaml
    import jsonschema

    schema = json.loads(
        (Path(plugin_loader.__file__).parent / "plugin_manifest.schema.json").read_text()
    )
    content = yaml.safe_load(Path(manifest_path).read_text())
    try:
        jsonschema.validate(content, schema)
    except jsonschema.ValidationError as exc:
        raise click.ClickException(f"validation failed: {exc.message}")
    click.echo(f"ok: {manifest_path}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_plugins_cli.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/commands/plugins_cli.py hvantk/tests/test_plugins_cli.py
git commit -m "feat(cli): add hvantk plugins {list,describe,errors,validate} commands"
```

---

## Task 8: Implement `hvantk drift` Click commands

**Files:**
- Create: `hvantk/commands/drift_cli.py`
- Test: `hvantk/tests/test_drift_cli.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_drift_cli.py
"""Tests for `hvantk drift …` Click commands."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.core import plugin_loader
from hvantk.commands.drift_cli import drift_cmd


FIXTURE_ROOT = Path(__file__).parent / "testdata" / "raw" / "plugins"


@pytest.fixture(autouse=True)
def reset_registry(monkeypatch):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)


def test_drift_clean_exit_zero():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["fake:default"])
    assert result.exit_code == 0
    assert "clean" in result.output.lower()


def test_drift_all_runs_every_dataset():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--all"])
    assert result.exit_code == 0
    assert "fake:default" in result.output


def test_drift_json_output_is_parseable():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--json", "fake:default"])
    assert result.exit_code == 0
    parsed = json.loads(result.output)
    assert parsed["dataset_name"] == "fake:default"
    assert parsed["status"] == "clean"


def test_drift_regenerate_overwrites_fingerprint(tmp_path: Path, monkeypatch):
    # Point the fixture at a tmpdir-copy so we don't mutate the test asset.
    import shutil
    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    old = json.loads(fp_path.read_text())
    # Mutate the expected file so a regenerate visibly changes it.
    fp_path.write_text(json.dumps({"probe_version": 1, "stale": True}))
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--regenerate", "fake:default"])
    assert result.exit_code == 0
    new = json.loads(fp_path.read_text())
    assert "stale" not in new
    assert new["probe_version"] == old["probe_version"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_drift_cli.py -v`
Expected: `ImportError`.

- [ ] **Step 3: Write the CLI**

```python
# hvantk/commands/drift_cli.py
"""`hvantk drift …` command for fingerprint-based drift detection."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import click

from hvantk.core import drift_runner, plugin_loader


EXIT_CLEAN = 0
EXIT_DRIFTED = 1
EXIT_PROBE_FAILED = 2
EXIT_REGISTRY_ERROR = 3


@click.command(name="drift")
@click.argument("dataset", required=False)
@click.option("--all", "all_flag", is_flag=True, help="Run drift check for every dataset")
@click.option("--domain", default=None, help="Filter by domain (with --all)")
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON")
@click.option("--regenerate", is_flag=True, help="Overwrite drift_fingerprint.json with observed probe output")
@click.option("--timeout", default=60, show_default=True, help="Probe timeout in seconds")
def drift_cmd(dataset, all_flag, domain, as_json, regenerate, timeout):
    """Compare a plugin's live drift-probe fingerprint against the expected file."""
    if all_flag and dataset:
        raise click.UsageError("pass either a dataset name or --all, not both")
    if not all_flag and not dataset:
        raise click.UsageError("specify a dataset name or --all")

    reg = plugin_loader.get_registry()

    if regenerate:
        if all_flag:
            raise click.UsageError("--regenerate requires a specific dataset name")
        _regenerate_fingerprint(reg, dataset)
        click.echo(f"regenerated: {dataset}")
        return

    targets = (
        reg.list_datasets(domain=domain) if all_flag else [reg.get_dataset(dataset)]
    )
    results = [drift_runner._run_drift_check_with_spec(spec, timeout=timeout) for spec in targets]

    if as_json:
        click.echo(json.dumps([_serialize(r) for r in results], indent=2, default=str))
    else:
        for r in results:
            click.echo(f"{r.dataset_name}: {r.status}")
            if r.diff:
                click.echo(json.dumps(r.diff, indent=2, default=str))

    exit_codes = {EXIT_CLEAN}
    for r in results:
        if r.status == "drifted":
            exit_codes.add(EXIT_DRIFTED)
        elif r.status == "probe_failed":
            exit_codes.add(EXIT_PROBE_FAILED)
    raise SystemExit(max(exit_codes))


def _regenerate_fingerprint(reg, dataset_name: str) -> None:
    spec = reg.get_dataset(dataset_name)
    observed = spec.drift_probe()
    Path(spec.test_paths.drift_fingerprint).write_text(
        json.dumps(observed, indent=2, default=str)
    )


def _serialize(result: drift_runner.DriftResult) -> dict:
    d = dataclasses.asdict(result)
    if d["probe_error"] is not None:
        d["probe_error"] = str(result.probe_error)
    return d
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_drift_cli.py -v`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add hvantk/commands/drift_cli.py hvantk/tests/test_drift_cli.py
git commit -m "feat(cli): add hvantk drift command (single/all/json/regenerate)"
```

---

## Task 9: Wire `plugins` and `drift` subgroups into the top-level CLI

**Files:**
- Modify: `hvantk/hvantk.py` (top-level Click group `cli`; pyproject.toml's `hvantk = "hvantk.hvantk:main"` script invokes `main()` which calls `cli()`)
- Create: `hvantk/tests/test_top_level_cli.py`

The top-level group is `hvantk.hvantk.cli`. New subcommands are attached via the existing block of `cli.add_command(...)` calls at module scope (around `hvantk/hvantk.py:80-93`).

- [ ] **Step 1: Write the failing test**

```python
# hvantk/tests/test_top_level_cli.py
"""Smoke test: top-level CLI exposes plugins and drift subcommands."""

from click.testing import CliRunner

from hvantk.hvantk import cli


def test_plugins_subcommand_is_attached():
    runner = CliRunner()
    result = runner.invoke(cli, ["plugins", "--help"])
    assert result.exit_code == 0
    assert "list" in result.output


def test_drift_subcommand_is_attached():
    runner = CliRunner()
    result = runner.invoke(cli, ["drift", "--help"])
    assert result.exit_code == 0
    assert "--all" in result.output
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/tests/test_top_level_cli.py -v`
Expected: both tests fail with "No such command 'plugins'" and "No such command 'drift'".

- [ ] **Step 3: Attach the subgroups**

In `hvantk/hvantk.py`, add the imports near the other command imports at the top of the file:

```python
from hvantk.commands.plugins_cli import plugins_group
from hvantk.commands.drift_cli import drift_cmd
```

Then add the two `add_command` calls inside the existing block (around line 80-93):

```python
cli.add_command(plugins_group)
cli.add_command(drift_cmd)
```

Place them adjacent to `cli.add_command(mktable_group)` for grouping symmetry.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest hvantk/tests/test_top_level_cli.py -v`
Expected: 2 passed.

Also smoke-test from the shell:

```bash
poetry run hvantk plugins list
poetry run hvantk drift --help
```

- [ ] **Step 5: Commit**

```bash
git add hvantk/hvantk.py hvantk/tests/test_top_level_cli.py
git commit -m "feat(cli): attach plugins and drift subgroups to top-level hvantk CLI"
```

---

## Task 10: Update `pyproject.toml` — entry points, pytest paths, package data

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Read current state**

```bash
grep -nE '^\[tool\.poetry|^\[tool\.pytest|^\[build-system' pyproject.toml
```

Identify where `[tool.poetry]`, `[tool.poetry.dependencies]`, `[tool.poetry.scripts]`, and `[tool.pytest.ini_options]` live.

- [ ] **Step 2: Add `jsonschema` and `PyYAML` to dependencies (if missing)**

Inspect `[tool.poetry.dependencies]`. If `jsonschema` is absent:

```toml
[tool.poetry.dependencies]
# … existing deps …
jsonschema = "^4.0"
PyYAML = "^6.0"
```

Run `poetry lock --no-update && poetry install` to refresh the lockfile.

- [ ] **Step 3: Add the entry-point group**

```toml
[tool.poetry.plugins."hvantk.providers"]
# In-tree plugins register themselves the same way as external ones.
# Only providers migrated to the plugin layout appear here.
hgnc = "hvantk.skills.hgnc"
```

- [ ] **Step 4: Extend pytest testpaths**

```toml
[tool.pytest.ini_options]
# preserve any existing keys
testpaths = ["hvantk/tests", "hvantk/skills"]
```

- [ ] **Step 5: Include manifests in the wheel, exclude tests**

```toml
[tool.poetry]
# preserve existing keys
include = [
    "hvantk/skills/**/plugin.yaml",
    "hvantk/skills/**/*.py",
    "hvantk/skills/**/SKILL.md",
]
exclude = [
    "hvantk/skills/**/tests/**",
]
```

- [ ] **Step 6: Verify**

```bash
poetry lock --no-update
poetry install
pytest hvantk/tests/test_plugin_api.py hvantk/tests/test_plugin_loader.py -v
poetry build  # build the wheel
unzip -l dist/hvantk-*.whl | grep -E 'skills/hgnc/(plugin\.yaml|tests/)' | head
# Expected: plugin.yaml present, tests/ absent.
```

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml poetry.lock
git commit -m "build: add jsonschema/PyYAML deps, hvantk.providers entry-points, wheel filters"
```

---

## Task 11: Create the HGNC stub plugin

**Files:**
- Create: `hvantk/skills/hgnc/plugin.yaml`
- Create: `hvantk/skills/hgnc/__init__.py`
- Create: `hvantk/skills/hgnc/drift_probe.py`
- Create: `hvantk/skills/hgnc/tests/__init__.py`
- Create: `hvantk/skills/hgnc/tests/drift_fingerprint.json` (seeded by running the probe)
- Test: `hvantk/skills/hgnc/tests/test_drift_probe.py`

- [ ] **Step 1: Write the failing test**

```python
# hvantk/skills/hgnc/tests/test_drift_probe.py
"""HGNC drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE — it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live HGNC endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.skills.hgnc.drift_probe import fetch_fingerprint, HGNC_DOWNLOAD_URL


def test_fetch_fingerprint_shape():
    fake_first_line = "hgnc_id\tsymbol\tname\tstatus\tlocus_type\n"
    with requests_mock.Mocker() as m:
        m.head(HGNC_DOWNLOAD_URL, headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"})
        m.get(HGNC_DOWNLOAD_URL, text=fake_first_line)
        fp = fetch_fingerprint()
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"]["hgnc_complete_set.txt"] == [
        "hgnc_id", "symbol", "name", "status", "locus_type"
    ]
    assert "checksums" in fp
    assert "fetched_at" in fp
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest hvantk/skills/hgnc/tests/test_drift_probe.py -v`
Expected: `ImportError` (probe file doesn't exist yet).

- [ ] **Step 3: Write the probe**

```python
# hvantk/skills/hgnc/drift_probe.py
"""HGNC drift probe: HEAD + first-line fetch against the complete-set TSV."""

from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.core.plugin_api import DriftProbeError
from hvantk.core.constants import HGNC_DOWNLOAD_URL  # existing constant

PROBE_VERSION = 1
_FILENAME = "hgnc_complete_set.txt"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    try:
        head = requests.head(HGNC_DOWNLOAD_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")
        # Stream a small chunk to capture just the header line.
        with requests.get(
            HGNC_DOWNLOAD_URL, timeout=_TIMEOUT_S, stream=True, allow_redirects=True
        ) as resp:
            resp.raise_for_status()
            buf = io.StringIO()
            for chunk in resp.iter_content(chunk_size=4096, decode_unicode=True):
                buf.write(chunk)
                if "\n" in buf.getvalue():
                    break
            first_line = buf.getvalue().splitlines()[0]
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    columns = first_line.split("\t")
    checksum = hashlib.sha256(first_line.encode("utf-8")).hexdigest()
    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {_FILENAME: columns},
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
```

- [ ] **Step 4: Write the stub manifest**

```yaml
# hvantk/skills/hgnc/plugin.yaml
api_version: 1
name: hgnc
version: 0.1.0
status: provisional
description: HGNC complete-set TSV → Hail Table

source:
  catalog_ref: hgnc

datasets:
  - name: lookup
    domain: mapping
    backend: hail
    builder:
      # During Phase 0 the builder still lives in table_builders.py.
      # Migrating it into this folder happens in a follow-up plan.
      module: hvantk.tables.table_builders
      function: create_hgnc_gene_tb
    drift_probe:
      module: hvantk.skills.hgnc.drift_probe
      function: fetch_fingerprint
    skill: SKILL.md
    tests:
      command: pytest hvantk/tests/test_hgnc_table_hail.py -m hail
      fixture: ../../tests/testdata/raw/hgnc
      schema_snapshot: ../../tests/snapshots/hgnc/schema.json
      row_snapshot: ../../tests/snapshots/hgnc/sample_rows.json
      drift_fingerprint: tests/drift_fingerprint.json
```

Note the `../../` in `fixture` and `_snapshot` paths: during Phase 0 the test artifacts have not been moved yet, so the manifest points back at their current location. A follow-up plan removes those `../../` segments when the test files relocate into the plugin folder.

- [ ] **Step 5: Write the empty __init__.py files**

```bash
touch hvantk/skills/hgnc/__init__.py
touch hvantk/skills/hgnc/tests/__init__.py
```

- [ ] **Step 6: Seed the expected drift fingerprint**

Run the probe once against the live HGNC endpoint:

```bash
poetry run python -c "
import json
from hvantk.skills.hgnc.drift_probe import fetch_fingerprint
print(json.dumps(fetch_fingerprint(), indent=2))
" > hvantk/skills/hgnc/tests/drift_fingerprint.json
```

Inspect the result — sanity-check that `headers["hgnc_complete_set.txt"]` starts with `hgnc_id`, `symbol`, `name`, `locus_group`.

- [ ] **Step 7: Add `requests_mock` to dev dependencies (if missing)**

```bash
grep requests_mock pyproject.toml || poetry add --group dev requests-mock
```

- [ ] **Step 8: Run the probe test offline**

Run: `pytest hvantk/skills/hgnc/tests/test_drift_probe.py -v`
Expected: 1 passed.

- [ ] **Step 9: Commit**

```bash
git add hvantk/skills/hgnc/plugin.yaml hvantk/skills/hgnc/__init__.py \
        hvantk/skills/hgnc/drift_probe.py hvantk/skills/hgnc/tests/__init__.py \
        hvantk/skills/hgnc/tests/test_drift_probe.py \
        hvantk/skills/hgnc/tests/drift_fingerprint.json \
        pyproject.toml poetry.lock
git commit -m "feat(hgnc): add stub plugin manifest, drift probe, expected fingerprint"
```

---

## Task 12: End-to-end smoke + packaging verification

**Files:**
- Create: `hvantk/tests/test_plugin_packaging.py`

- [ ] **Step 1: Write the packaging test**

```python
# hvantk/tests/test_plugin_packaging.py
"""Asserts that plugin.yaml files are accessible as package resources after
install (i.e., the Poetry include/exclude rules don't accidentally drop the
manifest from the wheel).
"""

from __future__ import annotations

from importlib.resources import files


def test_hgnc_plugin_yaml_is_packaged():
    res = files("hvantk.skills.hgnc").joinpath("plugin.yaml")
    assert res.is_file()
    content = res.read_text()
    assert "name: hgnc" in content
```

- [ ] **Step 2: Run the packaging test**

Run: `pytest hvantk/tests/test_plugin_packaging.py -v`
Expected: 1 passed.

- [ ] **Step 3: Run the full test suite**

Run: `pytest hvantk/ -q`
Expected: all pre-existing tests still pass, all new tests pass.

- [ ] **Step 4: Smoke-test the CLI end-to-end**

```bash
poetry run hvantk plugins list
# Expected output includes a row for `hgnc 0.1.0` with 1 dataset.

poetry run hvantk plugins describe hgnc
# Expected output shows hgnc:lookup with mapping domain, hail backend, skill_path,
# and the drift_fingerprint path.

poetry run hvantk drift hgnc:lookup
# Expected: "hgnc:lookup: clean" if the live HGNC endpoint hasn't drifted since
# we seeded the fingerprint. Exit code 0.

poetry run hvantk plugins errors
# Expected: "(no load errors)"
```

If any of those four commands report unexpected output, STOP and investigate before committing. Likely cause: a path in `hvantk/skills/hgnc/plugin.yaml` doesn't resolve from the plugin folder.

- [ ] **Step 5: Verify TABLE_BUILDERS coexistence**

```bash
poetry run python -c "
from hvantk.tables.registry import TABLE_BUILDERS
print('legacy clinvar present:', 'clinvar' in TABLE_BUILDERS)
print('plugin hgnc:lookup present:', 'hgnc:lookup' in TABLE_BUILDERS)
"
```

Expected: both report `True`.

- [ ] **Step 6: Commit**

```bash
git add hvantk/tests/test_plugin_packaging.py
git commit -m "test: verify plugin.yaml is packaged and end-to-end discovery works"
```

---

## Task 13: Update HGNC SKILL.md to point at the new plugin.yaml

**Files:**
- Modify: `hvantk/skills/hgnc/SKILL.md`

- [ ] **Step 1: Add a one-line pointer under § 6 (hvantk integration points)**

Open `hvantk/skills/hgnc/SKILL.md`. Under the existing `## 6. hvantk integration points` section, add as a new bullet near the top:

```markdown
- Plugin manifest: `hvantk/skills/hgnc/plugin.yaml` (drives loader registration and `hvantk drift hgnc:lookup`).
```

- [ ] **Step 2: Commit**

```bash
git add hvantk/skills/hgnc/SKILL.md
git commit -m "docs(hgnc): point SKILL.md at the new plugin.yaml manifest"
```

---

## Plan-level verification

Before declaring Phase 0 complete, verify:

- [ ] `pytest hvantk/ -q` exits 0
- [ ] `poetry run hvantk plugins list` shows `hgnc` with 1 dataset
- [ ] `poetry run hvantk plugins describe hgnc` resolves every field without error
- [ ] `poetry run hvantk drift hgnc:lookup` returns `clean` (or `drifted` with a believable diff if HGNC genuinely updated since we seeded the fingerprint)
- [ ] `poetry run hvantk plugins errors` shows no entries
- [ ] `TABLE_BUILDERS["clinvar"]` is still callable (legacy registration intact)
- [ ] `TABLE_BUILDERS["hgnc:lookup"]` is callable (plugin-driven registration works)
- [ ] `poetry build && unzip -l dist/hvantk-*.whl | grep skills/hgnc/plugin.yaml` shows the manifest IS in the wheel
- [ ] `unzip -l dist/hvantk-*.whl | grep skills/hgnc/tests/` shows tests are NOT in the wheel

---

## Out of scope (follow-up plans)

This plan deliberately stops before relocating any real code. Follow-up plans cover:

- **Plan 2**: migrate HGNC code (builder split from `table_builders.py`, constants relocated, test files moved, snapshots seeded with `--regenerate-snapshots`). After Plan 2, HGNC's plugin.yaml references `hvantk.skills.hgnc.builder.build_hgnc_lookup_tb` instead of the legacy path.
- **Plans 3–11**: migrate the remaining 9 providers (clinvar, gtex-eqtl, gwas-catalog, msigdb, insider, expression-atlas, peptideatlas, cptac, ucsc-cellbrowser), one plan per provider, following the same template as Plan 2.
- **Plan 12**: cleanup — rewrite `_conventions/SKILL.md`, update `docs_site/architecture.md` and `CONTRIBUTING.md`, remove now-empty `hvantk/datasets/` modules and old test paths, finalize CHANGELOG.

Each follow-up plan starts on this same branch (`feat/data-handlers-refactoring`) and the whole stack lands as one PR per the spec's migration-scope decision.
