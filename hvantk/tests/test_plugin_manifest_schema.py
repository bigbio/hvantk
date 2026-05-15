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
    m["cli"] = [{"command": "BAD UPPER", "module": "x", "function": "y"}]
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


def test_invalid_semver_rejected(schema: dict):
    m = _minimal_manifest()
    m["version"] = "1.0"  # missing patch
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)

    m["version"] = "v1.0.0"  # leading v not allowed
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_empty_datasets_array_rejected(schema: dict):
    m = _minimal_manifest()
    m["datasets"] = []
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_unknown_top_level_property_rejected(schema: dict):
    m = _minimal_manifest()
    m["extra_unexpected_key"] = "should not be allowed"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)
