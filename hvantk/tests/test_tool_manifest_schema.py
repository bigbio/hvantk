"""Tests that tool.yaml content validates against the JSON schema."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

SCHEMA_PATH = (
    Path(__file__).resolve().parents[1] / "core" / "tool" / "manifest.schema.json"
)


@pytest.fixture(scope="module")
def schema() -> dict:
    return json.loads(SCHEMA_PATH.read_text())


def _minimal_manifest() -> dict:
    return {
        "api_version": 1,
        "name": "drift",
        "domain": "plugins",
        "type": "command",
        "description": "Run plugin drift probes.",
        "cli": {
            "module": "hvantk.tools.plugins.drift_cli",
            "function": "drift_cmd",
        },
        "purpose": {
            "short": "Detect upstream data-source drift.",
        },
    }


def test_minimal_manifest_validates(schema: dict):
    jsonschema.validate(_minimal_manifest(), schema)


def test_missing_name_rejected(schema: dict):
    m = _minimal_manifest()
    del m["name"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_missing_purpose_rejected(schema: dict):
    m = _minimal_manifest()
    del m["purpose"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_missing_cli_rejected(schema: dict):
    m = _minimal_manifest()
    del m["cli"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_unknown_domain_enum_rejected(schema: dict):
    m = _minimal_manifest()
    m["domain"] = "weather"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_unknown_type_enum_rejected(schema: dict):
    m = _minimal_manifest()
    m["type"] = "plugin"  # only command/command_group allowed
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_api_version_2_rejected(schema: dict):
    m = _minimal_manifest()
    m["api_version"] = 2
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_command_group_with_subcommands_validates(schema: dict):
    m = _minimal_manifest()
    m["name"] = "plugins"
    m["type"] = "command_group"
    m["subcommands"] = [
        {
            "name": "list",
            "description": "List loaded providers",
            "outputs": ["stdout: tabular"],
        },
        {
            "name": "describe",
            "description": "Describe one provider",
            "inputs": [
                {
                    "name": "provider",
                    "type": "str",
                    "required": True,
                    "description": "Provider name",
                }
            ],
            "outputs": ["stdout: provider details"],
        },
    ]
    jsonschema.validate(m, schema)


def test_command_with_subcommands_is_allowed_just_ignored(schema: dict):
    """Schema does not forbid subcommands on a 'command' type; it's just metadata."""
    m = _minimal_manifest()
    assert m["type"] == "command"
    m["subcommands"] = [{"name": "noop", "description": "ignored"}]
    jsonschema.validate(m, schema)


def test_invalid_tool_name_rejected(schema: dict):
    m = _minimal_manifest()
    m["name"] = "BAD UPPER"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_unknown_top_level_property_rejected(schema: dict):
    m = _minimal_manifest()
    m["extra_unexpected_key"] = "nope"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_requires_block_validates(schema: dict):
    m = _minimal_manifest()
    m["requires"] = {"hail": True, "network": True, "extras": ["foo>=1.0"]}
    jsonschema.validate(m, schema)


def test_subcommand_missing_description_rejected(schema: dict):
    m = _minimal_manifest()
    m["type"] = "command_group"
    m["subcommands"] = [{"name": "list"}]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)


def test_subcommand_input_missing_required_field_rejected(schema: dict):
    m = _minimal_manifest()
    m["type"] = "command_group"
    m["subcommands"] = [
        {
            "name": "describe",
            "description": "Describe one provider",
            "inputs": [
                # Missing "required" field
                {"name": "provider", "type": "str", "description": "Provider name"}
            ],
        }
    ]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(m, schema)
