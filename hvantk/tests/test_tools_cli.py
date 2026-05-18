"""Tests for `hvantk tools` Click commands."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.core import tool_loader
from hvantk.tools.plugins.tools_cli import tools_group


TOOLS_ROOT = Path(__file__).resolve().parents[1] / "tools"


@pytest.fixture(autouse=True)
def reset_registry(monkeypatch):
    tool_loader.reset_registry_for_tests()
    reg = tool_loader.ToolRegistry()
    reg.load_from_tools_root()
    monkeypatch.setattr(tool_loader, "get_registry", lambda: reg)
    yield
    tool_loader.reset_registry_for_tests()


def test_list_command_shows_three_plugin_tools():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["list"])
    assert result.exit_code == 0
    assert "plugins" in result.output
    assert "drift" in result.output
    assert "reprocess" in result.output


def test_list_command_with_domain_filter():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["list", "--domain", "plugins"])
    assert result.exit_code == 0
    assert "drift" in result.output
    assert "plugins" in result.output


def test_describe_drift_outputs_expected_fields():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["describe", "drift"])
    assert result.exit_code == 0
    assert "name:" in result.output
    assert "drift" in result.output
    assert "domain:" in result.output
    assert "plugins" in result.output
    assert "type:" in result.output
    assert "command" in result.output
    assert "cli:" in result.output
    assert "hvantk.tools.plugins.drift_cli" in result.output
    assert "purpose:" in result.output


def test_describe_command_group_shows_subcommands():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["describe", "plugins"])
    assert result.exit_code == 0
    assert "subcommands:" in result.output
    assert "list" in result.output
    assert "describe" in result.output
    assert "validate" in result.output


def test_describe_unknown_tool_errors():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["describe", "nope"])
    assert result.exit_code != 0
    assert "unknown" in result.output.lower()


def test_errors_command_lists_load_errors(monkeypatch, tmp_path):
    bad = tmp_path / "broken.tool.yaml"
    bad.write_text("api_version: 1\nname: BAD UPPER\n")
    reg = tool_loader.ToolRegistry()
    reg.load_from_tools_root(tmp_path)
    monkeypatch.setattr(tool_loader, "get_registry", lambda: reg)

    runner = CliRunner()
    result = runner.invoke(tools_group, ["errors"])
    assert result.exit_code == 0
    assert "broken.tool.yaml" in result.output


def test_errors_command_with_no_errors():
    runner = CliRunner()
    result = runner.invoke(tools_group, ["errors"])
    assert result.exit_code == 0
    assert "no load errors" in result.output.lower()


def test_validate_command_accepts_valid_manifest():
    runner = CliRunner()
    manifest = TOOLS_ROOT / "plugins" / "drift_cli.tool.yaml"
    result = runner.invoke(tools_group, ["validate", str(manifest)])
    assert result.exit_code == 0
    assert "ok" in result.output.lower()


def test_validate_command_rejects_invalid_manifest(tmp_path):
    bad = tmp_path / "bad.tool.yaml"
    bad.write_text("api_version: 1\nname: BAD UPPER\n")
    runner = CliRunner()
    result = runner.invoke(tools_group, ["validate", str(bad)])
    assert result.exit_code != 0
