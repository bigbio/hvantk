"""Tests for `hvantk plugins` Click commands."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.core.plugin import loader as plugin_loader
from hvantk.tools.plugins.plugins_cli import plugins_group


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
