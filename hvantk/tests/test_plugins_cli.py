"""Tests for `hvantk plugins` Click commands."""

from __future__ import annotations

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


_VALID_DATASET_YAML = """\
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
"""


def test_validate_rejects_catalog_entry_missing_required_field(tmp_path):
    from click.testing import CliRunner
    from hvantk.tools.plugins.plugins_cli import plugins_group
    (tmp_path / "catalog").mkdir()
    (tmp_path / "catalog" / "datasets.json").write_text(
        '[{"title": "x", "description": "d", "data_source": "ClinGen", '
        '"organism": "Homo sapiens", "files": []}]'  # missing "accession"
    )
    (tmp_path / "plugin.yaml").write_text(
        "api_version: 2\nname: tmp-plug\nversion: 0.1.0\n"
        "catalog: catalog/datasets.json\n"
        "datasets:\n" + _VALID_DATASET_YAML
    )
    res = CliRunner().invoke(plugins_group, ["validate", str(tmp_path / "plugin.yaml")])
    assert res.exit_code != 0
    assert "accession" in res.output


def test_validate_rejects_duplicate_accession_within_catalog(tmp_path):
    from click.testing import CliRunner
    from hvantk.tools.plugins.plugins_cli import plugins_group
    (tmp_path / "catalog").mkdir()
    entry = ('{"accession": "DUP", "title": "x", "description": "d", '
             '"data_source": "ClinGen", "organism": "Homo sapiens", "files": []}')
    (tmp_path / "catalog" / "datasets.json").write_text(f"[{entry}, {entry}]")
    (tmp_path / "plugin.yaml").write_text(
        "api_version: 2\nname: tmp-plug\nversion: 0.1.0\n"
        "catalog: catalog/datasets.json\n"
        "datasets:\n" + _VALID_DATASET_YAML
    )
    res = CliRunner().invoke(plugins_group, ["validate", str(tmp_path / "plugin.yaml")])
    assert res.exit_code != 0
    assert "DUP" in res.output
