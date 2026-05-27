"""Loader unit + integration tests for hvantk tool manifests."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.tool import loader as tool_loader
from hvantk.core.tool.api import ToolLoadError, ToolSpec
from hvantk.core.tool.loader import ToolRegistry


TOOLS_ROOT = Path(__file__).resolve().parents[1] / "tools"


@pytest.fixture(autouse=True)
def reset_registry():
    tool_loader.reset_registry_for_tests()
    yield
    tool_loader.reset_registry_for_tests()


def test_empty_registry_lists_nothing():
    reg = ToolRegistry()
    assert reg.list_tools() == []
    assert reg.load_errors() == []


def test_loads_three_worked_manifests_from_tools_root():
    reg = ToolRegistry()
    reg.load_from_tools_root()
    names = {t.name for t in reg.list_tools()}
    assert {"plugins", "drift", "reprocess"}.issubset(names)


def test_get_tool_returns_correct_spec_for_drift():
    reg = ToolRegistry()
    reg.load_from_tools_root()
    t = reg.get_tool("drift")
    assert isinstance(t, ToolSpec)
    assert t.name == "drift"
    assert t.domain == "plugins"
    assert t.type == "command"
    assert t.cli_module == "hvantk.tools.plugins.drift_cli"
    assert t.cli_callable == "drift_cmd"
    assert t.requires.hail is True
    assert t.requires.network is True


def test_list_tools_filtered_by_domain_returns_plugin_tools():
    reg = ToolRegistry()
    reg.load_from_tools_root()
    plugin_tools = {t.name for t in reg.list_tools(domain="plugins")}
    assert {"plugins", "drift", "reprocess"}.issubset(plugin_tools)
    # Non-existent domain returns nothing.
    assert reg.list_tools(domain="does-not-exist") == []


def test_broken_manifest_records_load_error(tmp_path: Path):
    bad = tmp_path / "broken.tool.yaml"
    bad.write_text("api_version: 1\nname: BAD_UPPER\n")  # invalid name + missing required fields
    reg = ToolRegistry()
    reg.load_from_tools_root(tmp_path)
    assert reg.list_tools() == []
    errs = reg.load_errors()
    assert len(errs) == 1
    manifest_path, err = errs[0]
    assert "broken.tool.yaml" in manifest_path
    assert isinstance(err, ToolLoadError)


def test_collision_records_load_error(tmp_path: Path):
    """Loading two manifests with the same name records the second as a load error."""
    common_body = """\
api_version: 1
name: dup
domain: infra
type: command
description: dup tool
cli:
  module: hvantk.tools.infra.fake
  function: fake_cmd
purpose:
  short: dup
"""
    (tmp_path / "a.tool.yaml").write_text(common_body)
    (tmp_path / "b.tool.yaml").write_text(common_body)
    reg = ToolRegistry()
    reg.load_from_tools_root(tmp_path)
    assert len(reg.list_tools()) == 1
    errs = reg.load_errors()
    assert len(errs) == 1
    _, err = errs[0]
    assert "collision" in str(err)


def test_get_registry_lazy_singleton_is_stable():
    a = tool_loader.get_registry()
    b = tool_loader.get_registry()
    assert a is b


def test_reset_registry_for_tests_drops_cached_registry():
    a = tool_loader.get_registry()
    tool_loader.reset_registry_for_tests()
    b = tool_loader.get_registry()
    assert a is not b


def test_load_from_missing_tools_root_is_noop(tmp_path: Path):
    reg = ToolRegistry()
    reg.load_from_tools_root(tmp_path / "does-not-exist")
    assert reg.list_tools() == []
    assert reg.load_errors() == []


def test_command_group_manifest_exposes_subcommands():
    reg = ToolRegistry()
    reg.load_from_tools_root()
    t = reg.get_tool("plugins")
    assert t.type == "command_group"
    sub_names = {s.name for s in t.subcommands}
    assert {"list", "describe", "errors", "validate"}.issubset(sub_names)


def test_manifest_path_is_absolute():
    reg = ToolRegistry()
    reg.load_from_tools_root()
    t = reg.get_tool("drift")
    assert Path(t.manifest_path).is_absolute()
    assert t.manifest_path.endswith("drift_cli.tool.yaml")
