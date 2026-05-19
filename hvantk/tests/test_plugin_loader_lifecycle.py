"""Loader tests for the optional lifecycle.download / lifecycle.parse callables."""

from __future__ import annotations

import shutil
from pathlib import Path

import yaml

from hvantk.core.plugin.loader import PluginRegistry


FIXTURE_ROOT = Path(__file__).parent / "testdata" / "raw" / "plugins"


def test_lifecycle_absent_means_none():
    """The bundled fake_plugin has no lifecycle block; the fields default to None."""
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    ds = reg.get_dataset("fake:default")
    assert ds.download_fn is None
    assert ds.parse_fn is None


def test_lifecycle_fields_resolved_from_manifest(tmp_path: Path):
    """A plugin that declares lifecycle.download/parse has them resolved on the spec."""
    plugin_dir = tmp_path / "lifecycle_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)

    manifest_path = plugin_dir / "plugin.yaml"
    manifest = yaml.safe_load(manifest_path.read_text())
    manifest["api_version"] = 2
    manifest["datasets"][0]["lifecycle"] = {
        "download": {
            "module": "hvantk.tests.testdata.raw.plugins.fake_plugin.builder",
            "function": "build",
        },
        "parse": {
            "module": "hvantk.tests.testdata.raw.plugins.fake_plugin.builder",
            "function": "build",
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest))

    reg = PluginRegistry()
    reg.load_from_directory(plugin_dir)
    assert reg.load_errors() == []
    ds = reg.get_dataset("fake:default")
    assert callable(ds.download_fn)
    assert callable(ds.parse_fn)


def test_lifecycle_api_version_1_still_loads():
    """Backward compat: api_version 1 manifests load without any lifecycle fields."""
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    assert reg.load_errors() == []
    ds = reg.get_dataset("fake:default")
    assert ds.download_fn is None
    assert ds.parse_fn is None
