"""Loader unit + integration tests using fake plugin fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.plugin.api import DatasetSpec, PluginLoadError, Provider
from hvantk.core.plugin.loader import PluginRegistry


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
    import shutil
    a = tmp_path / "plugin_a"
    b = tmp_path / "plugin_b"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", a)
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", b)
    reg = PluginRegistry()
    reg.load_from_directory(a)
    with pytest.raises(PluginLoadError, match="collision"):
        reg.load_from_directory(b)


def test_loading_same_directory_twice_is_idempotent():
    """Loading the same plugin directory twice must not raise a collision.

    The skills-root scan and the entry-point scan can both surface the same
    in-tree plugin once it is listed in pyproject.toml's
    `[project.entry-points."hvantk.providers"]` table. The loader must dedupe
    so this discovery overlap does not crash every CLI invocation.
    """
    reg = PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    providers = reg.list_providers()
    assert [p.name for p in providers] == ["fake"]
    assert reg.load_errors() == []


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


def test_load_from_skills_root_scans_subdirectories(tmp_path: Path):
    """load_from_skills_root finds plugin.yaml in non-underscore subdirs only."""
    import shutil
    # Set up: tmp_path/plugins/{normal_plugin, _skipped_plugin}
    skills_root = tmp_path / "skills"
    normal = skills_root / "normal_plugin"
    skipped = skills_root / "_skipped_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", normal)
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", skipped)
    reg = PluginRegistry()
    reg.load_from_skills_root(skills_root)
    # Both copies declare provider name "fake"; only one (the non-underscored
    # one) should be loaded.
    assert [p.name for p in reg.list_providers()] == ["fake"]


def test_load_from_skills_root_missing_dir_is_noop(tmp_path: Path):
    """If the skills root doesn't exist, load_from_skills_root is a no-op."""
    reg = PluginRegistry()
    reg.load_from_skills_root(tmp_path / "does-not-exist")
    assert reg.list_providers() == []
    assert reg.load_errors() == []


def test_dedicated_collision_exception_class():
    """PluginNameCollision is a PluginLoadError subclass."""
    from hvantk.core.plugin.api import PluginLoadError, PluginNameCollision
    assert issubclass(PluginNameCollision, PluginLoadError)
    import shutil
    # Reuse the collision setup to verify the dedicated class is raised.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        a = Path(tmp) / "a"
        b = Path(tmp) / "b"
        shutil.copytree(FIXTURE_ROOT / "fake_plugin", a)
        shutil.copytree(FIXTURE_ROOT / "fake_plugin", b)
        reg = PluginRegistry()
        reg.load_from_directory(a)
        with pytest.raises(PluginNameCollision):
            reg.load_from_directory(b)
