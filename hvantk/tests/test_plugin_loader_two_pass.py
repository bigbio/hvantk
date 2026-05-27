"""H2 + MH4 fixes: two-pass plugin discovery and manifest-driven downloader registration.

H2 — descriptive manifest pass survives missing optional runtimes.
MH4 — manifest-driven downloader registration via cli: blocks in plugin.yaml.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.api import DatasetManifest


# ---------------------------------------------------------------------------
# H2 — two-pass discovery
# ---------------------------------------------------------------------------


def test_list_manifests_includes_all_plugins():
    """list_manifests() returns every plugin's datasets even if bindings would fail."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    manifests = reg.list_manifests()
    # 20 plugins, some with multiple datasets (ucsc-cellbrowser has 3, cptac has 2).
    assert len(manifests) >= 20
    names = {m.name for m in manifests}
    assert "clinvar:variants" in names
    assert "gevir:metrics" in names
    assert "ucsc-cellbrowser:default" in names


def test_list_manifests_returns_dataset_manifest_instances():
    """All entries returned by list_manifests() are DatasetManifest objects."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    manifests = reg.list_manifests()
    assert all(isinstance(m, DatasetManifest) for m in manifests)


def test_manifest_exposes_descriptive_fields():
    """DatasetManifest carries the right descriptive fields for clinvar:variants."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    m = next(x for x in reg.list_manifests() if x.name == "clinvar:variants")
    assert m.domain == "genomics"
    assert m.backend == "hail"
    assert m.artifact_type_name == "AnnotationTable"
    assert m.schema_id == "clinvar-variants-v1"
    assert m.plugin_name == "clinvar"
    # References are (module, function) string tuples — no callables imported yet.
    assert isinstance(m.builder_ref, tuple)
    assert len(m.builder_ref) == 2
    assert isinstance(m.builder_ref[0], str)
    assert isinstance(m.builder_ref[1], str)


def test_get_dataset_lazily_resolves_callables():
    """First call to get_dataset() resolves callables; second call returns cached object."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec1 = reg.get_dataset("clinvar:variants")
    spec2 = reg.get_dataset("clinvar:variants")
    assert spec1 is spec2  # cache hit — same object
    assert callable(spec1.builder)
    assert callable(spec1.drift_probe)


def test_get_dataset_raises_for_unknown_name():
    """get_dataset() raises KeyError for an unregistered dataset name."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    with pytest.raises(KeyError):
        reg.get_dataset("does:not:exist")


_FAKE_BROKEN_MANIFEST = """
api_version: 2
name: fake-broken
version: 0.1.0
status: provisional
description: A fake plugin with a broken builder import

source:
  catalog_ref: fake-broken

datasets:
  - name: rows
    domain: transcriptomics
    backend: pandas
    builder:
      module: hvantk.skills._nonexistent.fake_broken_builder
      function: build
    drift_probe:
      module: hvantk.skills._nonexistent.fake_broken_drift
      function: fetch_fingerprint
    skill: SKILL.md
    tests:
      command: pytest
      fixture: tests/testdata
      schema_snapshot: tests/snapshots/schema.json
      row_snapshot: tests/snapshots/sample_rows.json
      drift_fingerprint: tests/drift_fingerprint.json
""".strip()


def _make_broken_plugin_dir(tmp_path):
    plugin_dir = tmp_path / "fake-broken"
    plugin_dir.mkdir()
    (plugin_dir / "plugin.yaml").write_text(_FAKE_BROKEN_MANIFEST)
    (plugin_dir / "SKILL.md").write_text("# fake-broken")
    return plugin_dir


def test_missing_optional_runtime_doesnt_drop_manifest(tmp_path):
    """If a builder import would fail, the manifest still appears in list_manifests()."""
    plugin_dir = _make_broken_plugin_dir(tmp_path)

    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)

    # The manifest IS in the descriptive index even though the builder
    # module doesn't exist.
    manifests = reg.list_manifests()
    names = {m.name for m in manifests}
    assert "fake-broken:rows" in names

    # But trying to get the executable spec raises.
    with pytest.raises(Exception):
        reg.get_dataset("fake-broken:rows")


def test_failed_dataset_resolution_is_cached_after_pass2(tmp_path):
    """Pass-2 failures are cached so get_dataset() raises consistently.

    Regression guard for F14: prior to the fix, Pass 2 caught the
    PluginLoadError into _load_errors but did not cache a failure marker,
    so a later get_dataset() call would re-attempt the import. If the
    underlying failure was transient (e.g. a flaky network probe), the
    second attempt could succeed -- leaving the same registry returning
    two different answers in one session. The fix records failures in
    ``_failed_datasets`` and re-raises the cached error.
    """
    plugin_dir = _make_broken_plugin_dir(tmp_path)

    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)

    # Pass 2 cached the failure.
    assert "fake-broken:rows" in reg._failed_datasets
    cached_error = reg._failed_datasets["fake-broken:rows"]

    # Spy on _resolve_spec to prove get_dataset uses the cache (no second
    # resolution attempt). Replacing it with a sentinel that would fail the
    # test if called.
    resolve_calls = []
    original_resolve = reg._resolve_spec

    def _spy(dm, *args, **kwargs):
        resolve_calls.append(dm.name)
        return original_resolve(dm, *args, **kwargs)

    reg._resolve_spec = _spy  # type: ignore[method-assign]

    # First get_dataset call after Pass 2: raises the cached error without
    # re-resolving.
    with pytest.raises(plugin_loader.PluginLoadError) as exc1:
        reg.get_dataset("fake-broken:rows")
    assert exc1.value is cached_error
    assert resolve_calls == [], "get_dataset re-attempted resolution instead of using the cache"

    # Second get_dataset call: still cached, still no resolve.
    with pytest.raises(plugin_loader.PluginLoadError) as exc2:
        reg.get_dataset("fake-broken:rows")
    assert exc2.value is cached_error
    assert resolve_calls == [], "get_dataset re-attempted resolution on second call"


def test_provider_manifests_field_populated():
    """Provider.manifests contains DatasetManifest objects for all declared datasets."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    prov = reg.get_provider("ucsc-cellbrowser")
    # ucsc-cellbrowser declares 3 datasets.
    assert len(prov.manifests) == 3
    assert all(isinstance(m, DatasetManifest) for m in prov.manifests)
    manifest_names = {m.name for m in prov.manifests}
    assert "ucsc-cellbrowser:default" in manifest_names
    assert "ucsc-cellbrowser:adult-ctx" in manifest_names
    assert "ucsc-cellbrowser:dev-ctx" in manifest_names


def test_list_datasets_still_returns_specs():
    """list_datasets() backward-compat: returns DatasetSpec objects (not manifests)."""
    from hvantk.core.plugin.api import DatasetSpec

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    datasets = reg.list_datasets()
    assert len(datasets) >= 20
    assert all(isinstance(d, DatasetSpec) for d in datasets)


# ---------------------------------------------------------------------------
# MH4 — manifest-driven downloader registration
# ---------------------------------------------------------------------------


def test_apply_plugin_downloaders_wires_known_plugins():
    """Known plugins with cli:-download entries get wired into the click group."""
    import click

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    grp = click.Group("download")
    reg.apply_plugin_downloaders(grp)

    expected = {
        "clinvar",
        "hgnc",
        "clingen",
        "gencc",
        "ucsc",
        "expression-atlas",
        "uniprot-ptm",
        "peptideatlas-phospho",
        "cptac-phospho",
    }
    found = set(grp.commands.keys())
    missing = expected - found
    assert not missing, f"manifest-driven wiring missed: {missing}"


def test_apply_plugin_downloaders_idempotent():
    """Calling apply_plugin_downloaders twice is safe (skips duplicates)."""
    import click

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    grp = click.Group("download")
    reg.apply_plugin_downloaders(grp)
    n = len(grp.commands)
    reg.apply_plugin_downloaders(grp)
    assert len(grp.commands) == n


def test_apply_plugin_downloaders_skips_non_download_commands(tmp_path):
    """Commands that don't end in -download are not wired into the download group."""
    import click

    plugin_dir = tmp_path / "fake-toplevel"
    plugin_dir.mkdir()
    (plugin_dir / "plugin.yaml").write_text(
        """
api_version: 1
name: fake-toplevel
version: 0.1.0
status: provisional
description: Plugin with a top-level (non-download) CLI entry

datasets:
  - name: data
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

cli:
  - command: fake-toplevel-cmd
    module: hvantk.tests.testdata.raw.plugins.fake_plugin.builder
    function: build
""".strip()
    )
    (plugin_dir / "SKILL.md").write_text("# fake-toplevel")

    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)

    grp = click.Group("download")
    reg.apply_plugin_downloaders(grp)
    # "fake-toplevel-cmd" does NOT end in "-download", so it must not appear.
    assert "fake-toplevel" not in grp.commands
    assert "fake-toplevel-cmd" not in grp.commands
