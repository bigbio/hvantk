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


def test_a_dataset_that_fails_to_bind_is_logged_not_only_collected(tmp_path, caplog):
    """Anything dropped from the registry must say so (#351).

    `_load_errors` alone is not enough: `plugins list` prints the survivors with no
    "N failed" line, and `drift-health.yml` greps the report for `probe_failed` -- so a
    dataset that never registered emits no row at all and reads as healthy. A dataset can
    go from monitored to unmonitored with every signal still green.

    The dataset-binding path specifically, because it is the likelier failure and was
    still silent after the first fix: only the two directory-load paths had been routed
    through the logging helper.
    """
    plugin = tmp_path / "broken"
    plugin.mkdir()
    (plugin / "plugin.yaml").write_text(
        "api_version: 2\n"
        "name: brokenprov\n"
        "version: 0.1.0\n"
        "datasets:\n"
        "  - name: thing\n"
        "    domain: genomics\n"
        "    backend: hail\n"
        "    builder: {module: hvantk.nope.missing, function: build_thing}\n"
        "    drift_probe: {module: hvantk.nope.missing, function: fetch_fingerprint}\n"
        "    skill: SKILL.md\n"
        "    tests: {command: pytest, fixture: f, schema_snapshot: s,\n"
        "            row_snapshot: r, drift_fingerprint: d}\n"
    )
    (plugin / "SKILL.md").write_text("---\nname: x\ndescription: y\n---\n# x\n")

    reg = PluginRegistry()
    with caplog.at_level("WARNING", logger="hvantk.core.plugin.loader"):
        reg.load_from_directory(plugin)

    assert "brokenprov:thing" in [unit for unit, _ in reg.load_errors()]
    logged = [r.getMessage() for r in caplog.records]
    assert any("brokenprov:thing" in m for m in logged), (
        f"the failure was collected but never logged: {logged}"
    )


def test_every_load_error_is_recorded_through_the_one_helper():
    """Enforce the invariant `_record_load_error`'s docstring asserts.

    That docstring says "every path that drops something from the registry must come
    through here", and nothing checked it. Reverting individual call sites to a bare
    `self._load_errors.append(...)` -- which collects the error but emits no warning --
    passed the whole suite for four of the five sites. Per-site tests cannot scale here:
    each new drop path is a new silent hole, and the failure is invisible by
    construction (the registry just gets smaller).

    A structural check is the right shape: it costs nothing, it covers sites that do not
    exist yet, and it fails loudly with the line number when someone appends directly.
    """
    import ast

    from hvantk.core.plugin import loader as loader_module

    source = Path(loader_module.__file__).read_text()
    tree = ast.parse(source)

    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute) or node.attr != "append":
            continue
        target = node.value
        if (
            isinstance(target, ast.Attribute)
            and target.attr == "_load_errors"
            and isinstance(target.value, ast.Name)
            and target.value.id == "self"
        ):
            offenders.append(node.lineno)

    # Exactly one: the append inside _record_load_error itself.
    recorder = next(
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef) and n.name == "_record_load_error"
    )
    allowed = range(recorder.lineno, recorder.end_lineno + 1)

    outside = [ln for ln in offenders if ln not in allowed]
    assert outside == [], (
        "self._load_errors.append(...) appears outside _record_load_error at line(s) "
        f"{outside} of loader.py. Route it through the helper so the failure is logged "
        "as well as collected -- a silently-collected error is the #351 failure mode."
    )
    assert offenders, (
        "no self._load_errors.append(...) found at all -- this guard is looking for "
        "something that no longer exists and has stopped guarding anything"
    )


def _lazy_manifest(name: str):
    """A manifest injected straight into the registry, bypassing the eager bind pass, so
    the lazy get_dataset path is the FIRST to try importing the callables."""
    from hvantk.core.plugin.api import DatasetManifest, TestPaths

    return DatasetManifest(
        name=name,
        domain="genomics",
        backend="hail",
        skill_path="/abs/SKILL.md",
        test_paths=TestPaths(
            command="pytest",
            fixture="/f",
            schema_snapshot="/s",
            row_snapshot="/r",
            drift_fingerprint="/d",
        ),
        plugin_name=name.split(":")[0],
        builder_ref=("hvantk.nope.missing", "build"),
        drift_probe_ref=("hvantk.nope.missing", "fetch_fingerprint"),
    )


def test_lazy_get_dataset_records_a_binding_failure_once(caplog):
    """`get_dataset` wrote `_failed_datasets` without recording (#364 item 2). Benign only
    because the eager pass pre-populates both dicts, and the docstring promised otherwise."""
    reg = PluginRegistry()
    reg._manifests["lazyprov:thing"] = _lazy_manifest("lazyprov:thing")

    with caplog.at_level("WARNING", logger="hvantk.core.plugin.loader"):
        with pytest.raises(PluginLoadError):
            reg.get_dataset("lazyprov:thing")
        with pytest.raises(PluginLoadError):
            reg.get_dataset("lazyprov:thing")  # cached: must not record twice

    assert [u for u, _ in reg.load_errors()] == ["lazyprov:thing"]
    assert any("lazyprov:thing" in r.getMessage() for r in caplog.records)


def test_lazy_get_dataset_wraps_a_non_plugin_error_like_the_eager_pass(monkeypatch):
    """The eager pass wraps any Exception into PluginLoadError and caches it; the lazy
    path must give the same answer, or the same registry returns two different errors."""
    reg = PluginRegistry()
    reg._manifests["lazyprov:thing"] = _lazy_manifest("lazyprov:thing")

    def explode(*_a, **_k):
        raise RuntimeError("module body blew up")

    monkeypatch.setattr(reg, "_resolve_spec", explode)
    with pytest.raises(PluginLoadError, match="module body blew up"):
        reg.get_dataset("lazyprov:thing")
    assert [u for u, _ in reg.load_errors()] == ["lazyprov:thing"]


def test_list_datasets_skips_a_failed_dataset_without_recording_it_again():
    reg = PluginRegistry()
    reg._manifests["lazyprov:thing"] = _lazy_manifest("lazyprov:thing")
    assert reg.list_datasets() == []
    assert reg.list_datasets() == []


# --- item 2 (#374 review): a provider-level load failure recorded via
# load_from_skills_root must not be recorded AGAIN when load_from_entry_points
# resolves to the SAME directory ------------------------------------------------------
#
# get_registry() runs load_from_skills_root() then load_from_entry_points(); 13 of 23
# in-tree providers are ALSO declared as entry points pointing at their own skills/
# directory. `load_from_directory` only added to `_loaded_dirs` on SUCCESS, and the
# missing-plugin.yaml branch of `load_from_skills_root` never marked the dir at all --
# so a broken provider's directory was never marked "attempted" and the second loader
# tried it again, appending a second (differently-worded) error for the same unit.


def test_broken_provider_directory_records_one_error_across_skills_root_then_entry_point(
    tmp_path: Path,
):
    """Simulates the real get_registry() overlap: `load_from_skills_root` finds the
    broken directory first (via plugin.yaml), then `load_from_entry_points` resolves
    its installed entry point back to the SAME directory and calls
    `load_from_directory` a second time -- exactly what `load_from_entry_points` does
    after `module.__file__`'s parent is computed.
    """
    skills_root = tmp_path / "skills"
    plugin_dir = skills_root / "brokenprov"
    plugin_dir.mkdir(parents=True)
    (plugin_dir / "plugin.yaml").write_text(
        "api_version: 1\nname: BROKEN UPPERCASE\nversion: not-a-semver\ndatasets: []\n"
    )

    reg = PluginRegistry()
    reg.load_from_skills_root(skills_root)
    assert len(reg.load_errors()) == 1

    # The entry-point path: resolve to the identical directory and load it again.
    reg.load_from_directory(plugin_dir)

    errors = reg.load_errors()
    assert len(errors) == 1, errors


def test_missing_manifest_directory_loaded_twice_records_one_error(tmp_path: Path):
    """A directory with no plugin.yaml at all: `load_from_skills_root` records it via
    its own "no plugin.yaml in ..." branch, which must ALSO mark the directory as
    attempted -- otherwise a second discovery pass (the entry-point path) calls
    `load_from_directory` on it directly, which fails again with a differently-worded
    "missing plugin.yaml at ..." message, doubling the row for one broken unit.
    """
    skills_root = tmp_path / "skills"
    plugin_dir = skills_root / "ghost"
    plugin_dir.mkdir(parents=True)

    reg = PluginRegistry()
    reg.load_from_skills_root(skills_root)
    assert len(reg.load_errors()) == 1

    reg.load_from_directory(plugin_dir)

    errors = reg.load_errors()
    assert len(errors) == 1, errors


# --- a downloader whose module will not import is a load error, not a log line -------
#
# `apply_plugin_downloaders` warned and skipped the entry, so `hvantk download <p>`
# answered "No such command" while `plugins errors` stayed clean: a unit dropped from
# the registry with no row, which is the #351 failure mode again.


def test_a_downloader_that_fails_to_import_is_recorded_as_a_load_error(
    tmp_path, caplog
):
    import click

    plugin = tmp_path / "wonky"
    plugin.mkdir()
    (plugin / "plugin.yaml").write_text(
        "api_version: 2\n"
        "name: wonky\n"
        "version: 0.1.0\n"
        "datasets:\n"
        "  - name: thing\n"
        "    domain: genomics\n"
        "    backend: hail\n"
        "    builder: {module: hvantk.tests.testdata.raw.plugins.fake_plugin.builder,\n"
        "              function: build}\n"
        "    drift_probe: {module: hvantk.tests.testdata.raw.plugins.fake_plugin.drift_probe,\n"
        "                  function: fetch_fingerprint}\n"
        "    skill: SKILL.md\n"
        "    tests: {command: pytest, fixture: f, schema_snapshot: s,\n"
        "            row_snapshot: r, drift_fingerprint: d}\n"
        "cli:\n"
        "  - command: wonky-download\n"
        "    module: hvantk.nope.missing\n"
        "    function: download_cmd\n"
    )
    (plugin / "SKILL.md").write_text("---\nname: x\ndescription: y\n---\n# x\n")
    reg = PluginRegistry()
    reg.load_from_directory(plugin)
    assert reg.load_errors() == [], "the dataset binds; only the downloader is broken"

    grp = click.Group("download")
    with caplog.at_level("WARNING", logger="hvantk.core.plugin.loader"):
        reg.apply_plugin_downloaders(grp)

    assert "wonky" not in grp.commands
    assert [u for u, _ in reg.load_errors()] == ["wonky"]
    _, exc = reg.load_errors()[0]
    assert "wonky-download" in str(exc) and "hvantk.nope.missing" in str(exc)
    assert any("wonky" in r.getMessage() for r in caplog.records), "warning kept"
