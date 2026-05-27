"""Tests for hvantk.core.plugin_api dataclasses and error types."""

import dataclasses

import pytest

from hvantk.core.plugin.api import (
    DatasetSpec,
    DriftProbeError,
    PluginLoadError,
    Provider,
    TestPaths,
)


def _make_test_paths() -> TestPaths:
    return TestPaths(
        command="pytest fake",
        fixture="fake/fixture",
        schema_snapshot="fake/schema.json",
        row_snapshot="fake/rows.json",
        drift_fingerprint="fake/fp.json",
    )


def test_dataset_spec_is_frozen():
    spec = DatasetSpec(
        name="hgnc:lookup",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: {"probe_version": 1},
        skill_path="/abs/SKILL.md",
        test_paths=_make_test_paths(),
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        spec.name = "other"  # type: ignore[misc]


def test_provider_holds_datasets_tuple():
    spec = DatasetSpec(
        name="hgnc:lookup",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: {"probe_version": 1},
        skill_path="/abs/SKILL.md",
        test_paths=_make_test_paths(),
    )
    provider = Provider(name="hgnc", version="0.1.0", datasets=(spec,))
    assert provider.datasets == (spec,)
    assert isinstance(provider.datasets, tuple)


def test_provider_structural_equality():
    spec = DatasetSpec(
        name="hgnc:lookup",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: {"probe_version": 1},
        skill_path="/abs/SKILL.md",
        test_paths=_make_test_paths(),
    )
    p1 = Provider(name="hgnc", version="0.1.0", datasets=(spec,))
    p2 = Provider(name="hgnc", version="0.1.0", datasets=(spec,))
    assert p1 == p2
    assert hash(p1) == hash(p2)


def test_plugin_load_error_is_exception():
    with pytest.raises(PluginLoadError, match="boom"):
        raise PluginLoadError("boom")


def test_drift_probe_error_is_exception():
    with pytest.raises(DriftProbeError, match="network"):
        raise DriftProbeError("network")
