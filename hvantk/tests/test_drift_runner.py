"""Drift runner tests with stubbed probes."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hvantk.core.plugin.drift_runner import DriftResult, run_drift_check
from hvantk.core.plugin.api import DatasetSpec, DriftProbeError, TestPaths


def _make_spec(
    *,
    probe_return,
    fingerprint_path: Path,
) -> DatasetSpec:
    return DatasetSpec(
        name="fake:default",
        domain="genomics",
        backend="hail",
        builder=lambda **kw: None,
        drift_probe=lambda: probe_return() if callable(probe_return) else probe_return,
        skill_path="/abs/SKILL.md",
        test_paths=TestPaths(
            command="pytest -q",
            fixture="/abs/fixture",
            schema_snapshot="/abs/schema.json",
            row_snapshot="/abs/rows.json",
            drift_fingerprint=str(fingerprint_path),
        ),
    )


def _write_fingerprint(path: Path, fp: dict) -> None:
    path.write_text(json.dumps(fp))


def _run_with_spec(spec: DatasetSpec) -> DriftResult:
    """Bypass the registry by calling the internal runner directly."""
    from hvantk.core.plugin.drift_runner import _run_drift_check_with_spec
    return _run_drift_check_with_spec(spec)


def test_clean_when_observed_matches_expected(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    expected = {
        "probe_version": 1,
        "source_version": "v1",
        "headers": {"a.tsv": ["col1"]},
        "checksums": {"a.tsv": "deadbeef"},
        "fetched_at": "2026-01-01T00:00:00Z",
    }
    _write_fingerprint(fp_path, expected)
    observed = dict(expected)
    observed["fetched_at"] = "2099-12-31T00:00:00Z"  # should be ignored
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "clean"
    assert result.diff is None


def test_drifted_when_headers_change(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path,
        {"probe_version": 1, "headers": {"a.tsv": ["col1"]}, "checksums": {"a.tsv": "x"}},
    )
    observed = {
        "probe_version": 1,
        "headers": {"a.tsv": ["col1", "col2"]},
        "checksums": {"a.tsv": "x"},
    }
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "drifted"
    assert result.diff is not None


def test_probe_failed_when_probe_raises(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, {"probe_version": 1, "headers": {}, "checksums": {}})

    def boom():
        raise DriftProbeError("network down")

    spec = _make_spec(probe_return=boom, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "probe_failed"
    assert "network down" in str(result.probe_error)


def test_fetched_at_is_excluded_from_comparison(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path, {"probe_version": 1, "fetched_at": "A", "headers": {}, "checksums": {}}
    )
    observed = {"probe_version": 1, "fetched_at": "B", "headers": {}, "checksums": {}}
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "clean"


def test_run_drift_check_resolves_dataset_from_registry(monkeypatch, tmp_path: Path):
    from hvantk.core.plugin import drift_runner, loader as plugin_loader

    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path, {"probe_version": 1, "headers": {}, "checksums": {}}
    )
    spec = _make_spec(
        probe_return={"probe_version": 1, "headers": {}, "checksums": {}},
        fingerprint_path=fp_path,
    )

    class FakeReg:
        def get_dataset(self, name: str) -> DatasetSpec:
            assert name == "fake:default"
            return spec

    monkeypatch.setattr(plugin_loader, "get_registry", lambda: FakeReg())
    result = drift_runner.run_drift_check("fake:default")
    assert result.status == "clean"


def test_drift_diff_reports_added_keys(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path, {"probe_version": 1, "headers": {}, "checksums": {}}
    )
    observed = {
        "probe_version": 1,
        "headers": {},
        "checksums": {},
        "source_version": "v2",  # new key
    }
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "drifted"
    assert result.diff["added"] == {"source_version": "v2"}
    assert result.diff["removed"] == {}
    assert result.diff["changed"] == {}


def test_drift_diff_reports_removed_keys(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path,
        {
            "probe_version": 1,
            "headers": {},
            "checksums": {},
            "source_version": "v1",
        },
    )
    observed = {"probe_version": 1, "headers": {}, "checksums": {}}
    spec = _make_spec(probe_return=observed, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "drifted"
    assert result.diff["added"] == {}
    assert result.diff["removed"] == {"source_version": "v1"}
    assert result.diff["changed"] == {}


def test_stub_probe_reported_as_stub_not_false_green(tmp_path: Path):
    """A documentation-only stub probe must surface as status='stub' — never a
    false-green 'clean' and never a misleading 'probe_failed' — even with no
    committed baseline fingerprint (issue #177)."""
    from hvantk.core.plugin.api import stub_fingerprint

    # Stub plugins ship no committed baseline, so point at a nonexistent file.
    spec = _make_spec(
        probe_return=lambda: stub_fingerprint("doc-only; no probeable URL"),
        fingerprint_path=tmp_path / "does-not-exist.json",
    )
    result = _run_with_spec(spec)
    assert result.status == "stub"
    assert result.status not in ("clean", "probe_failed")
    assert result.observed["reason"] == "doc-only; no probeable URL"


def test_probe_returning_non_mapping_surfaces_clear_error(tmp_path: Path):
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, {"probe_version": 1, "headers": {}, "checksums": {}})
    spec = _make_spec(
        probe_return=lambda: ["not", "a", "dict"],
        fingerprint_path=fp_path,
    )
    result = _run_with_spec(spec)
    assert result.status == "probe_failed"
    # The error message must mention the actual type so debuggers know what
    # the probe returned, not a cryptic dict() TypeError.
    assert "non-mapping" in str(result.probe_error) or "list" in str(result.probe_error)
