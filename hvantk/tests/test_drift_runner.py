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


def test_probe_failure_with_missing_baseline_reports_both(tmp_path: Path):
    """When the probe fails AND no baseline fingerprint is committed, the error
    must surface both the probe failure and the missing-fingerprint path — the
    probe-first ordering must not mask the missing-baseline config error
    (PR #188 review: Copilot)."""
    fp_path = tmp_path / "does-not-exist.json"  # no baseline written

    def boom():
        raise DriftProbeError("network down")

    spec = _make_spec(probe_return=boom, fingerprint_path=fp_path)
    result = _run_with_spec(spec)
    assert result.status == "probe_failed"
    msg = str(result.probe_error)
    assert "network down" in msg  # underlying probe error preserved
    assert "missing" in msg.lower() and str(fp_path) in msg  # baseline flagged


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


def test_placeholder_checksum_baseline_is_probe_failed_not_drifted(tmp_path: Path):
    """A hand-seeded baseline is a missing baseline, not upstream movement.

    Six committed fingerprints carried a placeholder checksum or a Unix-epoch
    `fetched_at`. Neither can equal a live observation, so the comparator called
    them `drifted` on every run and the scheduled bot reopened the same no-op PRs
    nightly. Classifying them `probe_failed` names the actual defect and keeps the
    bot quiet (it never PRs a probe failure).
    """
    from hvantk.core.plugin.api import PLACEHOLDER_CHECKSUM

    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, {
        "probe_version": 1,
        "source_version": None,
        "headers": {"a.tsv": ["col1"]},
        "checksums": {"a.tsv": PLACEHOLDER_CHECKSUM},
        "fetched_at": "2026-05-16T00:00:00+00:00",
    })
    spec = _make_spec(
        probe_return={
            "probe_version": 1,
            "source_version": "Mon, 27 Jul 2026 09:52:02 GMT",
            "headers": {"a.tsv": ["col1"]},
            "checksums": {"a.tsv": "1bf26b3670c0d7ff"},
            "fetched_at": "2026-07-27T09:52:02+00:00",
        },
        fingerprint_path=fp_path,
    )
    result = _run_with_spec(spec)

    assert result.status == "probe_failed"
    assert "--regenerate" in str(result.probe_error)


def test_epoch_fetched_at_baseline_is_probe_failed_not_drifted(tmp_path: Path):
    """The other seeding marker: cptac's two baselines were stamped at the epoch."""
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, {
        "probe_version": 1,
        "source_version": None,
        "headers": {"protein_expression": ["installed_cptac_version"]},
        "checksums": {"protein_expression": ""},
        "fetched_at": "1970-01-01T00:00:00Z",
    })
    spec = _make_spec(
        probe_return={
            "probe_version": 1,
            "source_version": "1.5.13",
            "headers": {"protein_expression": ["installed_cptac_version"]},
            "checksums": {"protein_expression": "abc123"},
            "fetched_at": "2026-07-27T17:00:00+00:00",
        },
        fingerprint_path=fp_path,
    )
    result = _run_with_spec(spec)

    assert result.status == "probe_failed"
    assert "epoch" in str(result.probe_error)


def test_empty_checksums_is_not_treated_as_a_placeholder(tmp_path: Path):
    """ensembl-gene:structure fingerprints HTTP validators, not a body hash.

    Its `checksums` map is legitimately empty. Detection must key on markers that
    cannot occur in genuine probe output, or fixing the false-drifted class would
    just create a false-probe_failed one.
    """
    fp = {
        "probe_version": 2,
        "source_version": "Sun, 18 Aug 2024 22:02:07 GMT",
        "headers": {"g.gtf.gz": {"etag": '"3d2b9d9"', "content_length": "64141785"}},
        "checksums": {},
        "extras": {"release": "113"},
        "fetched_at": "2026-07-27T17:17:16+00:00",
    }
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, fp)
    spec = _make_spec(probe_return=dict(fp), fingerprint_path=fp_path)

    assert _run_with_spec(spec).status == "clean"


def test_empty_checksum_string_baseline_is_probe_failed(tmp_path: Path):
    """An empty string where a digest belongs is seeding, even with a real timestamp.

    cptac's two baselines carried an empty checksum *and* an epoch `fetched_at`, so the
    epoch marker alone covered every real case. A baseline hand-written with a genuine
    timestamp would otherwise have reported drifted forever.
    """
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(
        fp_path,
        {
            "probe_version": 1,
            "headers": {"x": ["a"]},
            "checksums": {"x": ""},
            "fetched_at": "2026-07-27T12:00:00+00:00",
        },
    )
    spec = _make_spec(
        probe_return={
            "probe_version": 1,
            "headers": {"x": ["a"]},
            "checksums": {"x": "a28a4982e5f9"},
            "fetched_at": "2026-07-27T13:00:00+00:00",
        },
        fingerprint_path=fp_path,
    )
    result = _run_with_spec(spec)

    assert result.status == "probe_failed"
    assert "empty" in str(result.probe_error)


def test_empty_checksums_map_is_still_not_a_placeholder(tmp_path: Path):
    """peptideatlas ships ``checksums: {}`` legitimately.

    An empty MAP is not an empty VALUE: widening the marker to cover the former would
    turn a working comparator into a permanent probe_failed, trading one false alarm
    for another.
    """
    fp = {
        "probe_version": 1,
        "headers": {"atlas.zip": {"build_id": "606"}},
        "checksums": {},
        "fetched_at": "2026-07-27T12:00:00+00:00",
    }
    fp_path = tmp_path / "fp.json"
    _write_fingerprint(fp_path, fp)
    spec = _make_spec(probe_return=dict(fp), fingerprint_path=fp_path)

    assert _run_with_spec(spec).status == "clean"
