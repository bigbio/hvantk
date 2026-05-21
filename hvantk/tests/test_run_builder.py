"""Tests for the run_builder() orchestrator.

Builds a tiny in-memory pandas-backed dataset through the new contract.
Does NOT touch real plugins — exercises the platform plumbing only.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from hvantk.core.models import AnnotationTable, BuildContext, Provenance
from hvantk.core.plugin.api import DatasetSpec, TestPaths
from hvantk.core.plugin.run_builder import (
    BuilderContractError,
    run_builder_for_spec,
)


def _make_test_paths() -> TestPaths:
    """Minimal valid TestPaths instance (every field required)."""
    return TestPaths(
        command="pytest",
        fixture="x",
        schema_snapshot="x",
        row_snapshot="x",
        drift_fingerprint="x",
    )


def _make_spec(build_fn) -> DatasetSpec:
    return DatasetSpec(
        name="test:rows",
        domain="transcriptomics",
        backend="pandas",
        builder=build_fn,
        drift_probe=lambda: {"fingerprint": "sha256:test-fingerprint"},
        skill_path="hvantk.skills.test",
        test_paths=_make_test_paths(),
        artifact_type=AnnotationTable,
        schema_id="test-rows-v1",
    )


def test_run_builder_returns_provenance(tmp_path):
    def build(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1, 2, 3]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="test-rows-v1"))

    spec = _make_spec(build)
    out = tmp_path / "rows.parquet"
    prov = run_builder_for_spec(
        spec, parsed_input=None, output_path=out, plugin_version="0.1.0"
    )

    assert isinstance(prov, Provenance)
    assert prov.plugin_version == "0.1.0"
    assert prov.source_fingerprint == "sha256:test-fingerprint"
    assert prov.schema_id == "test-rows-v1"
    assert out.exists()
    assert out.with_name("rows.parquet.provenance.json").exists()


def test_run_builder_validates_artifact_type(tmp_path):
    def build_wrong(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return df  # NOT an AnnotationTable

    spec = _make_spec(build_wrong)
    with pytest.raises(BuilderContractError, match="returned"):
        run_builder_for_spec(
            spec, parsed_input=None, output_path=tmp_path / "r.parquet", plugin_version="0.1.0"
        )


def test_run_builder_validates_schema_id_matches_manifest(tmp_path):
    """If the builder stamps a different schema_id than the manifest declares, raise."""
    def build_mismatch(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="wrong-v1"))

    spec = _make_spec(build_mismatch)
    with pytest.raises(BuilderContractError, match="schema_id"):
        run_builder_for_spec(
            spec, parsed_input=None, output_path=tmp_path / "r.parquet", plugin_version="0.1.0"
        )


def test_run_builder_requires_artifact_type_on_spec(tmp_path):
    """If the manifest hasn't been migrated (artifact_type is None), raise clearly."""
    def build_ok(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="test-rows-v1"))

    spec = DatasetSpec(
        name="legacy:rows",
        domain="transcriptomics",
        backend="pandas",
        builder=build_ok,
        drift_probe=lambda: {"fingerprint": "x"},
        skill_path="hvantk.skills.legacy",
        test_paths=_make_test_paths(),
        # artifact_type intentionally omitted (None)
    )
    with pytest.raises(BuilderContractError, match="no artifact_type"):
        run_builder_for_spec(
            spec, parsed_input=None, output_path=tmp_path / "r.parquet", plugin_version="0.1.0"
        )


def test_run_builder_coerces_real_probe_dict(tmp_path):
    """Real probes return source_version, headers, checksums — not 'fingerprint'."""
    def real_probe():
        return {
            "probe_version": "1",
            "source_version": "Wed, 01 Jan 2026 00:00:00 GMT",
            "headers": {"file": {"content_length": "12345"}},
            "checksums": {},
            "fetched_at": "2026-01-01T00:00:00+00:00",
        }

    def build(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="test-rows-v1"))

    spec = _make_spec(build)
    object.__setattr__(spec, "drift_probe", real_probe)
    out = tmp_path / "rows.parquet"
    prov = run_builder_for_spec(spec, parsed_input=None, output_path=out, plugin_version="0.1.0")

    assert prov.source_fingerprint.startswith("sha256:")
    assert prov.source_fingerprint != "sha256:"
    assert "<no-fingerprint>" not in prov.source_fingerprint


def test_run_builder_honors_explicit_fingerprint_key(tmp_path):
    """If the probe returns a 'fingerprint' key, use it verbatim."""
    def explicit_probe():
        return {"fingerprint": "sha256:explicit-value", "fetched_at": "ignored"}

    def build(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="test-rows-v1"))

    spec = _make_spec(build)
    object.__setattr__(spec, "drift_probe", explicit_probe)
    out = tmp_path / "rows.parquet"
    prov = run_builder_for_spec(spec, parsed_input=None, output_path=out, plugin_version="0.1.0")

    assert prov.source_fingerprint == "sha256:explicit-value"


def test_run_builder_rejects_non_dict_probe(tmp_path):
    """Probe returning something other than a dict raises BuilderContractError."""
    def bad_probe():
        return "not a dict"

    def build(parsed, ctx, **params):
        df = pd.DataFrame({"x": [1]})
        return AnnotationTable.from_pandas(df, provenance=ctx.provenance(schema_id="test-rows-v1"))

    spec = _make_spec(build)
    object.__setattr__(spec, "drift_probe", bad_probe)
    with pytest.raises(BuilderContractError, match="drift probe"):
        run_builder_for_spec(spec, parsed_input=None, output_path=tmp_path / "r.parquet", plugin_version="0.1.0")
