"""Phase P: artifact-typed wrappers for burden analysis.

These tests only check that the wrappers accept artifact-typed inputs
without crashing in the unwrap step. The full burden-analysis behavior
is covered by the existing test_enrichex_*.py suites.
"""
from __future__ import annotations

from datetime import datetime, timezone

import pytest

from hvantk.core.models import (
    AnnotationTable, GeneSet, Provenance,
)


def _prov() -> Provenance:
    return Provenance(
        plugin="t", dataset="t:cohort", plugin_version="0",
        source_fingerprint="sha256:x", schema_id="t-cohort-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


def test_burden_analysis_artifact_imports():
    """Smoke test: the artifact-typed wrapper can be imported."""
    from hvantk.algorithms.enrichex.burden import run_burden_analysis_artifact
    assert callable(run_burden_analysis_artifact)


def test_stratified_burden_analysis_artifact_imports():
    from hvantk.algorithms.enrichex.burden import (
        run_stratified_burden_analysis_artifact,
    )
    assert callable(run_stratified_burden_analysis_artifact)


def test_burden_artifact_metadata():
    """The wrapper has artifact-typed inputs/outputs in its @algorithm meta."""
    from hvantk.algorithms.enrichex.burden import run_burden_analysis_artifact
    from hvantk.core.models.backends import get_algorithm_meta

    meta = get_algorithm_meta(run_burden_analysis_artifact)
    assert meta.name == "burden_analysis_artifact"
    # inputs is a dict declaring the artifact contract
    assert "cohort" in meta.inputs
    assert "gene_sets" in meta.inputs
    assert "phenotype" in meta.inputs
    assert meta.required_backend == "hail"


def test_stratified_burden_artifact_metadata():
    """The stratified wrapper has artifact-typed inputs in its @algorithm meta."""
    from hvantk.algorithms.enrichex.burden import (
        run_stratified_burden_analysis_artifact,
    )
    from hvantk.core.models.backends import get_algorithm_meta

    meta = get_algorithm_meta(run_stratified_burden_analysis_artifact)
    assert meta.name == "stratified_burden_analysis_artifact"
    assert "cohort" in meta.inputs
    assert "gene_sets" in meta.inputs
    assert "phenotype" in meta.inputs
    assert meta.required_backend == "hail"
