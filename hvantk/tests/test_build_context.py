"""Tests for BuildContext + extended DatasetSpec fields."""
from __future__ import annotations

import pytest

from hvantk.core.models.build_context import BuildContext
from hvantk.core.models.provenance import Provenance


def test_build_context_provenance_builds_correctly():
    ctx = BuildContext(
        plugin="clinvar",
        dataset="clinvar:variants",
        plugin_version="0.2.0",
        source_fingerprint="sha256:abc",
        builder_commit="deadbeef",
    )
    prov = ctx.provenance(schema_id="clinvar-variants-v2")
    assert isinstance(prov, Provenance)
    assert prov.plugin == "clinvar"
    assert prov.dataset == "clinvar:variants"
    assert prov.plugin_version == "0.2.0"
    assert prov.source_fingerprint == "sha256:abc"
    assert prov.schema_id == "clinvar-variants-v2"
    assert prov.builder_commit == "deadbeef"
    assert prov.parents == ()


def test_build_context_requires_schema_id_on_provenance():
    ctx = BuildContext(
        plugin="x", dataset="x:y", plugin_version="0",
        source_fingerprint="sha256:x", builder_commit=None,
    )
    with pytest.raises(TypeError):
        ctx.provenance()  # schema_id is required


def test_dataset_spec_accepts_new_fields():
    """DatasetSpec gains artifact_type and schema_id (both Optional, default None)."""
    from hvantk.core.models.annotation_table import AnnotationTable
    from hvantk.core.plugin.api import DatasetSpec, Domain, Backend, TestPaths

    spec = DatasetSpec(
        name="x:y",
        domain="transcriptomics",
        backend="pandas",
        builder=lambda **_: None,
        drift_probe=lambda: {},
        skill_path="hvantk.skills.x",
        test_paths=TestPaths(
            command="pytest fake",
            fixture="fake/fixture",
            schema_snapshot="fake/schema.json",
            row_snapshot="fake/rows.json",
            drift_fingerprint="fake/fp.json",
        ),
        artifact_type=AnnotationTable,
        schema_id="x-y-v1",
    )
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "x-y-v1"


def test_dataset_spec_defaults_artifact_fields_to_none():
    """Backward compat: existing call sites that don't set artifact_type still work."""
    from hvantk.core.plugin.api import DatasetSpec, Domain, Backend, TestPaths

    spec = DatasetSpec(
        name="x:y",
        domain="transcriptomics",
        backend="pandas",
        builder=lambda **_: None,
        drift_probe=lambda: {},
        skill_path="hvantk.skills.x",
        test_paths=TestPaths(
            command="pytest fake",
            fixture="fake/fixture",
            schema_snapshot="fake/schema.json",
            row_snapshot="fake/rows.json",
            drift_fingerprint="fake/fp.json",
        ),
    )
    assert spec.artifact_type is None
    assert spec.schema_id is None
