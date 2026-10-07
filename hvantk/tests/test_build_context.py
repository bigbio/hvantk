"""Tests for BuildContext + extended DatasetSpec fields."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

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
        plugin="x",
        dataset="x:y",
        plugin_version="0",
        source_fingerprint="sha256:x",
        builder_commit=None,
    )
    with pytest.raises(TypeError):
        ctx.provenance()  # schema_id is required


def test_build_context_records_and_round_trips_build_parameters():
    from hvantk.core.io._manifest import _from_dict, _to_dict

    ctx = BuildContext(
        plugin="alphagenome",
        dataset="alphagenome:predictions",
        plugin_version="0.2.0",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )
    provenance = ctx.provenance(
        schema_id="alphagenome-v2",
        build_parameters={
            "output_types": ["RNA_SEQ"],
            "ontology_curies": ["UBERON:0006566"],
        },
    )

    assert _from_dict(_to_dict(provenance)) == provenance
    assert provenance.build_parameters["output_types"] == ["RNA_SEQ"]


def test_provenance_manifest_without_build_parameters_remains_readable():
    from hvantk.core.io._manifest import _from_dict, _to_dict

    provenance = BuildContext(
        plugin="legacy",
        dataset="legacy:data",
        plugin_version="0.1.0",
        source_fingerprint="sha256:legacy",
        builder_commit=None,
    ).provenance(schema_id="legacy-v1")
    payload = _to_dict(provenance)
    payload.pop("build_parameters")

    assert _from_dict(payload).build_parameters == {}


def test_build_parameters_are_stamped_as_a_json_copy(tmp_path):
    """A non-JSON value fails at stamping, before save() can write data with no sidecar.

    What is stored is a JSON copy: a later change to the caller's dict does not reach
    it, a tuple is stored as the list the sidecar reloads, and the provenance stays
    hashable while build_parameters still takes part in equality.
    """
    from hvantk.core.io._manifest import read_manifest, write_manifest

    ctx = BuildContext(
        plugin="x",
        dataset="x:y",
        plugin_version="0",
        source_fingerprint="sha256:x",
        builder_commit=None,
    )
    with pytest.raises(TypeError, match="not JSON serializable"):
        ctx.provenance(schema_id="x-y-v1", build_parameters={"input": Path("in.tsv")})
    with pytest.raises(TypeError, match="str keys"):
        ctx.provenance(schema_id="x-y-v1", build_parameters={1: "a", "b": 2})

    params = {"curies": ("UBERON:0006566",)}
    provenance = ctx.provenance(schema_id="x-y-v1", build_parameters=params)
    params["curies"] = ()
    write_manifest(provenance, tmp_path / "rows.parquet")
    assert read_manifest(tmp_path / "rows.parquet") == provenance
    assert hash(provenance) == hash(replace(provenance))
    assert provenance != replace(provenance, build_parameters={})


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
