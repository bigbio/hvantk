"""Conformance test for the gnomad-metrics plugin (Phase K)."""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.run_builder import run_builder_for_spec


def test_gnomad_metrics_metrics_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gnomad-metrics:metrics")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gnomad-metrics-v1"


@pytest.mark.hail
def test_gnomad_metrics_metrics_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gnomad-metrics:metrics")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gnomad-metrics-v1"

    fixture = Path(
        "hvantk/tests/testdata/raw/gnomad/"
        "gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz"
    )
    assert fixture.exists(), f"Fixture not found: {fixture}"

    out = tmp_path / "metrics.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "gnomad-metrics"
    assert prov.schema_id == "gnomad-metrics-v1"

    from hvantk.core import io as core_io
    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0
