"""Conformance test for the ensembl-gene plugin (Phase K)."""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.run_builder import run_builder_for_spec


def test_ensembl_gene_genes_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("ensembl-gene:genes")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "ensembl-gene-v1"


@pytest.mark.hail
def test_ensembl_gene_genes_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("ensembl-gene:genes")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "ensembl-gene-v1"

    fixture = Path(
        "hvantk/tests/testdata/raw/ensembl/ensembl_gene_biomart.tsv.bgz"
    )
    assert fixture.exists(), f"Fixture not found: {fixture}"

    out = tmp_path / "genes.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "ensembl-gene"
    assert prov.schema_id == "ensembl-gene-v1"

    from hvantk.core import io as core_io
    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0
