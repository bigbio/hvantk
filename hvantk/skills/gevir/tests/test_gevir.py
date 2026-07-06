"""Round-trip / snapshot tests for the gevir plugin (Phase K)."""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader
from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
    phase_b_snapshot_adapter,
)

# Aliased to avoid shadowing the fixture name `regenerate_snapshots` in the signature.
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = "hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz"
SNAPSHOT_DIR = Path("hvantk/skills/gevir/tests/snapshots")

# Stable gene_id keys picked from the fixture; the table is keyed by gene_id
# (Ensembl ENSG, str).
SAMPLE_KEYS = [
    {"gene_id": "ENSG00000092607"},  # TBX15
    {"gene_id": "ENSG00000162300"},  # ZFPL1
    {"gene_id": "ENSG00000109062"},  # SLC9A3R1
]


def test_gevir_metrics_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gevir:metrics")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gevir-metrics-v1"


@pytest.mark.hail
def test_gevir_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build GeVIR metrics from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.gevir.builder import build_gevir_metrics

    builder = phase_b_snapshot_adapter(build_gevir_metrics, "gevir:metrics")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "gevir.ht")
    fixture_uri = Path(FIXTURE).resolve().as_uri()
    builder(input_path=fixture_uri, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "GeVIR schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "GeVIR sample rows drifted from snapshot"
