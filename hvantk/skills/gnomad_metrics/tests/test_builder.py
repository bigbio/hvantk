"""Snapshot round-trip test for the gnomAD constraint metrics builder.

Builds from the committed fixture and asserts the schema and a sample of rows against
committed snapshots, so a change in builder behaviour or upstream field layout shows up
as an explicit diff rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/gnomad_metrics/tests/test_builder.py --regenerate-snapshots
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
    phase_b_snapshot_adapter,
)
# Aliased to avoid shadowing the fixture name `regenerate_snapshots` in the test signature
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = "hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz"
SNAPSHOT_DIR = Path("hvantk/skills/gnomad_metrics/tests/snapshots")

# Real keys, taken from an actual build of the fixture -- collect_sample_rows raises
# KeyError for any key absent from the table, so these cannot be invented.
SAMPLE_KEYS = [
    {"gene_id": "ENSG00000000419"},
    {"gene_id": "ENSG00000019186"},
    {"gene_id": "ENSG00000020256"},
    {"gene_id": "ENSG00000022277"},
    {"gene_id": "ENSG00000025293"},
    {"gene_id": "ENSG00000025772"},
]


@pytest.mark.hail
def test_gnomad_metrics_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build gnomAD constraint metrics from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.gnomad_metrics.builder import build_gnomad_metrics_metrics

    builder = phase_b_snapshot_adapter(build_gnomad_metrics_metrics, "gnomad-metrics:metrics")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "gnomad_metrics.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "gnomAD constraint metrics schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "gnomAD constraint metrics sample rows drifted from snapshot"
