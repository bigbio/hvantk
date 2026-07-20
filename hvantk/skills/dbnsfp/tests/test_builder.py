"""Snapshot round-trip test for the dbNSFP builder.

Builds from the committed fixture and asserts the schema and a sample of rows against
committed snapshots, so a change in builder behaviour or upstream field layout shows up
as an explicit diff rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/dbnsfp/tests/test_builder.py --regenerate-snapshots
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

FIXTURE = "hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz"
SNAPSHOT_DIR = Path("hvantk/skills/dbnsfp/tests/snapshots")

# Real keys, taken from an actual build of the fixture -- collect_sample_rows raises
# KeyError for any key absent from the table, so these cannot be invented.
SAMPLE_KEYS = [
    {"alleles": ["C", "A"], "locus": "chr10:47057"},
    {"alleles": ["C", "G"], "locus": "chr10:47057"},
    {"alleles": ["T", "A"], "locus": "chr10:47058"},
    {"alleles": ["T", "C"], "locus": "chr10:47058"},
    {"alleles": ["T", "G"], "locus": "chr10:47058"},
    {"alleles": ["A", "C"], "locus": "chr10:47059"},
]


@pytest.mark.hail
def test_dbnsfp_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build dbNSFP from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.dbnsfp.builder import build_dbnsfp_variants

    builder = phase_b_snapshot_adapter(build_dbnsfp_variants, "dbnsfp:variants")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "dbnsfp.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "dbNSFP schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "dbNSFP sample rows drifted from snapshot"
