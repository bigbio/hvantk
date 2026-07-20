"""Snapshot round-trip test for the GenCC builder.

Builds from the committed fixture and asserts the schema and a sample of rows against
committed snapshots, so a change in builder behaviour or upstream field layout shows up
as an explicit diff rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/gencc/tests/test_builder_snapshot.py --regenerate-snapshots
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

FIXTURE = "hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv"
SNAPSHOT_DIR = Path("hvantk/skills/gencc/tests/snapshots")

# Real keys, taken from an actual build of the fixture -- collect_sample_rows raises
# KeyError for any key absent from the table, so these cannot be invented.
SAMPLE_KEYS = [
    {"hgnc_id": "1100", "mondo_id": "0005144", "submitter": "ClinGen"},
    {"hgnc_id": "1100", "mondo_id": "0005144", "submitter": "Genomics England PanelApp"},
    {"hgnc_id": "1100", "mondo_id": "0013683", "submitter": "G2P"},
    {"hgnc_id": "1101", "mondo_id": "0005145", "submitter": "ClinGen"},
    {"hgnc_id": "1101", "mondo_id": "0005145", "submitter": "Genomics England PanelApp"},
    {"hgnc_id": "11998", "mondo_id": "0007903", "submitter": "ClinGen"},
]


@pytest.mark.hail
def test_gencc_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build GenCC from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.gencc.builder import build_gencc_submissions

    builder = phase_b_snapshot_adapter(build_gencc_submissions, "gencc:submissions")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "gencc.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "GenCC schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "GenCC sample rows drifted from snapshot"
