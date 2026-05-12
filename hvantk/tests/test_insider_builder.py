"""Round-trip and update tests for the INSIDER (BED) builder skill."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
)

# Aliased to avoid shadowing the fixture name `regenerate_snapshots` in the test signature
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = "hvantk/tests/testdata/raw/insider/insider_sample.bed"
SNAPSHOT_DIR = Path("hvantk/tests/snapshots/insider")

# `interval` is unique-in-table after .distinct() (per skill §5), so keys are
# inlined here. Picked from observed fixture rows on chr3 (the first track was
# chr11 with zero-length intervals that get skipped). Intervals stored as
# {"start": "<contig>:<pos>", "end": "<contig>:<pos>"} per the new tinterval
# branch in _snapshot_utils.
SAMPLE_KEYS = [
    {"interval": {"start": "chr3:9801671", "end": "chr3:9801673"}},
    {"interval": {"start": "chr3:9801704", "end": "chr3:9801706"}},
    {"interval": {"start": "chr3:9812986", "end": "chr3:9812988"}},
]


@pytest.mark.hail
def test_insider_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build INSIDER BED from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.tables.table_builders import create_interactome_tb

    builder_kwargs = {
        "reference_genome": "GRCh38",
        "overwrite": True,
    }

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=create_interactome_tb,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "insider.ht")
    create_interactome_tb(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    # Idempotency: rebuild with overwrite=True should succeed.
    create_interactome_tb(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "INSIDER schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "INSIDER sample rows drifted from snapshot"
