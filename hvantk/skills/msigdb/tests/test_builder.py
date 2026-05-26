"""Round-trip and update tests for the MSigDB builder skill."""

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

_TESTS_DIR = Path(__file__).parent
FIXTURE = str(_TESTS_DIR / "testdata/raw/msigdb/c2.cp-sample.gmt")
SNAPSHOT_DIR = _TESTS_DIR / "snapshots"

# Set names are unique-in-table (per skill §5), so keys are inlined here
# rather than maintained in a separate sample_keys.json (per conventions §9).
# Picks: shortest set, a medium set, the longest set in the fixture, plus a
# typical-prefix set. All four exist in the fixture (see slicer).
SAMPLE_KEYS = [
    {"set_name": "BIOCARTA_ACETAMINOPHEN_PATHWAY"},  # 5 genes (minimum)
    {"set_name": "KEGG_APOPTOSIS"},                  # 87 genes (medium)
    {"set_name": "PID_TCR_PATHWAY"},                 # 64 genes (medium)
    {"set_name": "REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION"},  # 1497 genes (tail)
]


@pytest.mark.hail
def test_msigdb_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build MSigDB from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.msigdb.builder import create_msigdb_tb

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=create_msigdb_tb,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs={"overwrite": True},
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "msigdb.ht")
    create_msigdb_tb(
        input_path=FIXTURE,
        output_path=output_path,
        overwrite=True,
    )
    # Idempotency: rebuild with overwrite=True should succeed.
    create_msigdb_tb(
        input_path=FIXTURE,
        output_path=output_path,
        overwrite=True,
    )
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "MSigDB schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "MSigDB sample rows drifted from snapshot"
