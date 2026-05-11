"""Round-trip and update tests for the ClinVar builder skill."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
)
# Aliased to avoid shadowing the fixture name `regenerate_snapshots` in the test signature
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = "hvantk/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz"
SNAPSHOT_DIR = Path("hvantk/tests/snapshots/clinvar")


@pytest.mark.hail
def test_clinvar_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build ClinVar from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.tables.table_builders import create_clinvar_tb

    keys = json.loads((SNAPSHOT_DIR / "sample_keys.json").read_text())

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=create_clinvar_tb,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=keys,
            builder_kwargs={"reference_genome": "GRCh38", "overwrite": True},
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "clinvar.ht")
    create_clinvar_tb(
        input_path=FIXTURE,
        output_path=output_path,
        reference_genome="GRCh38",
        overwrite=True,
    )
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "ClinVar schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=keys)
    assert actual_rows == expected_rows, "ClinVar sample rows drifted from snapshot"
