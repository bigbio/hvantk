"""Round-trip and update tests for the GWAS Catalog builder skill."""

from __future__ import annotations

import json
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

FIXTURE = "hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv"
SNAPSHOT_DIR = Path("hvantk/skills/gwas_catalog/tests/snapshots")


@pytest.mark.hail
def test_gwas_catalog_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build GWAS Catalog from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.gwas_catalog.builder import build_gwas_catalog_associations

    builder = phase_b_snapshot_adapter(
        build_gwas_catalog_associations, "gwas_catalog:associations"
    )
    builder_kwargs = {"reference_genome": "GRCh38"}

    keys = json.loads((SNAPSHOT_DIR / "sample_keys.json").read_text())

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=keys,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "gwas_catalog.ht")
    builder(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "GWAS Catalog schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=keys)
    assert (
        actual_rows == expected_rows
    ), "GWAS Catalog sample rows drifted from snapshot"
