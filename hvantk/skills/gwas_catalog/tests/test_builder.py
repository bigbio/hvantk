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

FIXTURE = (
    "hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv"
)
SNAPSHOT_DIR = Path("hvantk/skills/gwas_catalog/tests/snapshots")


def test_a_raw_directory_resolves_to_its_associations_file(tmp_path):
    """Regression guard for #388: ``hvantk reprocess`` with a plain ``--raw-dir``
    (no ``--intermediate``) hands the builder the raw directory, not the file
    inside it. Hail-free on purpose, so it runs in the default suite rather
    than only on a Hail node."""
    # The shared resolver, with the exact patterns `build_gwas_catalog_associations`
    # passes it.
    from hvantk.core.utils.raw_inputs import resolve_raw_file
    from hvantk.skills.gwas_catalog.builder import _RAW_FILE_PATTERNS

    tsv = tmp_path / "gwas-catalog-download-associations-v1.0-full.tsv"
    tsv.write_text("dummy")

    assert resolve_raw_file(tmp_path, patterns=_RAW_FILE_PATTERNS) == str(tsv)


@pytest.mark.round_trip("gwas-catalog:associations")
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
    assert actual_rows == expected_rows, (
        "GWAS Catalog sample rows drifted from snapshot"
    )
