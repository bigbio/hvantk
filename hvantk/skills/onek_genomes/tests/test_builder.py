"""Snapshot round-trip tests for the 1000 Genomes builders.

Builds both datasets from the committed fixture and asserts schema and sample rows
against committed snapshots, so a change in builder behaviour or upstream field layout
shows up as an explicit diff rather than silently.

This plugin shipped four snapshot JSONs and two drift fingerprints but no test module at
all, so ``plugin.yaml``'s declared ``pytest hvantk/skills/onek_genomes/tests -m hail``
collected nothing and six committed artifacts were read by no one. The manifest loader
resolves those paths without requiring them (``TestPaths`` calls ``.resolve()``, which
does not check existence), so nothing failed.

Regenerate after an intentional change:
    pytest hvantk/skills/onek_genomes/tests -m hail --regenerate-snapshots

Note these snapshots predate ``_snapshot_utils.regenerate_snapshots`` and use a different
convention from the other plugins: per-dataset filenames (``variants_schema.json`` rather
than ``schema.json``) and a flat head-of-N row sample rather than the key-matched
``{"key": ..., "row": ...}`` shape. The regeneration path below preserves that convention
rather than rewriting four committed artifacts to match a helper.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    _to_jsonable,
    hail_schema_to_dict,
    load_snapshot,
    phase_b_snapshot_adapter,
)

FIXTURE = "hvantk/skills/onek_genomes/tests/testdata/raw/onek_genomes"
SNAPSHOT_DIR = Path("hvantk/skills/onek_genomes/tests/snapshots")

#: Rows sampled per dataset. The fixture is a single sorted chr22 VCF and a sorted
#: sample table, so a head sample is deterministic across runs.
N_SAMPLE_ROWS = 5


def _head_rows(table, n: int = N_SAMPLE_ROWS) -> list[dict]:
    """Collect the first ``n`` rows of a Hail Table as JSON-stable dicts."""
    return [
        {k: _to_jsonable(v) for k, v in dict(row).items()}
        for row in table.head(n).collect()
    ]


def _assert_or_regenerate(rows_table, schema_obj, stem: str, regenerate: bool) -> None:
    schema_path = SNAPSHOT_DIR / f"{stem}_schema.json"
    rows_path = SNAPSHOT_DIR / f"{stem}_sample_rows.json"

    schema = hail_schema_to_dict(schema_obj)
    rows = _head_rows(rows_table)

    if regenerate:
        SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
        schema_path.write_text(json.dumps(schema, indent=2, sort_keys=True))
        rows_path.write_text(json.dumps(rows, indent=2, sort_keys=True))
        return

    assert schema == load_snapshot(schema_path), f"{stem} schema drifted from snapshot"
    assert rows == load_snapshot(rows_path), f"{stem} sample rows drifted from snapshot"


@pytest.mark.hail
def test_onek_genomes_variants_snapshot_round_trip(
    hail_session, tmp_path, regenerate_snapshots
):
    """Build onek-genomes:variants from fixture; assert schema and sample rows."""
    import hail as hl
    from hvantk.skills.onek_genomes.builder import build_onek_genomes_variants

    builder = phase_b_snapshot_adapter(
        build_onek_genomes_variants, "onek-genomes:variants"
    )
    output_path = str(tmp_path / "onek_variants.mt")
    builder(input_path=FIXTURE, output_path=output_path, chromosomes=["chr22"])
    mt = hl.read_matrix_table(output_path)

    _assert_or_regenerate(mt.rows(), mt, "variants", regenerate_snapshots)
    if regenerate_snapshots:
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots.")


@pytest.mark.hail
def test_onek_genomes_samples_snapshot_round_trip(
    hail_session, tmp_path, regenerate_snapshots
):
    """Build onek-genomes:samples from fixture; assert schema and sample rows."""
    import hail as hl
    from hvantk.skills.onek_genomes.samples_builder import build_onek_genomes_samples

    builder = phase_b_snapshot_adapter(
        build_onek_genomes_samples, "onek-genomes:samples"
    )
    output_path = str(tmp_path / "onek_samples.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    _assert_or_regenerate(ht, ht, "samples", regenerate_snapshots)
    if regenerate_snapshots:
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots.")
