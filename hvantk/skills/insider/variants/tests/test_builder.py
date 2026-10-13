"""Round-trip and update tests for the INSIDER (BED) builder skill."""

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

FIXTURE = (
    "hvantk/skills/insider/variants/tests/testdata/raw/insider/"
    "insider_sample-synthetic.bed"
)
SNAPSHOT_DIR = Path("hvantk/skills/insider/variants/tests/snapshots")

# `interval` is unique-in-table after .distinct() (per skill §5), so keys are
# inlined here. The fixture is synthetic (#426 — no published INSIDER licence was
# found; see testdata/raw/insider/README.md). Picked from observed fixture rows:
# chr3:2000301-2000303 is shared across all four SYNTH0003_ppi_* tracks (exercises
# the multi-PPI-per-interval aggregation), chr11:1000101-1000102 is one of the
# first track's three surviving rows (its opening row is zero-length and gets
# skipped), and chr3:2000711-2000713 is from the last track. Intervals stored as
# {"start": "<contig>:<pos>", "end": "<contig>:<pos>"} per the tinterval branch in
# _snapshot_utils.
SAMPLE_KEYS = [
    {"interval": {"start": "chr3:2000301", "end": "chr3:2000303"}},
    {"interval": {"start": "chr11:1000101", "end": "chr11:1000102"}},
    {"interval": {"start": "chr3:2000711", "end": "chr3:2000713"}},
]


@pytest.mark.hail
def test_insider_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build INSIDER BED from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.insider.variants.builder import build_insider_interactome

    builder = phase_b_snapshot_adapter(build_insider_interactome, "insider:variants")
    builder_kwargs = {"reference_genome": "GRCh38"}

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "insider.ht")
    fixture_uri = Path(FIXTURE).resolve().as_uri()
    builder(input_path=fixture_uri, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "INSIDER schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "INSIDER sample rows drifted from snapshot"


@pytest.mark.hail
def test_insider_drops_invalid_loci(hail_session, tmp_path):
    """Regression (DEFECT 1): out-of-range / non-reference-contig rows are
    dropped instead of aborting the build.

    The builder previously called ``hl.locus_interval`` with the default
    ``invalid_missing=False``, so a single BED row whose locus fell outside the
    reference bounds or on a non-reference contig raised and aborted the whole
    build. With ``invalid_missing=True`` plus an ``hl.is_defined`` filter, such
    rows are silently dropped (matching the prior builder's
    ``skip_invalid_intervals=True`` intent) and the build succeeds.

    Hail-dependent: must run on Slurm
    (``pytest hvantk/skills/insider/variants/tests -m hail``); the login node OOMs on
    Hail init, so this test is intentionally not executed there.
    """
    import hail as hl
    from hvantk.skills.insider.variants.builder import build_insider_interactome

    # One in-bounds chr21 row, plus two rows that would abort the old builder:
    # a non-reference contig and an out-of-range coordinate (chr21 on GRCh38 is
    # 46,709,983 bp long).
    bed = tmp_path / "insider_invalid.bed"
    bed.write_text(
        'track name=P1_ppi_P2 description="x" visibility=dense itemRgb="On"\n'
        "chr21\t5030000\t5030002\t.\t0\t+\n"
        "chrBOGUS\t100\t102\t.\t0\t+\n"
        "chr21\t900000000\t900000002\t.\t0\t+\n"
    )

    builder = phase_b_snapshot_adapter(build_insider_interactome, "insider:variants")
    output_path = str(tmp_path / "insider_invalid.ht")
    builder(
        input_path=bed.resolve().as_uri(),
        output_path=output_path,
        reference_genome="GRCh38",
    )

    ht = hl.read_table(output_path)
    # Only the single in-bounds chr21 interval survives; the two invalid rows
    # are dropped rather than aborting the build.
    assert ht.count() == 1
    row = ht.collect()[0]
    assert row.interval.start.contig == "chr21"
    assert row.interval.start.position == 5030001
    assert row.ppi_ids == ["P1_ppi_P2"]
