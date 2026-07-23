"""Snapshot round-trip test for the UniProt PTM sites builder.

Builds from the committed fixture and asserts the schema and a sample of rows against
committed snapshots, so a change in builder behaviour or upstream field layout shows up
as an explicit diff rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/uniprot_ptm/tests/test_builder.py -m hail --regenerate-snapshots
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

FIXTURE = "hvantk/skills/uniprot_ptm/tests/testdata/raw/uniprot-ptm/ptm_sites.tsv"
SNAPSHOT_DIR = Path("hvantk/skills/uniprot_ptm/tests/snapshots")

# Real keys, taken from an actual build of the fixture -- collect_sample_rows raises
# KeyError for any key absent from the table, so these cannot be invented.
#
# The table is keyed by `locus` alone (builder.py), derived as
# hl.locus("chr" + chrom, codon_start) after the MT -> M remap. The fixture's
# GL000009.2 row is deliberately absent here: "chrGL000009.2" is not a GRCh38 contig, so
# the builder filters it out. That row exists to cover the filter, not to be sampled.
SAMPLE_KEYS = [
    {"locus": "chr1:11796321"},
    {"locus": "chr13:32316461"},
    {"locus": "chr17:7675994"},
    {"locus": "chr17:7676272"},
    {"locus": "chrM:3308"},
]


@pytest.mark.hail
def test_uniprot_ptm_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build UniProt PTM sites from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.uniprot_ptm.builder import build_uniprot_ptm_sites

    builder = phase_b_snapshot_adapter(build_uniprot_ptm_sites, "uniprot-ptm:sites")
    builder_kwargs = {"reference_genome": "GRCh38", "flanking_codons": 5}

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "uniprot_ptm.ht")
    builder(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    # The non-reference contig row must not survive the GRCh38 contig filter.
    assert ht.count() == 5, "expected the GL000009.2 row to be filtered out"

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "UniProt PTM schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "UniProt PTM sample rows drifted from snapshot"
