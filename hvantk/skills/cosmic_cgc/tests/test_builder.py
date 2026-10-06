"""Round-trip test for the COSMIC CGC submissions builder.

Builds from the synthetic fixture (see the README beside
``hvantk/skills/cosmic_cgc/tests/testdata/raw/cosmic-cgc/``) and asserts the
schema and a sample of rows against committed snapshots, so a change in
builder behaviour or upstream header layout shows up as an explicit diff
rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/cosmic_cgc/tests/test_builder.py -m hail --regenerate-snapshots
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

FIXTURE = (
    "hvantk/skills/cosmic_cgc/tests/testdata/raw/cosmic-cgc/cosmic-cgc-synthetic.tsv.gz"
)
SNAPSHOT_DIR = Path("hvantk/skills/cosmic_cgc/tests/snapshots")

# hgnc_id keys for SYNTHA, SYNTHB, SYNTHD (see the fixture README's row-purpose
# table). SYNTHC (every multi-value field empty -> [] arrays) is still built
# and exercised, just not individually snapshotted here; SYNTHB already
# carries one empty multi-value field, so the []-vs-missing path is covered
# by the three keys below. SYNTHE has no hgnc_id -- it is dropped by the stub
# catalog and asserted absent instead of sampled.
SAMPLE_KEYS = [
    {"hgnc_id": "900001"},  # SYNTHA
    {"hgnc_id": "900002"},  # SYNTHB
    {"hgnc_id": "900004"},  # SYNTHD
]


class _StubGeneCatalog:
    """Minimal gene_catalog stub -- duck-typed, not a GeneCatalogStreamer subclass.

    The builder only ever calls ``.map_ids(...)`` on whatever object is passed
    as ``gene_catalog``, so a plain stub is sufficient. Maps SYNTHA-D to
    invented HGNC ids; SYNTHE is deliberately absent so the builder's
    ``hl.is_defined(hgnc_id) & (hgnc_id != "")`` filter drops that row.
    """

    _FULL_MAP = {
        "SYNTHA": "HGNC:900001",
        "SYNTHB": "HGNC:900002",
        "SYNTHC": "HGNC:900003",
        "SYNTHD": "HGNC:900004",
    }

    def map_ids(self, ids, source_type, target_type):
        assert (source_type, target_type) == ("gene_symbol", "hgnc_id"), (
            "build_cosmic_cgc_submissions must resolve gene_symbol -> hgnc_id, "
            f"got {source_type!r} -> {target_type!r}"
        )
        # Restricted to the symbols actually passed, like a real catalog would be.
        return {k: v for k, v in self._FULL_MAP.items() if k in ids}


@pytest.mark.hail
def test_cosmic_cgc_submissions_round_trip(
    hail_session, tmp_path, regenerate_snapshots
):
    """Build COSMIC CGC from the synthetic fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.cosmic_cgc.builder import build_cosmic_cgc_submissions

    builder = phase_b_snapshot_adapter(
        build_cosmic_cgc_submissions, "cosmic-cgc:submissions"
    )
    builder_kwargs = {
        "mutation_context": "both",
        "gene_catalog": _StubGeneCatalog(),
    }

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

    output_path = str(tmp_path / "cosmic_cgc.ht")
    builder(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    # SYNTHE has no entry in the stub catalog's map -> dropped when keyed by hgnc_id.
    kept_symbols = ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol))
    assert "SYNTHE" not in kept_symbols, "unresolved SYNTHE row must be dropped"
    assert ht.count() == 4, "expected SYNTHA-D to survive, SYNTHE dropped"

    # The current (v103+) header's coordinate columns must be renamed and cast,
    # and the raw UPPER_CASE header names must not leak into the schema.
    row_fields = set(ht.row.dtype.fields)
    assert ht.row.dtype["genome_start"] == hl.tint32
    assert ht.row.dtype["genome_stop"] == hl.tint32
    raw_upper_case = {"COSMIC_GENE_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP"}
    assert not (raw_upper_case & row_fields), (
        f"raw upper-case columns leaked into the schema: {raw_upper_case & row_fields}"
    )

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, (
        "COSMIC CGC schema drifted from snapshot"
    )

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "COSMIC CGC sample rows drifted from snapshot"
