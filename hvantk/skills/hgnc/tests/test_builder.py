"""Snapshot round-trip test for the HGNC builder.

Builds from the committed fixture and asserts the schema and a sample of rows against
committed snapshots, so a change in builder behaviour or upstream field layout shows up
as an explicit diff rather than silently.

Regenerate after an intentional change:
    pytest hvantk/skills/hgnc/tests/test_builder.py -m hail --regenerate-snapshots
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

FIXTURE = "hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv"
SNAPSHOT_DIR = Path("hvantk/skills/hgnc/tests/snapshots")

# Real keys, taken from an actual build of the fixture -- collect_sample_rows raises
# KeyError for any key absent from the table, so these cannot be invented.
SAMPLE_KEYS = [
    {"hgnc_id": "HGNC:1100"},
    {"hgnc_id": "HGNC:1101"},
    {"hgnc_id": "HGNC:16376"},
    {"hgnc_id": "HGNC:4641"},
    {"hgnc_id": "HGNC:51839"},
]


@pytest.mark.hail
def test_hgnc_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build HGNC from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.hgnc.builder import build_hgnc_gene_lookup

    builder = phase_b_snapshot_adapter(build_hgnc_gene_lookup, "hgnc:lookup")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "hgnc.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert (
        hail_schema_to_dict(ht) == expected_schema
    ), "HGNC schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "HGNC sample rows drifted from snapshot"


@pytest.mark.hail
def test_map_from_hgnc_returns_empty_for_empty_input(hail_session, tmp_path):
    """map_from_hgnc([]) must return {} rather than crash on hl.literal(set())."""
    from hvantk.skills.hgnc.builder import build_hgnc_gene_lookup
    from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

    builder = phase_b_snapshot_adapter(build_hgnc_gene_lookup, "hgnc:lookup")
    out = str(tmp_path / "hgnc.ht")
    builder(input_path=FIXTURE, output_path=out)
    streamer = HGNCGeneCatalogStreamer.from_path(out)

    assert streamer.map_from_hgnc([], "ensembl_gene_id") == {}


@pytest.mark.hail
def test_quoted_multivalue_fields_are_unquoted_before_splitting(hail_session, tmp_path):
    """Regression: the real HGNC dump QUOTES multi-value fields; the fixture does not.

    ``hgnc_complete_set.txt`` writes e.g. ``prev_symbol`` as ``"H1F4|HIST1H1E"``. Splitting
    on ``|`` without stripping the quotes yields ``['"H1F4', 'HIST1H1E"']`` -- every first
    and last element carries a stray quote, so ``prev_symbols``/``alias_symbols`` lookups
    for the clean symbol miss and the gene silently fails to map. The committed fixture
    has zero quotes, which is why this never surfaced.

    Real-world impact this reproduces: HIST1H1E (previous symbol of H1-4) did not resolve,
    so 5 CHD cohort genes dropped out at mapping time. ``uniprot_ids`` and ``gene_group``
    are corrupted the same way.
    """
    import hail as hl
    from hvantk.skills.hgnc.builder import build_hgnc_gene_lookup

    header = (
        "hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\t"
        "location_sortable\talias_symbol\talias_name\tprev_symbol\tprev_name\t"
        "gene_group\tgene_group_id\tdate_approved_reserved\tdate_symbol_changed\t"
        "date_name_changed\tdate_modified\tentrez_id\tensembl_gene_id\tvega_id\t"
        "ucsc_id\tena\trefseq_accession\tccds_id\tuniprot_ids\n"
    )
    # Quoted exactly as the real dump writes them.
    row = (
        "HGNC:4718\tH1-4\tH1.4 linker histone\tprotein-coding gene\t"
        "gene with protein product\tApproved\t6p22.2\t06p22.2\t"
        '"H1.4|H1e"\t\t"H1F4|HIST1H1E"\t\t"grp1|grp2"\t"1|2"\t\t\t\t\t\t'
        'ENSG00000168298\t\t\t\t\t\t"P10412|Q4VB24"\n'
    )
    src = tmp_path / "hgnc_quoted.tsv"
    src.write_text(header + row)

    builder = phase_b_snapshot_adapter(build_hgnc_gene_lookup, "hgnc:lookup")
    out = str(tmp_path / "hgnc.ht")
    builder(input_path=Path(src).resolve().as_uri(), output_path=out)
    r = hl.read_table(out).collect()[0]

    assert list(r.prev_symbols) == ["H1F4", "HIST1H1E"], list(r.prev_symbols)
    assert list(r.alias_symbols) == ["H1.4", "H1e"], list(r.alias_symbols)
    assert list(r.uniprot_ids) == ["P10412", "Q4VB24"], list(r.uniprot_ids)
    assert list(r.gene_group) == ["grp1", "grp2"], list(r.gene_group)
