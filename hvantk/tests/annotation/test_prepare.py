"""prepare_source: map a built source onto the spine's gene_id and select columns."""
from __future__ import annotations

import pytest

from hvantk.algorithms.annotation.spec import SourceEntry


def _source_ht():
    import hail as hl

    # 3 source genes: ENSG1/ENSG2 on the spine, ENSG_OFF not on it.
    return hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "mis_z": 2.5, "pLI": 0.99, "extra": "drop-me"},
            {"gene_id": "ENSG2", "mis_z": -1.0, "pLI": 0.01, "extra": "drop-me"},
            {"gene_id": "ENSG_OFF", "mis_z": 9.9, "pLI": 0.5, "extra": "drop-me"},
        ],
        hl.tstruct(gene_id=hl.tstr, mis_z=hl.tfloat64, pLI=hl.tfloat64, extra=hl.tstr),
        key=["gene_id"],
    )


ENTRY = SourceEntry(
    axis="constraint",
    source="gnomad-metrics:metrics",
    key="gene_id",
    columns=("mis_z", "pLI"),
    min_mapping_rate=0.5,
)


@pytest.mark.hail
def test_prepared_table_is_gene_id_keyed_with_only_declared_columns(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    prepared, report = prepare_source(_source_ht(), {"ENSG1", "ENSG2"}, ENTRY)
    assert list(prepared.key) == ["gene_id"]
    assert set(prepared.row) == {"gene_id", "mis_z", "pLI"}  # "extra" dropped


@pytest.mark.hail
def test_rows_off_the_spine_are_dropped(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    prepared, report = prepare_source(_source_ht(), {"ENSG1", "ENSG2"}, ENTRY)
    assert sorted(prepared.gene_id.collect()) == ["ENSG1", "ENSG2"]  # ENSG_OFF gone


@pytest.mark.hail
def test_mapping_report_counts_the_off_spine_row_as_unmapped(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    prepared, report = prepare_source(_source_ht(), {"ENSG1", "ENSG2"}, ENTRY)
    assert report.n_in == 3
    assert report.n_mapped == 2
    assert report.n_unmapped == 1
    assert "ENSG_OFF" in report.unmapped
    assert report.rate == pytest.approx(2 / 3)


@pytest.mark.hail
def test_a_non_gene_id_key_is_rejected_in_p2a(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    bad = SourceEntry("x", "s:d", "hgnc_id", ("mis_z",))
    with pytest.raises(ValueError, match="gene_id"):
        prepare_source(_source_ht(), {"ENSG1"}, bad)
