"""The gene spine: one row per protein-coding gene, keyed on Ensembl gene_id."""
from __future__ import annotations

import pytest


def _genes_ht():
    import hail as hl

    return hl.Table.parallelize(
        [
            {
                "gene_id": "ENSG1",
                "gene_name": "MYL7",
                "chromosome": "7",
                "gene_start": 100,
                "gene_end": 200,
            },
            {
                "gene_id": "ENSG2",
                "gene_name": "LNC1",
                "chromosome": "2",
                "gene_start": 300,
                "gene_end": 400,
            },
            {
                "gene_id": "ENSG3",
                "gene_name": "ORPH",
                "chromosome": "3",
                "gene_start": 500,
                "gene_end": 600,
            },
        ],
        hl.tstruct(
            gene_id=hl.tstr,
            gene_name=hl.tstr,
            chromosome=hl.tstr,
            gene_start=hl.tint32,
            gene_end=hl.tint32,
        ),
        key=["gene_id"],
    )


def _structure_ht():
    import hail as hl

    return hl.Table.parallelize(
        [
            {
                "gene_id": "ENSG1",
                "gene_biotype": "protein_coding",
                "mane_select": "ENST1",
                "cds_length": 300,
                "n_coding_exons": 2,
                "n_transcripts": 2,
            },
            {
                "gene_id": "ENSG2",
                "gene_biotype": "lncRNA",
                "mane_select": "",
                "cds_length": 0,
                "n_coding_exons": 0,
                "n_transcripts": 1,
            },
            {
                "gene_id": "ENSG3",
                "gene_biotype": "protein_coding",
                "mane_select": "",
                "cds_length": 150,
                "n_coding_exons": 1,
                "n_transcripts": 1,
            },
        ],
        hl.tstruct(
            gene_id=hl.tstr,
            gene_biotype=hl.tstr,
            mane_select=hl.tstr,
            cds_length=hl.tint32,
            n_coding_exons=hl.tint32,
            n_transcripts=hl.tint32,
        ),
        key=["gene_id"],
    )


def _hgnc_ht():
    """ENSG3 is deliberately absent -- a protein-coding gene with no HGNC record."""
    import hail as hl

    return hl.Table.parallelize(
        [
            {
                "hgnc_id": "HGNC:1",
                "symbol": "MYL7",
                "ensembl_gene_id": "ENSG1",
                "gene_group": "Myosin light chains",
            },
            {
                "hgnc_id": "HGNC:2",
                "symbol": "LNC1",
                "ensembl_gene_id": "ENSG2",
                "gene_group": "",
            },
        ],
        hl.tstruct(
            hgnc_id=hl.tstr, symbol=hl.tstr, ensembl_gene_id=hl.tstr, gene_group=hl.tstr
        ),
        key=["hgnc_id"],
    )


@pytest.mark.hail
def test_spine_keeps_only_protein_coding(hail_session):
    from hvantk.algorithms.annotation.spine import build_spine

    spine = build_spine(_genes_ht(), _structure_ht(), _hgnc_ht())
    assert sorted(spine.gene_id.collect()) == ["ENSG1", "ENSG3"]


@pytest.mark.hail
def test_gene_id_is_the_key_and_is_unique(hail_session):
    from hvantk.algorithms.annotation.spine import build_spine

    spine = build_spine(_genes_ht(), _structure_ht(), _hgnc_ht())
    assert list(spine.key) == ["gene_id"]
    assert spine.count() == spine.distinct().count()


@pytest.mark.hail
def test_structure_and_hgnc_fields_are_carried(hail_session):
    from hvantk.algorithms.annotation.spine import build_spine

    spine = build_spine(_genes_ht(), _structure_ht(), _hgnc_ht())
    row = spine.filter(spine.gene_id == "ENSG1").collect()[0]

    assert row.cds_length == 300
    assert row.n_coding_exons == 2
    assert row.mane_select == "ENST1"
    assert row.hgnc_id == "HGNC:1"
    assert row.gene_group == "Myosin light chains"


@pytest.mark.hail
def test_a_gene_without_an_hgnc_record_is_kept_with_a_missing_hgnc_id(hail_session):
    """The join must not silently drop genes HGNC does not cover."""
    from hvantk.algorithms.annotation.spine import build_spine

    spine = build_spine(_genes_ht(), _structure_ht(), _hgnc_ht())
    row = spine.filter(spine.gene_id == "ENSG3").collect()[0]

    assert row.hgnc_id is None
    assert row.cds_length == 150


@pytest.mark.hail
def test_hgnc_mapping_rate_is_reported(hail_session):
    from hvantk.algorithms.annotation.spine import build_spine, spine_mapping_rate

    spine = build_spine(_genes_ht(), _structure_ht(), _hgnc_ht())
    assert spine_mapping_rate(spine) == pytest.approx(0.5)  # 1 of 2 kept genes


@pytest.mark.hail
def test_two_hgnc_records_for_one_gene_do_not_duplicate_the_row(hail_session):
    """Join arity must be invariant: HGNC can map two hgnc_ids to one Ensembl gene.

    A left join on a non-unique right side silently multiplies spine rows, which would
    then multiply every feature joined on afterwards. The pick must also be
    deterministic, so a rebuild does not change which hgnc_id a gene carries.
    """
    import hail as hl

    from hvantk.algorithms.annotation.spine import build_spine

    dup_hgnc = hl.Table.parallelize(
        [
            {
                "hgnc_id": "HGNC:9",
                "symbol": "MYL7B",
                "ensembl_gene_id": "ENSG1",
                "gene_group": "Z group",
            },
            {
                "hgnc_id": "HGNC:1",
                "symbol": "MYL7",
                "ensembl_gene_id": "ENSG1",
                "gene_group": "Myosin light chains",
            },
        ],
        hl.tstruct(
            hgnc_id=hl.tstr, symbol=hl.tstr, ensembl_gene_id=hl.tstr, gene_group=hl.tstr
        ),
        key=["hgnc_id"],
    )

    spine = build_spine(_genes_ht(), _structure_ht(), dup_hgnc)

    assert spine.count() == 2  # ENSG1 + ENSG3, not 3
    assert spine.count() == spine.distinct().count()
    row = spine.filter(spine.gene_id == "ENSG1").collect()[0]
    assert row.hgnc_id == "HGNC:1"  # lowest hgnc_id wins, deterministically
    assert row.gene_group == "Myosin light chains"
