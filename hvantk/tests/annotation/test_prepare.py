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
def test_a_source_that_maps_no_genes_raises_a_clear_error(hail_session):
    from hvantk.algorithms.annotation.mapping import MappingRateError
    from hvantk.algorithms.annotation.prepare import prepare_source

    # The spine shares no gene_id with the source -> every row is unmapped.
    with pytest.raises(MappingRateError, match="0 of 3"):
        prepare_source(_source_ht(), {"ENSG_NONE"}, ENTRY)


class _FakeHGNC:
    """Minimal stand-in for HGNCGeneCatalogStreamer: only the methods GeneIdMapper calls."""

    def __init__(self, hgnc_to_ensembl=None, symbol_to_hgnc=None, canonical=None):
        self._h2e = hgnc_to_ensembl or {}
        self._s2h = symbol_to_hgnc or {}
        self._canon = canonical or {}

    def map_from_hgnc(self, ids, field):
        return {i: self._h2e.get(i) for i in ids}

    def resolve_to_canonical(self, symbol):
        return self._canon.get(symbol, symbol)

    def map_to_hgnc(self, symbols, field):
        return {s: self._s2h.get(s) for s in symbols}


def _hgnc_source_ht():
    import hail as hl

    return hl.Table.parallelize(
        [
            {"hgnc_id": "HGNC:1", "score": 1.0},
            {"hgnc_id": "HGNC:2", "score": 2.0},
            {"hgnc_id": "HGNC:404", "score": 9.0},  # will not resolve onto the spine
        ],
        hl.tstruct(hgnc_id=hl.tstr, score=hl.tfloat64),
        key=["hgnc_id"],
    )


HGNC_ENTRY = SourceEntry(axis="x", source="s:hg", key="hgnc_id", columns=("score",))


@pytest.mark.hail
def test_hgnc_id_source_is_rekeyed_onto_gene_id(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    fake = _FakeHGNC(
        hgnc_to_ensembl={"HGNC:1": "ENSG1", "HGNC:2": "ENSG2", "HGNC:404": "ENSGX"}
    )
    prepared, report = prepare_source(
        _hgnc_source_ht(), {"ENSG1", "ENSG2"}, HGNC_ENTRY, hgnc=fake
    )
    assert list(prepared.key) == ["gene_id"]
    assert set(prepared.row) == {"gene_id", "score"}
    assert sorted(prepared.gene_id.collect()) == [
        "ENSG1",
        "ENSG2",
    ]  # HGNC:404 off-spine dropped
    assert report.n_in == 3 and report.n_mapped == 2


@pytest.mark.hail
def test_symbol_source_resolves_aliases_before_mapping(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.prepare import prepare_source

    src = hl.Table.parallelize(
        [{"symbol": "OLDNAME", "score": 1.0}],
        hl.tstruct(symbol=hl.tstr, score=hl.tfloat64),
        key=["symbol"],
    )
    # OLDNAME is a previous symbol for the gene whose current symbol is NEWNAME.
    fake = _FakeHGNC(
        canonical={"OLDNAME": "NEWNAME"},
        symbol_to_hgnc={"NEWNAME": "HGNC:1"},
        hgnc_to_ensembl={"HGNC:1": "ENSG1"},
    )
    entry = SourceEntry(axis="x", source="s:sym", key="symbol", columns=("score",))
    prepared, report = prepare_source(src, {"ENSG1"}, entry, hgnc=fake)
    assert prepared.gene_id.collect() == ["ENSG1"]
    assert report.n_mapped == 1


@pytest.mark.hail
def test_a_non_gene_id_key_without_hgnc_raises(hail_session):
    from hvantk.algorithms.annotation.prepare import prepare_source

    with pytest.raises(ValueError, match="hgnc"):
        prepare_source(_hgnc_source_ht(), {"ENSG1"}, HGNC_ENTRY, hgnc=None)


@pytest.mark.hail
def test_two_keys_mapping_to_one_gene_is_rejected(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.prepare import prepare_source

    src = hl.Table.parallelize(
        [{"hgnc_id": "HGNC:1", "score": 1.0}, {"hgnc_id": "HGNC:2", "score": 2.0}],
        hl.tstruct(hgnc_id=hl.tstr, score=hl.tfloat64),
        key=["hgnc_id"],
    )
    # Both HGNC ids resolve to the same gene -> the re-keyed table would have two rows for ENSG1.
    fake = _FakeHGNC(hgnc_to_ensembl={"HGNC:1": "ENSG1", "HGNC:2": "ENSG1"})
    with pytest.raises(ValueError, match="same gene_id"):
        prepare_source(
            src, {"ENSG1"}, SourceEntry("x", "s:hg", "hgnc_id", ("score",)), hgnc=fake
        )


def _tol_entry():
    from hvantk.algorithms.annotation.spec import AggregateSpec, ScoreSpec, SourceEntry

    agg = AggregateSpec(
        by="Ensembl_geneid",
        to="gene_id",
        filter="missense",
        reduce="max",
        scores=(ScoreSpec("revel", "REVEL_score", ("mean", "max", "frac_gt_0.5")),),
    )
    return SourceEntry(
        axis="tolerance",
        source="dbnsfp:variants",
        key="variant",
        columns=("revel_mean", "revel_max", "revel_frac_gt_0.5", "n_possible_missense"),
        aggregate=agg,
    )


@pytest.mark.hail
def test_prepare_variant_source_aggregates_and_maps_onto_spine(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.prepare import prepare_variant_source

    # variant table: 2 missense in ENSG_A (on spine), 1 in ENSG_OFF (not on spine).
    rows = [
        {
            "locus": hl.locus("chr1", 100, "GRCh38"),
            "alleles": ["A", "C"],
            "aaref": "M",
            "aaalt": "T",
            "Ensembl_geneid": "ENSG_A",
            "REVEL_score": {"t1": 0.9},
        },
        {
            "locus": hl.locus("chr1", 200, "GRCh38"),
            "alleles": ["G", "T"],
            "aaref": "R",
            "aaalt": "Q",
            "Ensembl_geneid": "ENSG_A",
            "REVEL_score": {"t1": 0.1},
        },
        {
            "locus": hl.locus("chr1", 300, "GRCh38"),
            "alleles": ["C", "A"],
            "aaref": "D",
            "aaalt": "N",
            "Ensembl_geneid": "ENSG_OFF",
            "REVEL_score": {"t1": 0.5},
        },
    ]
    src = hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            aaref=hl.tstr,
            aaalt=hl.tstr,
            Ensembl_geneid=hl.tstr,
            REVEL_score=hl.tdict(hl.tstr, hl.tfloat64),
        ),
        key=["locus", "alleles"],
    )
    prepared, report = prepare_variant_source(src, {"ENSG_A"}, _tol_entry())
    assert list(prepared.key) == ["gene_id"]
    assert set(prepared.row) == {
        "gene_id",
        "revel_mean",
        "revel_max",
        "revel_frac_gt_0.5",
        "n_possible_missense",
    }
    assert prepared.gene_id.collect() == ["ENSG_A"]  # ENSG_OFF dropped (off spine)
    row = prepared.collect()[0]
    assert row.n_possible_missense == 2 and row.revel_max == pytest.approx(0.9)
    assert report.n_mapped == 1 and report.n_in == 2  # 2 aggregated genes; 1 on spine


@pytest.mark.hail
def test_prepare_source_gene_id_path_unchanged(hail_session):
    # Guard the refactor: the direct gene_id path still behaves as in P2c-2.
    from hvantk.algorithms.annotation.prepare import prepare_source

    prepared, report = prepare_source(_source_ht(), {"ENSG1", "ENSG2"}, ENTRY)
    assert list(prepared.key) == ["gene_id"]
    assert sorted(prepared.gene_id.collect()) == ["ENSG1", "ENSG2"]


@pytest.mark.hail
def test_prepare_matrix_source_reduces_and_maps_onto_spine(hail_session):
    import anndata as ad
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.annotation.prepare import prepare_matrix_source
    from hvantk.algorithms.annotation.spec import (
        MatrixSpec,
        SpecificitySpec,
        SourceEntry,
    )

    genes = ["TNNT2", "ACTB", "OFFGENE"]
    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": ["CM", "Other"]}, index=["CM", "Other"]),
        var=pd.DataFrame(index=genes),
        layers={"mean": np.array([[10.0, 5.0, 1.0], [0.0, 5.0, 1.0]])},
    )
    entry = SourceEntry(
        axis="expr",
        source="ucsc-cellbrowser:asp_2019",
        key="symbol",
        columns=("asp_cm_spec",),
        matrix=MatrixSpec(
            group_axis="celltype",
            atlas="asp",
            specificity=SpecificitySpec("ewce_fraction", ("CM",), "max", "cm_spec"),
        ),
    )
    # HGNC fake: TNNT2/ACTB resolve onto the spine; OFFGENE resolves off-spine.
    fake = _FakeHGNC(
        canonical={"TNNT2": "TNNT2", "ACTB": "ACTB", "OFFGENE": "OFFGENE"},
        symbol_to_hgnc={"TNNT2": "HGNC:T", "ACTB": "HGNC:A", "OFFGENE": "HGNC:O"},
        hgnc_to_ensembl={"HGNC:T": "ENSG_T", "HGNC:A": "ENSG_A", "HGNC:O": "ENSG_OFF"},
    )
    prepared, report = prepare_matrix_source(a, {"ENSG_T", "ENSG_A"}, entry, hgnc=fake)
    assert list(prepared.key) == ["gene_id"]
    assert set(prepared.row) == {"gene_id", "asp_cm_spec"}
    assert sorted(prepared.gene_id.collect()) == ["ENSG_A", "ENSG_T"]  # OFFGENE dropped
    d = {r.gene_id: r.asp_cm_spec for r in prepared.collect()}
    assert d["ENSG_T"] == pytest.approx(1.0)  # TNNT2: CM-specific
    assert d["ENSG_A"] == pytest.approx(0.5)  # ACTB: ubiquitous
