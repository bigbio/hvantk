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
    # The per-group VECTOR survives prepare alongside the declared roll-up. It used to be
    # dropped here -- entry.columns gated the collapse, so emit="vector" (the default) was
    # a no-op through the real pipeline no matter what the reducer emitted.
    assert set(prepared.row) == {"gene_id", "asp_cm", "asp_other", "asp_cm_spec"}
    assert sorted(prepared.gene_id.collect()) == ["ENSG_A", "ENSG_T"]  # OFFGENE dropped
    rows = {r.gene_id: r for r in prepared.collect()}
    assert rows["ENSG_T"].asp_cm_spec == pytest.approx(1.0)  # TNNT2: CM-specific
    assert rows["ENSG_A"].asp_cm_spec == pytest.approx(0.5)  # ACTB: ubiquitous
    # Vector columns carry the same EWCE fractions, per group rather than rolled up.
    assert rows["ENSG_T"].asp_cm == pytest.approx(1.0)
    assert rows["ENSG_T"].asp_other == pytest.approx(0.0)
    assert rows["ENSG_A"].asp_cm == pytest.approx(0.5)
    assert rows["ENSG_A"].asp_other == pytest.approx(0.5)


@pytest.mark.hail
def test_prepare_matrix_source_rejects_undeclared_column(hail_session):
    """A declared column the reducer never produces is an error, not a silent absence.

    entry.columns no longer gates which columns survive, so its remaining job is to catch
    a typo or an atlas whose group labels moved out from under the spec.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.annotation.prepare import prepare_matrix_source
    from hvantk.algorithms.annotation.spec import (
        MatrixSpec,
        SpecificitySpec,
        SourceEntry,
    )

    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": ["CM", "Other"]}, index=["CM", "Other"]),
        var=pd.DataFrame(index=["TNNT2"]),
        layers={"mean": np.array([[10.0], [0.0]])},
    )
    entry = SourceEntry(
        axis="expr",
        source="ucsc-cellbrowser:asp_2019",
        key="symbol",
        columns=("asp_typo_spec",),  # reducer emits asp_cm_spec, not this
        matrix=MatrixSpec(
            group_axis="celltype",
            atlas="asp",
            specificity=SpecificitySpec("ewce_fraction", ("CM",), "max", "cm_spec"),
        ),
    )
    fake = _FakeHGNC(
        canonical={"TNNT2": "TNNT2"},
        symbol_to_hgnc={"TNNT2": "HGNC:T"},
        hgnc_to_ensembl={"HGNC:T": "ENSG_T"},
    )
    with pytest.raises(ValueError, match="did not produce"):
        prepare_matrix_source(a, {"ENSG_T"}, entry, hgnc=fake)


@pytest.mark.hail
def test_prepare_matrix_source_collapses_symbol_collisions_onto_one_gene(hail_session):
    # Two distinct atlas symbols (an alias and its current symbol) resolve to ONE gene_id.
    # Unlike a plain gene table, a matrix must combine them (by max), not raise.
    import anndata as ad
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.annotation.prepare import prepare_matrix_source
    from hvantk.algorithms.annotation.spec import (
        MatrixSpec,
        SpecificitySpec,
        SourceEntry,
    )

    genes = ["TNNT2", "TNNT2_ALIAS"]  # both map to ENSG_T
    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": ["CM", "Other"]}, index=["CM", "Other"]),
        var=pd.DataFrame(index=genes),
        # TNNT2 fully CM-specific (spec 1.0); the alias row is ubiquitous (spec 0.5).
        layers={"mean": np.array([[10.0, 5.0], [0.0, 5.0]])},
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
    fake = _FakeHGNC(
        canonical={"TNNT2": "TNNT2", "TNNT2_ALIAS": "TNNT2_ALIAS"},
        symbol_to_hgnc={"TNNT2": "HGNC:T", "TNNT2_ALIAS": "HGNC:T"},
        hgnc_to_ensembl={"HGNC:T": "ENSG_T"},
    )
    prepared, report = prepare_matrix_source(a, {"ENSG_T"}, entry, hgnc=fake)
    assert list(prepared.key) == ["gene_id"]
    assert prepared.count() == 1  # collapsed to one row, not raised
    assert prepared.distinct().count() == 1
    d = {r.gene_id: r.asp_cm_spec for r in prepared.collect()}
    assert d["ENSG_T"] == pytest.approx(1.0)  # max(1.0 CM-specific, 0.5 ubiquitous)


def test_prepare_source_accepts_uniprot_id_key():
    """A uniprot_id-keyed gene-level source must reach the mapper, not be rejected early.

    `insider:interfaces` is protein-keyed: the reduction is per UniProt accession, and
    `GeneIdMapper.from_uniprot_ids` maps accession -> hgnc_id -> ensembl_gene_id. The
    id-space dispatch in `_resolve_to_gene_id` already handles it; this guards the guard
    in `prepare_source`, which listed the accepted keys separately and so silently
    excluded the space the schema advertises.
    """
    import pytest

    from hvantk.algorithms.annotation.prepare import prepare_source
    from hvantk.algorithms.annotation.spec import SourceEntry

    entry = SourceEntry(
        axis="ppi",
        source="insider:interfaces",
        key="uniprot_id",
        columns=("n_partners",),
    )
    # No hgnc streamer -> must fail on the MISSING STREAMER, not on the key space.
    with pytest.raises(ValueError, match="needs the HGNC streamer"):
        prepare_source(object(), [], entry, hgnc=None)


def test_prepare_source_still_rejects_an_unknown_key_space():
    import pytest

    from hvantk.algorithms.annotation.prepare import prepare_source
    from hvantk.algorithms.annotation.spec import SourceEntry

    entry = SourceEntry(axis="x", source="s", key="refseq_id", columns=("c",))
    with pytest.raises(ValueError, match="declares key"):
        prepare_source(object(), [], entry, hgnc=object())


def _by_gene_id_spec(collapse=None):
    from hvantk.algorithms.annotation.spec import AggregateSpec, ScoreSpec, SourceEntry

    return SourceEntry(
        axis="pqtl",
        source="pqtl:metrics",
        key="variant",
        columns=("n_pairs", "b_max"),
        collapse=collapse,
        aggregate=AggregateSpec(
            by="gene_id", to="gene_id", reduce="identity", count_name="n_pairs",
            scores=(ScoreSpec("b", "beta", ("max",)),),
        ),
    )


@pytest.mark.hail
def test_aggregate_by_gene_id_does_not_collide_with_the_key(hail_session):
    """A source that already carries gene_id must be aggregable by it.

    `_rekey_onto_gene_id` annotated `gene_id` onto the grouped table -- but when
    `aggregate.by` IS `gene_id`, group_by has already keyed the table on it, so Hail
    raises "cannot overwrite key field 'gene_id'". Both the pQTL and the GTEx eQTL
    sources are keyed that way, so this blocked two axes in all four configurations.
    """
    import hail as hl

    from hvantk.algorithms.annotation.prepare import prepare_variant_source

    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG_A", "beta": 0.4},
            {"gene_id": "ENSG_A", "beta": 0.9},
            {"gene_id": "ENSG_B", "beta": 0.1},
        ],
        hl.tstruct(gene_id=hl.tstr, beta=hl.tfloat64),
        key=["gene_id"],
    )
    prepared, _ = prepare_variant_source(ht, ["ENSG_A", "ENSG_B"], _by_gene_id_spec())
    d = {r.gene_id: r for r in prepared.collect()}
    assert d["ENSG_A"].n_pairs == 2
    assert d["ENSG_A"].b_max == pytest.approx(0.9)
    assert d["ENSG_B"].n_pairs == 1


@pytest.mark.hail
def test_many_to_one_mapping_raises_without_a_collapse_policy(hail_session):
    """Silent row loss must stay impossible unless the spec opts in."""
    import hail as hl

    from hvantk.algorithms.annotation.prepare import _rekey_onto_gene_id

    ht = hl.Table.parallelize(
        [{"uniprot_id": "P1", "v": 3}, {"uniprot_id": "P2", "v": 5}],
        hl.tstruct(uniprot_id=hl.tstr, v=hl.tint32),
        key=["uniprot_id"],
    )
    resolved = {"P1": "ENSG_A", "P2": "ENSG_A"}  # two accessions, one gene
    with pytest.raises(ValueError, match="many-to-one"):
        _rekey_onto_gene_id(ht, "uniprot_id", resolved, ("v",), "src", collapse=None)


@pytest.mark.hail
def test_collapse_max_reduces_a_many_to_one_mapping(hail_session):
    """`collapse: max` is the opt-in for id spaces that are legitimately many-to-one.

    One gene commonly has several UniProt accessions (isoforms, historical entries), so
    `insider:interfaces` and `uniprot-ptm:sites` both collapse onto the spine. `max` is
    the safe default reduction: the accessions describe the SAME protein, so it never
    double-counts the way `sum` would.
    """
    import hail as hl

    from hvantk.algorithms.annotation.prepare import _rekey_onto_gene_id

    ht = hl.Table.parallelize(
        [
            {"uniprot_id": "P1", "v": 3},
            {"uniprot_id": "P2", "v": 5},
            {"uniprot_id": "P3", "v": 7},
        ],
        hl.tstruct(uniprot_id=hl.tstr, v=hl.tint32),
        key=["uniprot_id"],
    )
    resolved = {"P1": "ENSG_A", "P2": "ENSG_A", "P3": "ENSG_B"}
    prepared = _rekey_onto_gene_id(
        ht, "uniprot_id", resolved, ("v",), "src", collapse="max"
    )
    d = {r.gene_id: r.v for r in prepared.collect()}
    assert d == {"ENSG_A": 5, "ENSG_B": 7}
