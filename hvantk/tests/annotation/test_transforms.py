"""transforms.py: stat-token vocabulary + gene aggregation."""
from __future__ import annotations

import pytest


def test_output_name_is_score_underscore_token():
    from hvantk.algorithms.annotation.transforms import output_name

    assert output_name("revel", "mean") == "revel_mean"
    assert output_name("revel", "frac_gt_0.5") == "revel_frac_gt_0.5"


def test_unknown_stat_token_is_rejected():
    import hail as hl  # only for a scalar expr to pass in

    from hvantk.algorithms.annotation.transforms import _agg_for

    with pytest.raises(ValueError, match="unknown stat token"):
        _agg_for("median", hl.float64(1.0))


def _variant_ht():
    import hail as hl

    # 3 variants: v1/v2 missense in ENSG_A; v3 stop-gain (aaalt X) -> filtered out.
    # REVEL_score is a per-transcript dict; the two transcript values are equal (broadcast).
    rows = [
        {
            "locus": hl.locus("chr1", 100, "GRCh38"),
            "alleles": ["A", "C"],
            "aaref": "M",
            "aaalt": "T",
            "Ensembl_geneid": "ENSG_A;ENSG_A",
            "REVEL_score": {"t1": 0.9, "t2": 0.9},
            "MPC_score": {"t1": 2.0, "t2": 2.0},
            "CADD_phred": {"t1": 25.0, "t2": 25.0},
        },
        {
            "locus": hl.locus("chr1", 200, "GRCh38"),
            "alleles": ["G", "T"],
            "aaref": "R",
            "aaalt": "Q",
            "Ensembl_geneid": "ENSG_A;ENSG_B",
            "REVEL_score": {"t1": 0.1, "t2": 0.1},
            "MPC_score": {"t1": 0.5, "t2": 0.5},
            "CADD_phred": {"t1": 5.0, "t2": 5.0},
        },
        {
            "locus": hl.locus("chr1", 300, "GRCh38"),
            "alleles": ["C", "A"],
            "aaref": "W",
            "aaalt": "X",
            "Ensembl_geneid": "ENSG_A",
            "REVEL_score": {"t1": 0.99},
            "MPC_score": {"t1": 3.0},
            "CADD_phred": {"t1": 40.0},
        },
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            aaref=hl.tstr,
            aaalt=hl.tstr,
            Ensembl_geneid=hl.tstr,
            REVEL_score=hl.tdict(hl.tstr, hl.tfloat64),
            MPC_score=hl.tdict(hl.tstr, hl.tfloat64),
            CADD_phred=hl.tdict(hl.tstr, hl.tfloat64),
        ),
        key=["locus", "alleles"],
    )


def _agg_spec():
    from hvantk.algorithms.annotation.spec import AggregateSpec, ScoreSpec

    return AggregateSpec(
        by="Ensembl_geneid",
        to="gene_id",
        filter="missense",
        reduce="max",
        scores=(
            ScoreSpec("revel", "REVEL_score", ("mean", "max", "frac_gt_0.5")),
            ScoreSpec("mpc", "MPC_score", ("mean", "max")),
            ScoreSpec("cadd", "CADD_phred", ("mean", "max")),
        ),
    )


@pytest.mark.hail
def test_aggregate_to_gene_collapses_missense_variants_per_gene(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.transforms import aggregate_to_gene

    grouped = aggregate_to_gene(_variant_ht(), _agg_spec())
    assert list(grouped.key) == ["Ensembl_geneid"]
    assert set(grouped.row) == {
        "Ensembl_geneid",
        "revel_mean",
        "revel_max",
        "revel_frac_gt_0.5",
        "mpc_mean",
        "mpc_max",
        "cadd_mean",
        "cadd_max",
        "n_possible_missense",
    }
    d = {r.Ensembl_geneid: r for r in grouped.collect()}
    # ENSG_A: v1(0.9) and v2(0.1) are missense; v3(stop-gain) filtered out.
    assert d["ENSG_A"].n_possible_missense == 2
    assert d["ENSG_A"].revel_max == pytest.approx(0.9)
    assert d["ENSG_A"].revel_mean == pytest.approx(0.5)
    assert d["ENSG_A"]["revel_frac_gt_0.5"] == pytest.approx(0.5)  # 1 of 2 > 0.5
    # ENSG_B: only v2 (0.1).
    assert d["ENSG_B"].n_possible_missense == 1
    assert d["ENSG_B"].revel_max == pytest.approx(0.1)


@pytest.mark.hail
def test_row_count_column_is_named_by_count_name(hail_session):
    """A non-dbNSFP source must be able to name its own row-count column.

    The count is "rows that survived the filter", which is only 'possible missense'
    for dbNSFP. A PTM-site or eQTL-pair source counts something else entirely, and
    two such axes composed together would otherwise collide on one hardcoded name.
    """
    from dataclasses import replace

    from hvantk.algorithms.annotation.transforms import aggregate_to_gene

    grouped = aggregate_to_gene(
        _variant_ht(), replace(_agg_spec(), count_name="n_ptm_sites")
    )
    assert "n_ptm_sites" in grouped.row
    assert "n_possible_missense" not in grouped.row
    d = {r.Ensembl_geneid: r for r in grouped.collect()}
    assert d["ENSG_A"].n_ptm_sites == 2
    assert d["ENSG_B"].n_ptm_sites == 1


@pytest.mark.hail
def test_stop_gain_and_stop_loss_are_excluded(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.transforms import PREDICATES

    ht = _variant_ht()
    missense = ht.filter(PREDICATES["missense"](ht))
    # v3 (aaalt == "X") is the only non-missense row.
    assert missense.count() == 2


def _scalar_ht():
    """A row-level source whose score columns are already scalars, not transcript dicts.

    dbNSFP broadcasts each score across a variant's transcripts, so its columns are
    ``dict<transcript_id, float>``. Most other row-level sources are not shaped that way:
    a UniProt PTM site has one ``n_observations``, a GTEx eQTL pair one ``slope``. Two
    sites/pairs per group here so the group_by has something to reduce.
    """
    import hail as hl

    rows = [
        {"uniprot_id": "P00001", "n_observations": 4, "score": 0.2},
        {"uniprot_id": "P00001", "n_observations": 10, "score": 0.8},
        {"uniprot_id": "P00002", "n_observations": 1, "score": 0.5},
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(uniprot_id=hl.tstr, n_observations=hl.tint32, score=hl.tfloat64),
        key=["uniprot_id"],
    )


@pytest.mark.hail
def test_identity_reduce_supports_sources_with_scalar_score_columns(hail_session):
    """``reduce: identity`` lets a non-dbNSFP source through the aggregate path.

    The default ``max`` reducer calls ``.values()`` on each score column, which only
    exists on a dict. Without an identity reducer the whole aggregate path is unusable
    for any source that is not shaped like dbNSFP -- PTM density, GTEx eQTL and pQTL all
    hit this.
    """
    from hvantk.algorithms.annotation.spec import AggregateSpec, ScoreSpec
    from hvantk.algorithms.annotation.transforms import aggregate_to_gene

    spec = AggregateSpec(
        by="uniprot_id",
        to="uniprot_id",
        reduce="identity",
        count_name="n_ptm_sites",
        scores=(ScoreSpec("obs", "n_observations", ("mean", "max")),
                ScoreSpec("sc", "score", ("max",))),
    )
    grouped = aggregate_to_gene(_scalar_ht(), spec)
    assert set(grouped.row) == {"uniprot_id", "obs_mean", "obs_max", "sc_max",
                                "n_ptm_sites"}
    d = {r.uniprot_id: r for r in grouped.collect()}
    assert d["P00001"].n_ptm_sites == 2
    assert d["P00001"].obs_max == 10
    assert d["P00001"].obs_mean == pytest.approx(7.0)
    assert d["P00001"].sc_max == pytest.approx(0.8)
    assert d["P00002"].n_ptm_sites == 1


def test_identity_reducer_is_registered():
    """Fast (no-Hail) guard that the token exists, so a spec typo fails loudly."""
    from hvantk.algorithms.annotation.transforms import REDUCERS

    assert "identity" in REDUCERS
    sentinel = object()
    assert REDUCERS["identity"](sentinel) is sentinel
