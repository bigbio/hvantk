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
def test_stop_gain_and_stop_loss_are_excluded(hail_session):
    import hail as hl

    from hvantk.algorithms.annotation.transforms import PREDICATES

    ht = _variant_ht()
    missense = ht.filter(PREDICATES["missense"](ht))
    # v3 (aaalt == "X") is the only non-missense row.
    assert missense.count() == 2
