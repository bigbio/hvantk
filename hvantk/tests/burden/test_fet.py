import math
import pandas as pd
from scipy.stats import fisher_exact

from hvantk.algorithms.burden.fet import (
    fisher_2x2,
    add_fisher,
    pick_min_p,
    finalize_reductions,
    apply_mtc,
    run_gene_burden_fet,
)


def test_fisher_2x2_matches_scipy():
    p, orr = fisher_2x2(8, 2, 100, 400)
    exp_or, exp_p = fisher_exact([[8, 2], [100, 400]], alternative="two-sided")
    assert math.isclose(p, exp_p, rel_tol=1e-12)
    assert math.isclose(orr, exp_or, rel_tol=1e-9)


def test_add_fisher_derives_cd():
    counts = pd.DataFrame(
        [{"gene": "G", "route": "lof", "a": 5, "b": 1, "n_case": 100, "n_control": 100}]
    )
    out = add_fisher(counts)
    assert out.loc[0, "c"] == 95  # n_case - a
    assert out.loc[0, "d"] == 99  # n_control - b
    exp_or, exp_p = fisher_exact([[5, 1], [95, 99]], alternative="two-sided")
    assert math.isclose(out.loc[0, "p"], exp_p, rel_tol=1e-12)
    assert math.isclose(out.loc[0, "odds_ratio"], exp_or, rel_tol=1e-9)


def test_pick_min_p_selects_winning_route():
    fdf = pd.DataFrame(
        [
            {"gene": "G", "route": "lof", "p": 0.2, "odds_ratio": 2.0},
            {"gene": "G", "route": "mis", "p": 0.01, "odds_ratio": 9.0},
        ]
    )
    out = pick_min_p(fdf).set_index("gene")
    assert out.loc["G", "route"] == "mis"
    assert math.isclose(out.loc["G", "minp"], 0.01)


def test_finalize_reductions_math():
    rdf = pd.DataFrame(
        [
            {
                "gene": "G",
                "route": "lof",
                "n_case_var": 4,
                "conc_num": 3,
                "conc_den": 6,
                "n_case_private": 3,
                "score_sum": 2.0,
                "score_n": 2,
                "drivers": [
                    {"cc": 3, "ctrl_freq": 0.004},
                    {"cc": 1, "ctrl_freq": 0.05},
                ],
            }
        ]
    )
    out = finalize_reductions(rdf).set_index(["gene", "route"])
    assert math.isclose(out.loc[("G", "lof"), "conc"], 0.5)
    assert math.isclose(
        out.loc[("G", "lof"), "driver_af"], 0.004
    )  # ctrl_freq of the max-cc driver
    assert math.isclose(out.loc[("G", "lof"), "frac_case_private"], 0.75)
    assert math.isclose(out.loc[("G", "lof"), "mean_score_case"], 1.0)


def test_apply_mtc_bh_monotone():
    df = pd.DataFrame({"minp": [0.001, 0.02, 0.5]})
    out = apply_mtc(df, "bh")
    assert list(out["p_adj"]) == sorted(out["p_adj"])  # BH-adjusted preserve order
    assert out["p_adj"].iloc[0] >= 0.001


def test_run_end_to_end_prior_and_columns():
    counts = pd.DataFrame(
        [
            {
                "gene": "G",
                "route": "lof",
                "a": 8,
                "b": 1,
                "n_case": 100,
                "n_control": 400,
            },
            {
                "gene": "G",
                "route": "mis",
                "a": 2,
                "b": 2,
                "n_case": 100,
                "n_control": 400,
            },
        ]
    )
    reductions = pd.DataFrame(
        [
            {
                "gene": "G",
                "route": "lof",
                "n_case_var": 5,
                "conc_num": 4,
                "conc_den": 8,
                "n_case_private": 5,
                "score_sum": 0.0,
                "score_n": 0,
                "drivers": [{"cc": 4, "ctrl_freq": 0.0}],
            },
            {
                "gene": "G",
                "route": "mis",
                "n_case_var": 2,
                "conc_num": 1,
                "conc_den": 2,
                "n_case_private": 1,
                "score_sum": 0.0,
                "score_n": 0,
                "drivers": [{"cc": 1, "ctrl_freq": 0.01}],
            },
        ]
    )
    out = run_gene_burden_fet(counts, reductions, mtc="bh").set_index("gene")
    assert out.loc["G", "route"] == "lof"  # lof is the min-p route
    for col in [
        "minp",
        "odds_ratio",
        "n_case_var",
        "conc",
        "driver_af",
        "mean_score_case",
        "frac_case_private",
        "p_adj",
    ]:
        assert col in out.columns
    assert out.loc["G", "n_case_var"] == 5  # reductions taken from the WINNING route


def test_apply_mtc_bh_unsorted_and_monotone_correction():
    # raw p unsorted; sorted=[0.01,0.03,0.04], q=n*p/rank=[0.03,0.045,0.04];
    # BH step-up (min from largest) -> sorted adj=[0.03,0.04,0.04];
    # mapped back to original order [0.04,0.01,0.03] -> [0.04,0.03,0.04]
    df = pd.DataFrame({"minp": [0.04, 0.01, 0.03]})
    out = apply_mtc(df, "bh")
    got = [round(x, 6) for x in out["p_adj"]]
    assert got == [0.04, 0.03, 0.04]


def test_run_default_no_mtc_keeps_raw_prior():
    counts = pd.DataFrame(
        [
            {
                "gene": "G",
                "route": "lof",
                "a": 8,
                "b": 1,
                "n_case": 100,
                "n_control": 400,
            },
            {
                "gene": "G",
                "route": "mis",
                "a": 2,
                "b": 2,
                "n_case": 100,
                "n_control": 400,
            },
        ]
    )
    reductions = pd.DataFrame(
        [
            {
                "gene": "G",
                "route": "lof",
                "n_case_var": 5,
                "conc_num": 4,
                "conc_den": 8,
                "n_case_private": 5,
                "score_sum": 0.0,
                "score_n": 0,
                "drivers": [{"cc": 4, "ctrl_freq": 0.0}],
            },
            {
                "gene": "G",
                "route": "mis",
                "n_case_var": 2,
                "conc_num": 1,
                "conc_den": 2,
                "n_case_private": 1,
                "score_sum": 0.0,
                "score_n": 0,
                "drivers": [{"cc": 1, "ctrl_freq": 0.01}],
            },
        ]
    )
    out = run_gene_burden_fet(counts, reductions)  # mtc defaults to None
    assert "p_adj" not in out.columns
    # winning route is lof; its raw min-p equals the standalone fisher p for that 2x2
    exp_p, _ = fisher_2x2(8, 1, 92, 399)
    assert math.isclose(out.set_index("gene").loc["G", "minp"], exp_p, rel_tol=1e-9)
