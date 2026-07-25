"""Fisher-exact gene-burden statistics on the small collected frames (no Hail).

Given per-(gene,route) 2x2 counts and reduction inputs from ``aggregate.py``, compute
the Fisher p/OR, take the min-p route as the gene's prior, attach the architecture
reductions of that winning route, and (optionally) a multiple-testing-corrected column.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return float(p_value), float(odds_ratio)


def add_fisher(counts_df: pd.DataFrame) -> pd.DataFrame:
    df = counts_df.copy()
    df["c"] = df["n_case"] - df["a"]
    df["d"] = df["n_control"] - df["b"]
    # itertuples avoids the per-row Series construction of apply(axis=1); the
    # per-row scipy fisher_exact call is unchanged, as is row order.
    stats = [
        fisher_2x2(int(r.a), int(r.b), int(r.c), int(r.d))
        for r in df.itertuples(index=False)
    ]
    df["p"] = [s[0] for s in stats]
    df["odds_ratio"] = [s[1] for s in stats]
    return df


def pick_min_p(fisher_df: pd.DataFrame) -> pd.DataFrame:
    df = fisher_df.sort_values(["gene", "p", "route"]).drop_duplicates(
        "gene", keep="first"
    )
    return df.rename(columns={"p": "minp"})[["gene", "route", "minp", "odds_ratio"]]


def _driver_af(drivers) -> float:
    # "af" here is the control-CARRIER frequency (carriers / control samples) of the
    # max-case-carrier ("cc") driver variant, not an allele frequency; the exact
    # allele-vs-carrier semantics are pinned by the CHD reproduction gate -- do not
    # change without re-running it.
    if not isinstance(drivers, (list, np.ndarray)) or len(drivers) == 0:
        return float("nan")
    top = max(drivers, key=lambda d: d["cc"])
    return float(top["ctrl_freq"])


def finalize_reductions(reductions_df: pd.DataFrame) -> pd.DataFrame:
    df = reductions_df.copy()
    df["conc"] = np.where(df["conc_den"] > 0, df["conc_num"] / df["conc_den"], np.nan)
    df["driver_af"] = df["drivers"].apply(_driver_af)
    df["frac_case_private"] = np.where(
        df["n_case_var"] > 0, df["n_case_private"] / df["n_case_var"], np.nan
    )
    df["mean_score_case"] = np.where(
        df["score_n"] > 0, df["score_sum"] / df["score_n"], np.nan
    )
    return df[
        [
            "gene",
            "route",
            "n_case_var",
            "conc",
            "driver_af",
            "mean_score_case",
            "frac_case_private",
        ]
    ]


def apply_mtc(df: pd.DataFrame, method: str) -> pd.DataFrame:
    out = df.copy()
    p = out["minp"].to_numpy(dtype=float)
    n = len(p)
    if method == "bonferroni":
        out["p_adj"] = np.minimum(p * n, 1.0)
    elif method == "bh":
        order = np.argsort(p)
        ranked = p[order] * n / (np.arange(1, n + 1))
        # enforce monotonicity from the largest p downward
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adj = np.empty(n)
        adj[order] = np.minimum(ranked, 1.0)
        out["p_adj"] = adj
    else:
        raise ValueError(f"unknown MTC method {method!r}; use 'bh' or 'bonferroni'")
    return out


def run_gene_burden_fet(
    counts_df: pd.DataFrame, reductions_df: pd.DataFrame, *, mtc: str | None = None
) -> pd.DataFrame:
    winners = pick_min_p(add_fisher(counts_df))
    reds = finalize_reductions(reductions_df)
    out = winners.merge(reds, on=["gene", "route"], how="left")
    if mtc is not None:
        out = apply_mtc(out, mtc)
    return out
