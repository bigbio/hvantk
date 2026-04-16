"""Phase-2 PTM constraint tests (LMM and binned-interaction LMM).

Factored from:

- notebook M (per-tissue constraint LMM on GTEx) -> :func:`run_lmm`
- notebook K (per-cell-type binned-interaction LMM on Velmeshev cortex) ->
  :func:`run_binned_interaction_lmm`

Both functions are per-stratum: the caller iterates over tissues / cell types
and is responsible for assembling results. This matches the notebook style
and keeps the API surface small.

Example
-------
    >>> from hvantk.ptm.test import run_lmm, run_binned_interaction_lmm
    >>> result = run_lmm(df_heart, stratum="Heart")
    >>> print(result.beta_ptm, result.p_ptm)
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from hvantk.ptm.constants import (
    LMM_BINNED_MIN_CELL_N,
    LMM_BINNED_MIN_POS_EXPR,
    LMM_MIN_MIXED_GENES,
    LMM_MIN_N_NONPTM,
    LMM_MIN_N_PTM,
    LOG_AF_EPSILON,
)

logger = logging.getLogger(__name__)


@dataclass
class LMMResult:
    """Per-stratum result of the constraint LMM (notebook M Cell 5).

    Fields are populated on fit success; skipped strata return NaN betas and
    an explanatory ``note`` string (never raise).
    """

    stratum: str
    n_variants: int
    n_ptm: int
    n_nonptm: int
    n_genes: int
    n_mixed_genes: int
    beta_ptm: float
    se_ptm: float
    p_ptm: float
    converged: bool
    note: str


@dataclass
class BinnedLMMResult:
    """Per-stratum result of the binned-interaction LMM (notebook K Cell 4d).

    ``bin_levels`` always starts with ``"b0_none"`` and ends with the last
    quantile bin actually realized by ``pd.qcut(..., duplicates='drop')``.
    ``bin_betas`` / ``bin_ses`` / ``bin_pvalues`` are keyed by bin label and
    may be empty dicts when the fit was skipped.
    """

    stratum: str
    bin_levels: List[str]
    bin_betas: Dict[str, float] = field(default_factory=dict)
    bin_ses: Dict[str, float] = field(default_factory=dict)
    bin_pvalues: Dict[str, float] = field(default_factory=dict)
    n_variants: int = 0
    n_genes: int = 0
    converged: bool = False
    note: str = ""


def run_lmm(
    df: pd.DataFrame,
    stratum: str,
    gene_col: str = "gene",
    af_col: str = "af_filled",
    is_ptm_col: str = "is_ptm",
    eps: float = LOG_AF_EPSILON,
    min_n_ptm: int = LMM_MIN_N_PTM,
    min_n_nonptm: int = LMM_MIN_N_NONPTM,
    min_mixed_genes: int = LMM_MIN_MIXED_GENES,
) -> LMMResult:
    """Per-stratum constraint LMM: ``log_af ~ is_ptm + (1|gene)``.

    Replicates notebook_m Cell 5 mixedlm usage exactly. Filters:

    - ``af > 0`` on ``af_col``;
    - ``n_ptm >= min_n_ptm`` and ``n_nonptm >= min_n_nonptm``;
    - ``n_mixed_genes >= min_mixed_genes`` where a mixed gene has both PTM and
      non-PTM variants.

    Returns an :class:`LMMResult` with NaN betas and an explanatory ``note``
    when the filters fail. Never raises.

    Parameters
    ----------
    df : pandas.DataFrame
        Variants already filtered to the target stratum (caller slices).
    stratum : str
        Label for this stratum (e.g. tissue name); carried through to the
        result object for plotting/report.
    gene_col, af_col, is_ptm_col : str
        Column names.
    eps : float
        Pseudocount for ``log10(af + eps)``. Default: :data:`LOG_AF_EPSILON`.
    min_n_ptm, min_n_nonptm, min_mixed_genes : int
        Filter thresholds (notebook M defaults).
    """
    # Drop rows missing any required column (match notebook M's dropna).
    sub = df.dropna(subset=[af_col, gene_col, is_ptm_col]).copy()
    sub = sub[sub[af_col] > 0].copy()
    sub["log_af"] = np.log10(sub[af_col] + eps)
    sub["is_ptm"] = sub[is_ptm_col].astype(int)

    n_variants = int(len(sub))
    n_ptm = int(sub["is_ptm"].sum())
    n_nonptm = n_variants - n_ptm
    n_genes = int(sub[gene_col].nunique()) if n_variants else 0
    if n_variants:
        mix = sub.groupby(gene_col)["is_ptm"].agg(["min", "max"])
        n_mixed = int((mix["max"] > mix["min"]).sum())
    else:
        n_mixed = 0

    def _skipped(note: str) -> LMMResult:
        return LMMResult(
            stratum=stratum,
            n_variants=n_variants,
            n_ptm=n_ptm,
            n_nonptm=n_nonptm,
            n_genes=n_genes,
            n_mixed_genes=n_mixed,
            beta_ptm=float("nan"),
            se_ptm=float("nan"),
            p_ptm=float("nan"),
            converged=False,
            note=note,
        )

    if n_ptm < min_n_ptm:
        return _skipped(f"skipped: n_ptm={n_ptm} < {min_n_ptm}")
    if n_nonptm < min_n_nonptm:
        return _skipped(f"skipped: n_nonptm={n_nonptm} < {min_n_nonptm}")
    if n_mixed < min_mixed_genes:
        return _skipped(f"skipped: n_mixed_genes={n_mixed} < {min_mixed_genes}")

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm("log_af ~ is_ptm", data=sub, groups=sub[gene_col])
            mf = md.fit(method="lbfgs", reml=True)
    except Exception as e:  # pragma: no cover - fit errors bubbled into note
        return _skipped(f"fit failed: {str(e)[:80]}")

    return LMMResult(
        stratum=stratum,
        n_variants=n_variants,
        n_ptm=n_ptm,
        n_nonptm=n_nonptm,
        n_genes=n_genes,
        n_mixed_genes=n_mixed,
        beta_ptm=float(mf.params["is_ptm"]),
        se_ptm=float(mf.bse["is_ptm"]),
        p_ptm=float(mf.pvalues["is_ptm"]),
        converged=bool(mf.converged),
        note="",
    )


def run_binned_interaction_lmm(
    df: pd.DataFrame,
    expr_series: pd.Series,
    stratum: str,
    gene_col: str = "gene",
    af_col: str = "af_filled",
    is_ptm_col: str = "is_ptm",
    eps: float = LOG_AF_EPSILON,
    n_quantiles: int = 4,
    min_pos_expr: int = LMM_BINNED_MIN_POS_EXPR,
    min_cell_n: int = LMM_BINNED_MIN_CELL_N,
) -> BinnedLMMResult:
    """Binned-interaction LMM: ``log_af ~ is_ptm * C(expr_bin) + (1|gene)``.

    Replicates notebook_k Cell 4d mixedlm usage exactly:

    - ``expr = expr_series[gene]``, ``expr_log = log2(expr + 1)``;
    - ``pd.qcut(expr_log[expr > 0], q=n_quantiles, duplicates='drop')`` yields
      dynamic bin labels ``b1_Q1, b2_Q2, ..., bK_Qk`` where ``K <= n_quantiles``;
    - zero-expression variants are assigned to reference bin ``"b0_none"``;
    - ``statsmodels`` fits with Treatment coding (reference = first category);
    - per-bin beta is extracted via sum rule
      ``beta_b = beta_ref + beta_(is_ptm:expr_bin[T.b])``;
    - per-bin SE is ``sqrt(Var(a) + Var(b) + 2*Cov(a,b))`` via ``mf.cov_params()``;
    - the reference bin p-value is the main ``is_ptm`` term; non-reference
      bins report the interaction-term p.

    Parameters
    ----------
    df : pandas.DataFrame
        Variants already filtered to the stratum's "base" set.
    expr_series : pandas.Series
        Gene-to-expression mapping (index = gene symbol, values = expression
        for this stratum). Missing genes are treated as zero.
    stratum : str
        Label for this stratum.
    eps : float
        Pseudocount for ``log10(af + eps)``.
    n_quantiles : int
        Upper bound on quantile bins (``pd.qcut`` with ``duplicates='drop'``).
    min_pos_expr : int
        Minimum count of positive-expression variants before attempting qcut.
    min_cell_n : int
        Minimum count per ``(expr_bin, is_ptm)`` cell to proceed with the fit.
    """
    # Drop rows missing AF/gene/is_ptm, then af > 0 (matches notebook K).
    dfx = df.dropna(subset=[af_col, gene_col, is_ptm_col]).copy()
    dfx = dfx[dfx[af_col] > 0].copy()
    dfx["log_af"] = np.log10(dfx[af_col] + eps)
    dfx["is_ptm"] = dfx[is_ptm_col].astype(int)

    # Map gene -> expression; missing genes default to 0.
    expr_map = dict(expr_series)
    dfx["expr"] = dfx[gene_col].map(expr_map).fillna(0.0)
    dfx["expr_log"] = np.log2(dfx["expr"] + 1)

    not_expr = dfx["expr"] == 0
    pos = dfx.loc[~not_expr, "expr_log"]

    def _skipped(note: str, bin_levels=None) -> BinnedLMMResult:
        return BinnedLMMResult(
            stratum=stratum,
            bin_levels=list(bin_levels or ["b0_none"]),
            n_variants=int(len(dfx)),
            n_genes=int(dfx[gene_col].nunique()) if len(dfx) else 0,
            converged=False,
            note=note,
        )

    if len(pos) < min_pos_expr:
        return _skipped(f"skipped: only {len(pos)} positive-expr variants < {min_pos_expr}")

    try:
        pos_bins = pd.qcut(pos, q=n_quantiles, duplicates="drop")
    except Exception as e:
        return _skipped(f"qcut failed: {str(e)[:80]}")

    n_bins_realized = pos_bins.cat.categories.size
    if n_bins_realized == 0:
        return _skipped("skipped: no quantile bins realized")

    bin_labels = [f"b{i + 1}_Q{i + 1}" for i in range(n_bins_realized)]
    bin_map = dict(zip(pos_bins.cat.categories, bin_labels))
    assigned = pd.Series("b0_none", index=dfx.index, dtype=object)
    assigned.loc[~not_expr] = pos_bins.map(bin_map)
    # Any residual NaN (shouldn't happen) falls back to the first positive bin
    # to match notebook K's defensive handling.
    assigned = assigned.fillna(bin_labels[0] if bin_labels else "b0_none")
    ordered = ["b0_none"] + bin_labels
    dfx["expr_bin"] = pd.Categorical(assigned, categories=ordered, ordered=True)

    # Sparsity gate: every (expr_bin, is_ptm) cell must have >= min_cell_n.
    bc = dfx.groupby(["expr_bin", "is_ptm"], observed=False).size().unstack(fill_value=0)
    if (bc < min_cell_n).any().any():
        return _skipped("skipped: sparse bins", bin_levels=ordered)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm("log_af ~ is_ptm * expr_bin", data=dfx, groups=dfx[gene_col])
            mf = md.fit(method="lbfgs", reml=True)
    except Exception as e:
        return _skipped(f"fit failed: {str(e)[:80]}", bin_levels=ordered)

    beta_ref = float(mf.params["is_ptm"])
    bin_betas: Dict[str, float] = {"b0_none": beta_ref}
    bin_ses: Dict[str, float] = {"b0_none": float(mf.bse["is_ptm"])}
    bin_pvalues: Dict[str, float] = {"b0_none": float(mf.pvalues["is_ptm"])}

    try:
        vc = mf.cov_params()
    except Exception:
        vc = None

    for b in bin_labels:
        int_name = f"is_ptm:expr_bin[T.{b}]"
        if int_name not in mf.params.index:
            continue
        int_beta = float(mf.params[int_name])
        bin_betas[b] = beta_ref + int_beta
        bin_pvalues[b] = float(mf.pvalues[int_name])
        if vc is not None and int_name in vc.index:
            v_a = vc.loc["is_ptm", "is_ptm"]
            v_b = vc.loc[int_name, int_name]
            c_ab = vc.loc["is_ptm", int_name]
            bin_ses[b] = float(np.sqrt(v_a + v_b + 2 * c_ab))
        else:
            bin_ses[b] = float(mf.bse.get(int_name, float("nan")))

    return BinnedLMMResult(
        stratum=stratum,
        bin_levels=ordered,
        bin_betas=bin_betas,
        bin_ses=bin_ses,
        bin_pvalues=bin_pvalues,
        n_variants=int(len(dfx)),
        n_genes=int(dfx[gene_col].nunique()),
        converged=bool(mf.converged),
        note="",
    )


__all__ = [
    "LMMResult",
    "BinnedLMMResult",
    "run_lmm",
    "run_binned_interaction_lmm",
]
