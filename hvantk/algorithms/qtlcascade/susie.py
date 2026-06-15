"""Pure-NumPy SuSiE-RSS + ``coloc.susie`` (issue #193 — drop the R bridge).

Faithful port of the two R routines the fine-mapping confirmation layer used to
shell out to:

- :func:`susie_rss` ports ``susieR::susie_rss(z, R, n, L)`` along its
  *n-provided* path: regularise the z-scores, reframe as sufficient statistics
  (``XtX = (n-1) R``, ``Xty = sqrt(n-1) z``, ``var_y = 1``) and fit SuSiE with
  the IBSS algorithm holding the residual variance fixed
  (``estimate_residual_variance = FALSE``, the RSS default).
- :func:`coloc_susie` ports ``coloc::coloc.susie``: a pairwise colocalization
  over every (trait-1 single-effect, trait-2 single-effect) credible-set pair,
  using each single effect's per-variant log Bayes factor vector. It reuses the
  validated coloc Bayes-factor algebra in :mod:`hvantk.algorithms.qtlcascade.coloc`.

The Python SuSiE is *not* bit-identical to ``susieR`` (the documented caveat for
any reimplementation), but on the project's positive/negative controls it
reproduces the R credible-set counts exactly and ``coloc.susie`` PP4 to within
~0.005, giving the same CONFIRMED/REFUTED verdicts (validated AF→MYOZ1 CONFIRMED,
CHD 17q21/NSF REFUTED).

References
----------
- Wang et al. (2020) JRSS-B 82(5):1273 — SuSiE (IBSS).
- Zou, Carbonetto, Wang, Stephens (2022) PLoS Genet 18(7):e1010299 — SuSiE-RSS.
- Wallace (2021) PLoS Genet 17(9):e1009440 — coloc with SuSiE (``coloc.susie``).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np

from hvantk.algorithms.qtlcascade.coloc import _logsumexp, _logdiff
from hvantk.algorithms.qtlcascade.constants import (
    DEFAULT_SUSIE_L,
    DEFAULT_SUSIE_COVERAGE,
    DEFAULT_SUSIE_MIN_ABS_CORR,
    DEFAULT_SUSIE_MAX_ITER,
    DEFAULT_COLOC_SUSIE_P1,
    DEFAULT_COLOC_SUSIE_P2,
    DEFAULT_COLOC_SUSIE_P12,
)

# Variants subsampled when computing credible-set purity (matches susieR's cap).
_PURITY_MAX_VARIANTS = 100


# ---------------------------------------------------------------------------
# Single-effect-regression prior-variance MLE (susieR estimate_prior_method="optim")
# ---------------------------------------------------------------------------


def _ser_loglik(V: float, z2: np.ndarray, shat2: np.ndarray) -> float:
    """Log marginal likelihood of a single-effect regression as a function of ``V``.

    ``lbf_j(V) = 0.5 log(shat2/(shat2+V)) + 0.5 z2_j V/(V+shat2)`` with a uniform
    prior ``pi = 1/p``. Returns ``logsumexp_j(lbf_j) - log(p)`` (the ``-log(p)``
    is constant in ``V`` and irrelevant to the argmax, but kept for honesty).
    """
    p = len(z2)
    if V <= 0:
        return float(_logsumexp(np.zeros(p)) - np.log(p))
    lbf = 0.5 * np.log(shat2 / (shat2 + V)) + 0.5 * z2 * (V / (V + shat2))
    return float(_logsumexp(lbf) - np.log(p))


def _optimize_V(z2: np.ndarray, shat2: np.ndarray) -> float:
    """Maximise the SER log-likelihood over ``V >= 0`` (grid + golden-section).

    Returns ``0.0`` when ``V = 0`` is not strictly improved (the effect is turned
    off), mirroring ``susieR``'s behaviour. NumPy-only (no scipy.optimize).
    """
    grid = np.concatenate([[0.0], np.logspace(-7.0, 0.5, 60)])
    ll = np.array([_ser_loglik(V, z2, shat2) for V in grid])
    i = int(np.argmax(ll))
    lo = grid[max(0, i - 1)]
    hi = grid[min(len(grid) - 1, i + 1)]
    if hi > lo:
        gr = (np.sqrt(5.0) - 1.0) / 2.0
        a, b = lo, hi
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        for _ in range(60):
            if _ser_loglik(c, z2, shat2) < _ser_loglik(d, z2, shat2):
                a = c
            else:
                b = d
            c = b - gr * (b - a)
            d = a + gr * (b - a)
            if abs(b - a) < 1e-9:
                break
        Vhat = 0.5 * (a + b)
    else:
        Vhat = grid[i]
    if _ser_loglik(Vhat, z2, shat2) <= _ser_loglik(0.0, z2, shat2) + 1e-12:
        return 0.0
    return float(Vhat)


# ---------------------------------------------------------------------------
# SuSiE-RSS fit
# ---------------------------------------------------------------------------


@dataclass
class SusieFit:
    """Result of :func:`susie_rss`.

    Attributes
    ----------
    alpha : ndarray (L, p)
        Posterior inclusion probability of each variant in each single effect.
    lbf_variable : ndarray (L, p)
        Per-effect, per-variant log Bayes factor (``susieR`` ``lbf_variable``);
        consumed by :func:`coloc_susie`.
    V : ndarray (L,)
        Estimated prior variance per effect (``0`` => effect off).
    cs : list[ndarray]
        Index arrays of the retained (purity-passing) 95% credible sets.
    cs_effect : list[int]
        The effect index backing each retained credible set.
    """

    alpha: np.ndarray
    lbf_variable: np.ndarray
    V: np.ndarray
    cs: List[np.ndarray] = field(default_factory=list)
    cs_effect: List[int] = field(default_factory=list)


def susie_rss(
    z: np.ndarray,
    R: np.ndarray,
    n: int,
    L: int = DEFAULT_SUSIE_L,
    coverage: float = DEFAULT_SUSIE_COVERAGE,
    min_abs_corr: float = DEFAULT_SUSIE_MIN_ABS_CORR,
    max_iter: int = DEFAULT_SUSIE_MAX_ITER,
    tol: float = 1e-3,
    scaled_prior_variance: float = 0.2,
) -> SusieFit:
    """Fit SuSiE-RSS from z-scores and a reference-LD matrix.

    Parameters
    ----------
    z : ndarray (p,)
        Per-variant z-scores (``beta / se``). Non-finite entries are zeroed.
    R : ndarray (p, p)
        Reference-LD correlation matrix (same variant order as ``z``).
    n : int
        Trait sample size.
    L : int
        Maximum number of single effects.
    coverage, min_abs_corr : float
        Credible-set coverage and purity (min absolute pairwise ``r``) threshold.
    max_iter, tol : int, float
        IBSS iteration cap and convergence tolerance on the total posterior mean.
    scaled_prior_variance : float
        Initial prior variance (fraction of ``var_y = 1``); re-estimated per effect.

    Returns
    -------
    SusieFit
    """
    z = np.asarray(z, dtype=float).copy()
    z[~np.isfinite(z)] = 0.0
    R = np.asarray(R, dtype=float)
    p = len(z)
    if p == 0:
        return SusieFit(alpha=np.zeros((0, 0)), lbf_variable=np.zeros((0, 0)),
                        V=np.zeros(0))

    # susie_rss z-score regularisation (n-provided path).
    adj = (n - 1) / (z ** 2 + n - 2)
    z = np.sqrt(adj) * z

    # Sufficient statistics for standardised X and y (var_y = 1).
    XtX = (n - 1) * R
    Xty = np.sqrt(n - 1) * z
    sigma2 = 1.0                       # residual variance, held fixed (RSS default)
    dj = np.diag(XtX).astype(float).copy()
    dj[dj <= 0] = n - 1                # guard a degenerate reference diagonal
    shat2 = sigma2 / dj

    L = min(L, p)
    alpha = np.zeros((L, p))
    mu = np.zeros((L, p))
    lbf_var = np.zeros((L, p))
    V = np.full(L, float(scaled_prior_variance))
    b = np.zeros((L, p))               # alpha*mu contribution per effect
    b_tot = np.zeros(p)
    Xtb = XtX @ b_tot
    logpi = -np.log(p)

    for _ in range(max_iter):
        b_tot_old = b_tot.copy()
        for eff in range(L):
            # Residualise: subtract every effect except `eff`.
            Xtr = Xty - (Xtb - XtX @ b[eff])
            bhat = Xtr / dj
            z2 = (bhat ** 2) / shat2
            Veff = _optimize_V(z2, shat2)
            V[eff] = Veff
            if Veff <= 0:
                alpha[eff] = 1.0 / p
                mu[eff] = 0.0
                lbf_var[eff] = 0.0
                b_new = np.zeros(p)
            else:
                lbf = (0.5 * np.log(shat2 / (shat2 + Veff))
                       + 0.5 * z2 * (Veff / (Veff + shat2)))
                lbf_var[eff] = lbf
                w = lbf + logpi
                w -= w.max()
                a = np.exp(w)
                a /= a.sum()
                alpha[eff] = a
                post_var = (shat2 * Veff) / (shat2 + Veff)
                mu[eff] = post_var * (Xtr / sigma2)
                b_new = a * mu[eff]
            Xtb += XtX @ (b_new - b[eff])
            b_tot += b_new - b[eff]
            b[eff] = b_new
        if np.max(np.abs(b_tot - b_tot_old)) < tol * max(1.0, np.max(np.abs(b_tot))):
            break

    fit = SusieFit(alpha=alpha, lbf_variable=lbf_var, V=V)
    _annotate_credible_sets(fit, R, coverage, min_abs_corr)
    return fit


def _annotate_credible_sets(fit: SusieFit, R: np.ndarray,
                            coverage: float, min_abs_corr: float) -> None:
    """``susie_get_cs``: per-effect 95% CS with a purity filter; drop duplicates."""
    seen: List[frozenset] = []
    rng = np.random.RandomState(1)     # deterministic purity subsampling
    for eff in range(fit.alpha.shape[0]):
        if fit.V[eff] <= 0:
            continue
        a = fit.alpha[eff]
        order = np.argsort(a)[::-1]
        csum = np.cumsum(a[order])
        k = int(np.searchsorted(csum, coverage) + 1)
        idx = np.sort(order[:k])
        sub = idx
        if len(sub) > _PURITY_MAX_VARIANTS:
            sub = np.sort(rng.choice(idx, _PURITY_MAX_VARIANTS, replace=False))
        if len(sub) == 1:
            purity = 1.0
        else:
            sub_r = np.abs(R[np.ix_(sub, sub)])
            purity = float(sub_r[~np.eye(len(sub), dtype=bool)].min())
        if purity < min_abs_corr:
            continue
        key = frozenset(int(i) for i in idx)
        if key in seen:
            continue
        seen.append(key)
        fit.cs.append(idx)
        fit.cs_effect.append(eff)


# ---------------------------------------------------------------------------
# coloc.susie
# ---------------------------------------------------------------------------


def _coloc_from_lbf(lbf1: np.ndarray, lbf2: np.ndarray,
                    p1: float, p2: float, p12: float) -> float:
    """PP.H4 for one pair of single-effect log-BF vectors (``coloc.bf_bf`` algebra).

    Same H0–H4 combination as :func:`coloc.coloc_abf`, but the per-variant log
    Bayes factors come straight from the two single effects rather than from
    ``beta/se``.
    """
    s1 = _logsumexp(lbf1)
    s2 = _logsumexp(lbf2)
    s_both = _logsumexp(lbf1 + lbf2)
    log_h = np.array([
        0.0,
        np.log(p1) + s1,
        np.log(p2) + s2,
        np.log(p1) + np.log(p2) + _logdiff(s1 + s2, s_both),
        np.log(p12) + s_both,
    ])
    post = np.exp(log_h - log_h.max())
    post /= post.sum()
    return float(post[4])


def coloc_susie(
    fit1: SusieFit,
    fit2: SusieFit,
    p1: float = DEFAULT_COLOC_SUSIE_P1,
    p2: float = DEFAULT_COLOC_SUSIE_P2,
    p12: float = DEFAULT_COLOC_SUSIE_P12,
) -> float:
    """Maximum PP.H4 over all credible-set pairs of two SuSiE fits.

    Returns ``0.0`` if either trait has no credible set (matching the R script's
    ``csg == 0 || cse == 0`` guard).
    """
    if not fit1.cs_effect or not fit2.cs_effect:
        return 0.0
    best = 0.0
    for i in fit1.cs_effect:
        for j in fit2.cs_effect:
            best = max(best, _coloc_from_lbf(
                fit1.lbf_variable[i], fit2.lbf_variable[j], p1, p2, p12))
    return best
