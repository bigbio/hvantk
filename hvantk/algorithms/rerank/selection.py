"""Within-axis feature selection, run on a TRAINING SLICE ONLY.

Three steps, in the order benchmarked by Perez-Riverol et al. (PLOS ONE 2017,
doi:10.1371/journal.pone.0189875): univariate filter -> redundancy filter -> wrapper RFE.
Their Table 1 is the reason for that order: prefiltering cut RFE runtime ~3x (35 min -> 11)
at equal accuracy, and inner CV3/CV7/CV10 all gave the same RMSE, so 3 inner folds suffice.

EVERY function here takes a training slice. Nothing in this module may see held-out rows.
Selection that has seen the test labels inflates the reported metric (Ambroise & McLachlan
2002), and here the reported delta-AUC IS the scientific claim. `test_selection_nesting.py`
is the guardrail that proves the caller respects this.

Departure from the reference workflow: it used a fixed correlation cutoff, appropriate to a
continuous outcome. With a binary label and a few hundred positives, a fixed |AUC-0.5|
threshold is roughly a 1-SE gate -- for 219 positives against 19,229 negatives the null SE
is about 0.020 -- so noise passes. The univariate step uses within-axis BH-FDR instead,
which adapts to both cohort size and axis width with no per-cohort tuning.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class UnivariateStat:
    """Per-column univariate screening result.

    For untestable columns (all-NaN, or one arm empty) ``p`` defaults to ``1.0``
    while ``z`` is ``nan`` -- they are not interchangeable stand-ins for "no result".
    Consumers must gate on ``.passed``, never on ``np.isfinite(p)`` (always true,
    even for untestable columns) or on ``p`` alone.
    """

    auc: float
    z: float
    p: float
    passed: bool
    n_pos: int
    n_neg: int


def univariate_auc(x, y):
    """Rank-based AUC of a single feature, using only rows where it is defined.

    Returns ``(auc, n_pos_eff, n_neg_eff)``. Directionless use is intended: an
    anti-predictive feature (AUC near 0) is informative, and the caller scores on
    ``|auc - 0.5|``. Rank-based so that rankscore, phred and raw scales compare fairly.
    """
    from scipy.stats import rankdata

    x = np.asarray(x, dtype=float)
    y = np.asarray(y)
    ok = ~np.isnan(x)
    xs, ys = x[ok], y[ok]
    n_pos = int((ys == 1).sum())
    n_neg = int((ys == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), n_pos, n_neg
    ranks = rankdata(xs)
    auc = (ranks[ys == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return float(auc), n_pos, n_neg


def _null_se(n_pos: int, n_neg: int) -> float:
    """SE of AUC under the null. Bamber's bound; adapts the gate to cohort size."""
    return float(np.sqrt((1.0 / 12.0) * (1.0 / n_pos + 1.0 / n_neg)))


def univariate_filter(X, y, columns, q: float = 0.10):
    """Score each column by |AUC-0.5| and keep those surviving BH-FDR at ``q``.

    FDR is computed WITHIN the axis, so a 50-column dbNSFP axis is automatically held to a
    stricter bar than a 1-column expression axis. Columns that are entirely missing, or
    whose effective arm is empty, are reported as not passing rather than raising.
    """
    from scipy.stats import norm

    stats = {}
    for col in columns:
        auc, n_pos, n_neg = univariate_auc(X[col].to_numpy(), y)
        if np.isnan(auc) or n_pos == 0 or n_neg == 0:
            stats[col] = UnivariateStat(auc, float("nan"), 1.0, False, n_pos, n_neg)
            continue
        z = (auc - 0.5) / _null_se(n_pos, n_neg)
        p = 2.0 * norm.sf(abs(z))
        stats[col] = UnivariateStat(auc, float(z), float(p), False, n_pos, n_neg)

    testable = [c for c in columns if np.isfinite(stats[c].z)]
    if not testable:
        return stats
    order = sorted(testable, key=lambda c: stats[c].p)
    m = len(order)
    passing = set()
    for rank, col in enumerate(order, start=1):
        if stats[col].p <= q * rank / m:
            passing.update(order[:rank])   # BH: everything up to the largest passing rank
    for col in passing:
        s = stats[col]
        stats[col] = UnivariateStat(s.auc, s.z, s.p, True, s.n_pos, s.n_neg)
    return stats
