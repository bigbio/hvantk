"""Headline metrics and the per-axis ablation table for one rerank run.

Three numbers describe the model as a whole (ROC AUC, average precision, Brier on
min-max-rescaled scores, plus a quantile-binned calibration curve), and then one row per
feature axis describes what that axis ADDS: the out-of-fold AUC of baseline+axis, with a
paired bootstrap interval on the difference from the baseline alone.

``_raw_oof`` is uncalibrated on purpose. The ablation compares discrimination between two
feature sets, and isotonic calibration is a monotone transform -- it cannot change AUC, and
fitting it 5x per fold would multiply the cost of the table by the calibration folds for no
change in the answer. ``ReRanker.score`` still calibrates, because the SCORES it produces
are read as probabilities; these deltas are not.

The per-axis interval here resamples GENES only. It is one of at least three variance
components; the CV partition and the best-of-N axis selection are the other two, and they
live in ``nulls.py`` and in ``seed_sweep`` below.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import cross_val_predict

from hvantk.algorithms.rerank.reranker import _cv, _gbm
from hvantk.algorithms.rerank.seeds import DEFAULT_SEED


ABLATION_FOLDS = 5
"""The ablation's fold count. ``Config.folds`` governs ``ReRanker.score`` only."""


def _raw_oof(
    matrix, cols, y, selector=None, groups=None, seed=DEFAULT_SEED, folds=ABLATION_FOLDS
):
    """Uncalibrated out-of-fold probabilities for one feature subset.

    Mirrors ReRanker.score's nesting contract: when a selector is given it runs per fold on
    the training slice only, so every ablation delta-AUC is as leakage-free as the headline.

    ``groups`` switches to paralogue-blocked folds. It is threaded to the ABLATION path as
    well as the headline deliberately: a blocked headline beside unblocked deltas would
    report a corrected level and an uncorrected claim.
    """
    y = np.asarray(y)
    cv = _cv(folds, seed, groups)
    if selector is None:
        return cross_val_predict(
            _gbm(), matrix[cols].values, y, cv=cv, groups=groups, method="predict_proba"
        )[:, 1]

    oof = np.full(len(y), np.nan)
    splits = (
        cv.split(matrix[cols].values, y)
        if groups is None
        else cv.split(matrix[cols].values, y, groups)
    )
    for train_idx, test_idx in splits:
        X_tr = matrix.iloc[train_idx]
        sel = list(selector(X_tr, y[train_idx], list(cols)))
        if not sel:
            oof[test_idx] = y[train_idx].mean()
            continue
        m = _gbm().fit(X_tr[sel].values, y[train_idx])
        oof[test_idx] = m.predict_proba(matrix.iloc[test_idx][sel].values)[:, 1]
    return oof


@dataclass
class EvalResult:
    auc: float
    pr_auc: float
    brier: float
    ablation: pd.DataFrame
    calibration: tuple


def _boot_ci(y, p1, p0, n=1000):
    """Paired bootstrap over GENES for the delta-AUC of ``p1`` over ``p0``.

    Draws that end up single-class are skipped rather than counted: AUC is undefined there,
    and substituting 0.5 would drag the interval toward no-difference.
    """
    rng = np.random.default_rng(42)
    idx = np.arange(len(y))
    d = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), True)
        if 0 < y[b].sum() < len(b):
            d.append(roc_auc_score(y[b], p1[b]) - roc_auc_score(y[b], p0[b]))
    return [round(v, 3) for v in np.percentile(d, [2.5, 50, 97.5])]


class Evaluator:
    def evaluate(
        self, matrix, feat_cols, y, scores, axis_groups, baseline_axis, selector=None,
        groups=None,
    ):
        y = np.asarray(y)
        scores = np.asarray(scores)
        auc = roc_auc_score(y, scores)
        pr = average_precision_score(y, scores)
        cs = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
        brier = brier_score_loss(y, cs)
        cal = calibration_curve(y, cs, n_bins=8, strategy="quantile")

        base_cols = axis_groups[baseline_axis]
        p_base = _raw_oof(matrix, base_cols, y, selector, groups=groups)
        base_auc = roc_auc_score(y, p_base)
        rows = [
            {
                "family": baseline_axis,
                "auc": round(base_auc, 3),
                "d_lo": 0.0,
                "d_md": 0.0,
                "d_hi": 0.0,
            }
        ]
        for fam, cols in axis_groups.items():
            if fam == baseline_axis:
                continue
            cc = list(dict.fromkeys(base_cols + cols))
            p = _raw_oof(matrix, cc, y, selector, groups=groups)
            lo, md, hi = _boot_ci(y, p, p_base)
            rows.append(
                {
                    "family": fam,
                    "auc": round(roc_auc_score(y, p), 3),
                    "d_lo": lo,
                    "d_md": md,
                    "d_hi": hi,
                }
            )
        return EvalResult(
            round(auc, 3), round(pr, 3), round(brier, 4), pd.DataFrame(rows), cal
        )
