# local/rerank_engine/evaluator.py
from dataclasses import dataclass
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.calibration import calibration_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from hvantk.algorithms.rerank.reranker import ReRanker, _gbm


def _raw_oof(matrix, cols, y, selector=None):
    """Uncalibrated out-of-fold probabilities for one feature subset.

    Mirrors ReRanker.score's nesting contract: when a selector is given it runs per fold on
    the training slice only, so every ablation delta-AUC is as leakage-free as the headline.
    """
    y = np.asarray(y)
    cv = StratifiedKFold(5, shuffle=True, random_state=42)
    if selector is None:
        return cross_val_predict(_gbm(), matrix[cols].values, y, cv=cv,
                                 method="predict_proba")[:, 1]

    oof = np.full(len(y), np.nan)
    for train_idx, test_idx in cv.split(matrix[cols].values, y):
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
    auc: float; pr_auc: float; brier: float; ablation: pd.DataFrame; calibration: tuple

def _boot_ci(y, p1, p0, n=1000):
    rng = np.random.default_rng(42); idx = np.arange(len(y)); d=[]
    for _ in range(n):
        b = rng.choice(idx, len(idx), True)
        if 0 < y[b].sum() < len(b): d.append(roc_auc_score(y[b], p1[b]) - roc_auc_score(y[b], p0[b]))
    return [round(v,3) for v in np.percentile(d, [2.5,50,97.5])]

class Evaluator:
    def evaluate(self, matrix, feat_cols, y, scores, axis_groups, baseline_axis, selector=None):
        y = np.asarray(y); scores = np.asarray(scores)
        auc = roc_auc_score(y, scores); pr = average_precision_score(y, scores)
        cs = (scores - scores.min())/(scores.max()-scores.min()+1e-12)
        brier = brier_score_loss(y, cs)
        cal = calibration_curve(y, cs, n_bins=8, strategy="quantile")
        base_cols = axis_groups[baseline_axis]
        p_base = _raw_oof(matrix, base_cols, y, selector); base_auc = roc_auc_score(y, p_base)
        rows = [{"family": baseline_axis, "auc": round(base_auc,3), "d_lo":0.0,"d_md":0.0,"d_hi":0.0}]
        for fam, cols in axis_groups.items():
            if fam == baseline_axis: continue
            cc = list(dict.fromkeys(base_cols + cols))
            p = _raw_oof(matrix, cc, y, selector); lo, md, hi = _boot_ci(y, p, p_base)
            rows.append({"family": fam, "auc": round(roc_auc_score(y,p),3), "d_lo":lo,"d_md":md,"d_hi":hi})
        return EvalResult(round(auc,3), round(pr,3), round(brier,4), pd.DataFrame(rows), cal)
