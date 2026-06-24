# local/rerank_engine/evaluator.py
from dataclasses import dataclass
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.calibration import calibration_curve
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from hvantk.algorithms.rerank.reranker import ReRanker


def _raw_oof(matrix, cols, y):
    gbm = HistGradientBoostingClassifier(max_depth=3, max_iter=250, learning_rate=0.05,
        l2_regularization=1.0, min_samples_leaf=20, class_weight="balanced", random_state=42)
    cv = StratifiedKFold(5, shuffle=True, random_state=42)
    return cross_val_predict(gbm, matrix[cols].values, y, cv=cv, method="predict_proba")[:, 1]

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
    def evaluate(self, matrix, feat_cols, y, scores, axis_groups, baseline_axis):
        y = np.asarray(y); scores = np.asarray(scores)
        auc = roc_auc_score(y, scores); pr = average_precision_score(y, scores)
        cs = (scores - scores.min())/(scores.max()-scores.min()+1e-12)
        brier = brier_score_loss(y, cs)
        cal = calibration_curve(y, cs, n_bins=8, strategy="quantile")
        base_cols = axis_groups[baseline_axis]
        p_base = _raw_oof(matrix, base_cols, y); base_auc = roc_auc_score(y, p_base)
        rows = [{"family": baseline_axis, "auc": round(base_auc,3), "d_lo":0.0,"d_md":0.0,"d_hi":0.0}]
        for fam, cols in axis_groups.items():
            if fam == baseline_axis: continue
            cc = list(dict.fromkeys(base_cols + cols))
            p = _raw_oof(matrix, cc, y); lo, md, hi = _boot_ci(y, p, p_base)
            rows.append({"family": fam, "auc": round(roc_auc_score(y,p),3), "d_lo":lo,"d_md":md,"d_hi":hi})
        return EvalResult(round(auc,3), round(pr,3), round(brier,4), pd.DataFrame(rows), cal)
