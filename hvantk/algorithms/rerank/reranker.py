# local/rerank_engine/reranker.py
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_val_predict


def _gbm():   # ported from chd_score_lib.gbm()
    return HistGradientBoostingClassifier(max_depth=3, max_iter=250, learning_rate=0.05,
        l2_regularization=1.0, min_samples_leaf=20, class_weight="balanced", random_state=42)


class ReRanker:
    """Calibrated GBM scorer: inner isotonic calibration nested inside an outer
    5-fold cross_val_predict for out-of-fold, no-leakage calibrated probabilities.
    Matches Phase-1 chd_calibrate.py exactly."""
    def __init__(self, calibration="isotonic", folds=5):
        self.calibration = calibration
        self.folds = folds

    def score(self, matrix, feat_cols, y):
        X = matrix[feat_cols].values
        cv = StratifiedKFold(self.folds, shuffle=True, random_state=42)
        clf = CalibratedClassifierCV(_gbm(), method=self.calibration, cv=self.folds)
        return cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
