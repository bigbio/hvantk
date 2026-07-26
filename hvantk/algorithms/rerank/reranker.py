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

    def score(self, matrix, feat_cols, y, selector=None):
        """Out-of-fold calibrated probabilities.

        With ``selector`` given, feature selection re-runs inside every fold on the
        TRAINING slice only. The fold loop is explicit because ``cross_val_predict``
        cannot express per-fold feature selection -- it takes a fixed design matrix.

        With ``selector=None`` this is the original single ``cross_val_predict`` call and
        the output is unchanged.
        """
        y = np.asarray(y)
        cv = StratifiedKFold(self.folds, shuffle=True, random_state=42)
        if selector is None:
            X = matrix[feat_cols].values
            clf = CalibratedClassifierCV(_gbm(), method=self.calibration, cv=self.folds)
            return cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]

        oof = np.full(len(y), np.nan)
        for train_idx, test_idx in cv.split(matrix[feat_cols].values, y):
            X_tr = matrix.iloc[train_idx]
            cols = list(selector(X_tr, y[train_idx], list(feat_cols)))
            if not cols:
                # Nothing survived on this fold: predict the training prevalence rather
                # than crash. A fold that selects nothing is information, not an error.
                oof[test_idx] = y[train_idx].mean()
                continue
            clf = CalibratedClassifierCV(_gbm(), method=self.calibration, cv=self.folds)
            clf.fit(X_tr[cols].values, y[train_idx])
            oof[test_idx] = clf.predict_proba(matrix.iloc[test_idx][cols].values)[:, 1]
        return oof
