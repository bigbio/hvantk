# local/rerank_engine/reranker.py
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold, cross_val_predict

from hvantk.algorithms.rerank.seeds import DEFAULT_SEED


def _gbm(seed: int = DEFAULT_SEED):   # ported from chd_score_lib.gbm()
    return HistGradientBoostingClassifier(max_depth=3, max_iter=250, learning_rate=0.05,
        l2_regularization=1.0, min_samples_leaf=20, class_weight="balanced",
        random_state=seed)


class _SeededStratifiedGroupKFold(StratifiedGroupKFold):
    """``StratifiedGroupKFold(shuffle=False)`` with a seeded random relabelling of the
    group ids.

    sklearn 1.7.2's own ``shuffle=True`` shuffles the per-group class-count ROWS before
    deciding each group's fold, but then records that decision under the group's
    UNshuffled index -- so a group's fold is chosen from another group's class counts,
    not its own, and the class balance across folds stops being reliable (with singleton
    groups the positives-per-fold spread becomes several times wider than with
    ``shuffle=False``). Relabelling the group IDS with a seeded permutation before calling
    ``shuffle=False`` sidesteps that defect entirely: the balancing step always sees each
    group's own, correct counts, while randomising the ID order still randomises the
    tie-break among groups whose class-distribution spread is equal -- which is what the
    upstream shuffle was meant to do. Partitions therefore will not change silently if
    sklearn's own shuffle is ever fixed.
    """

    def __init__(self, n_splits=5, seed=None):
        super().__init__(n_splits=n_splits, shuffle=False)
        self.seed = seed

    def split(self, X, y=None, groups=None):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        _, inv = np.unique(groups, return_inverse=True)
        relabel = np.random.default_rng(self.seed).permutation(inv.max() + 1)
        return super().split(X, y, relabel[inv])


def _cv(folds, seed, groups):
    """Blocked folds when ``groups`` is given, stratified random folds otherwise.

    ``StratifiedGroupKFold`` must place an entire group in one fold, so it trades fold
    BALANCE for the guarantee that no block straddles the split. Grouped folds are harder
    than random ones and every absolute AUC is expected to fall; the quantity of interest is
    whether a delta survives, not the level.

    Uses ``_SeededStratifiedGroupKFold`` rather than plain ``StratifiedGroupKFold(shuffle=
    True)``: the latter does not actually stratify (see its docstring above), so this is
    what makes the word "stratified" in "blocked folds" true.
    """
    if groups is None:
        return StratifiedKFold(folds, shuffle=True, random_state=seed)
    return _SeededStratifiedGroupKFold(folds, seed=seed)


def _grouped_splits(cv, X, y, groups):
    """The blocked splits, materialised once and checked.

    ``StratifiedGroupKFold`` balances per-fold class PROPORTIONS over whole blocks, but its
    own check operates on class sizes in SAMPLES rather than in blocks, so it has no way to
    notice when a class is confined to too few blocks for the requested fold count.
    Left unchecked that surfaces two different ways: a training fold with no example of a
    class -- a classifier fit on a single class predicts that class with probability 1, so
    every held-out positive silently scores 0.0, with nothing louder than a sklearn
    UserWarning -- or, with fewer blocks than folds, a test fold with no example at all,
    which surfaces many estimator calls later as an opaque "Found array with 0 sample(s)".
    Both are refused HERE, before either can reach an estimator.
    """
    y = np.asarray(y)
    n_folds = cv.get_n_splits()
    n_blocks = len(np.unique(groups))
    if n_blocks < n_folds:
        raise ValueError(
            f"only {n_blocks} paralogue block(s) for {n_folds} folds -- lower the fold "
            "count, or run without blocks (blocks come from the universe, not the labels, "
            "so widening the label set cannot add more of them)"
        )
    splits = list(cv.split(X, y, groups))
    for i, (train_idx, test_idx) in enumerate(splits):
        if len(test_idx) == 0:
            raise ValueError(
                f"blocked fold {i} has no test example: the paralogue blocks split too "
                f"unevenly across {n_folds} folds -- lower the fold count, or run without "
                "blocks"
            )
        train_y = y[train_idx]
        if train_y.size == 0:
            raise ValueError(
                f"blocked fold {i} has no training example at all: the paralogue blocks "
                f"split too unevenly across {n_folds} folds -- lower the fold count, or "
                "run without blocks"
            )
        if train_y.min() == train_y.max():
            label = "positive" if train_y.max() == 0 else "negative"
            raise ValueError(
                f"blocked fold {i} has no {label} training example: the {label}s are "
                f"confined to too few of the paralogue blocks to cover {n_folds} folds -- "
                "lower the fold count, widen the label set, or run without blocks"
            )
    return splits


class ReRanker:
    """Calibrated GBM scorer: inner isotonic calibration nested inside an outer
    5-fold cross_val_predict for out-of-fold, no-leakage calibrated probabilities.
    Matches Phase-1 chd_calibrate.py exactly."""
    def __init__(self, calibration="isotonic", folds=5, seed=DEFAULT_SEED):
        self.calibration = calibration
        self.folds = folds
        self.seed = seed

    def score(self, matrix, feat_cols, y, selector=None, groups=None):
        """Out-of-fold calibrated probabilities.

        With ``selector`` given, feature selection re-runs inside every fold on the
        TRAINING slice only. The fold loop is explicit because ``cross_val_predict``
        cannot express per-fold feature selection -- it takes a fixed design matrix.

        With ``selector=None`` this is the original single ``cross_val_predict`` call and
        the output is unchanged.

        ``groups`` blocks the OUTER folds: no held-out block can leak into any prediction.
        Known and deliberate limit: the inner isotonic calibration split inside
        ``CalibratedClassifierCV`` takes an integer ``cv`` and so cannot see the blocks,
        which means a family CAN straddle the calibration split WITHIN one training fold.
        The per-axis deltas, which are the scientific claim, go through
        ``evaluator._raw_oof`` instead, which is uncalibrated and fully blocked. With
        ``SelectionPolicy(wrapper="rfecv")``, RFECV's own inner ``StratifiedKFold`` inside
        each training slice is likewise unblocked -- an opt-in path with the same kind of
        limit. Documented here rather than silently accepted.
        """
        y = np.asarray(y)
        cv = _cv(self.folds, self.seed, groups)
        if selector is None:
            X = matrix[feat_cols].values
            clf = CalibratedClassifierCV(_gbm(self.seed), method=self.calibration, cv=self.folds)
            cv_arg = cv if groups is None else _grouped_splits(cv, X, y, groups)
            return cross_val_predict(
                clf, X, y, cv=cv_arg, groups=groups, method="predict_proba"
            )[:, 1]

        oof = np.full(len(y), np.nan)
        X = matrix[feat_cols].values
        splits = cv.split(X, y) if groups is None else _grouped_splits(cv, X, y, groups)
        for train_idx, test_idx in splits:
            X_tr = matrix.iloc[train_idx]
            cols = list(selector(X_tr, y[train_idx], list(feat_cols)))
            if not cols:
                # Nothing survived on this fold: predict the training prevalence rather
                # than crash. A fold that selects nothing is information, not an error.
                oof[test_idx] = y[train_idx].mean()
                continue
            clf = CalibratedClassifierCV(_gbm(self.seed), method=self.calibration, cv=self.folds)
            clf.fit(X_tr[cols].values, y[train_idx])
            oof[test_idx] = clf.predict_proba(matrix.iloc[test_idx][cols].values)[:, 1]
        return oof
