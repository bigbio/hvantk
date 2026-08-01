"""Proof that selection is nested inside cross-validation.

If the label-aware filter runs once on the whole matrix and the model is then
cross-validated on the survivors, the reported AUC is optimistically biased -- the
selection has already seen every label (Ambroise & McLachlan 2002). With enough pure-noise
features, that bias alone produces a healthy-looking AUC from data with no signal at all.

These tests fail loudly if anyone later hoists selection out of the fold loop for speed.
"""
import numpy as np
import pandas as pd


def _noise(n=600, p=300, seed=0):
    rng = np.random.default_rng(seed)
    y = np.repeat([0, 1], n // 2)
    X = pd.DataFrame({f"f{i}": rng.normal(0, 1, n) for i in range(p)})
    return X, y


def _top_k(X, y, columns, k=10):
    """Naive practice: keep the k strongest columns by |AUC-0.5|, with no significance gate.

    This is deliberately NOT the shipped policy. The bias demonstration below needs a
    selector that always returns something, and `SelectionPolicy` refuses to (see
    `test_the_shipped_fdr_gate_selects_nothing_from_pure_noise`). Top-k is how feature
    selection is most often done in practice, which is exactly why the bias is worth
    pinning.
    """
    from hvantk.algorithms.rerank.selection import univariate_filter

    stats = univariate_filter(X, y, list(columns), q=1.0)
    return sorted(columns, key=lambda c: -abs(stats[c].auc - 0.5))[:k]


def test_the_shipped_fdr_gate_selects_nothing_from_pure_noise():
    """The default policy's own protection against the bias, pinned.

    BH-FDR at q=0.50 over 300 null columns keeps zero of them -- which is why the two
    tests that follow have to use a naive top-k selector to demonstrate the bias at all.
    If this ever starts keeping columns, the univariate gate has lost its calibration and
    the nesting guardrail is carrying weight it was never meant to carry alone.
    """
    from hvantk.algorithms.rerank.selection import SelectionPolicy, select_axis

    X, y = _noise()
    kept = select_axis(X, y, list(X.columns), SelectionPolicy(wrapper="none", q=0.50)).kept
    assert kept == ()


def test_global_selection_on_pure_noise_looks_predictive():
    """Establishes that the bias this design guards against is real and large.

    Pure noise, no signal whatsoever: selecting the top 10 columns once on the full
    matrix and then cross-validating the survivors yields ~0.58, because the selection
    has already read every label.
    """
    from sklearn.metrics import roc_auc_score

    from hvantk.algorithms.rerank.reranker import ReRanker

    X, y = _noise()
    kept = _top_k(X, y, list(X.columns))
    scores = ReRanker(folds=5).score(X, kept, y)
    assert roc_auc_score(y, scores) > 0.55


def test_nested_selection_on_pure_noise_is_chance():
    """The same data and the same selector, moved inside the folds, must give ~0.5."""
    from sklearn.metrics import roc_auc_score

    from hvantk.algorithms.rerank.reranker import ReRanker

    X, y = _noise()
    scores = ReRanker(folds=5).score(X, list(X.columns), y, selector=_top_k)
    assert roc_auc_score(y, scores) < 0.56


def test_selector_never_receives_held_out_rows():
    """Direct assertion on the contract, independent of any metric."""
    from hvantk.algorithms.rerank.reranker import ReRanker

    X, y = _noise(n=200, p=10, seed=3)
    seen_sizes = []

    def selector(X_tr, y_tr, columns):
        seen_sizes.append(len(X_tr))
        return list(columns)

    ReRanker(folds=5).score(X, list(X.columns), y, selector=selector)
    assert seen_sizes, "selector was never called"
    assert max(seen_sizes) < len(X), "selector saw the full matrix -- selection is not nested"


def test_selection_none_is_unchanged():
    """Backward compatibility: no selector must reproduce the current code path."""
    from hvantk.algorithms.rerank.reranker import ReRanker

    X, y = _noise(n=300, p=20, seed=5)
    a = ReRanker(folds=5).score(X, list(X.columns), y)
    b = ReRanker(folds=5).score(X, list(X.columns), y, selector=None)
    assert np.allclose(a, b)


def test_empty_fold_selection_falls_back_to_training_prevalence():
    """A fold where nothing survives must yield a constant, not a crash or a NaN.

    Pinned separately because it is the one branch the noise tests never reach: top-k
    always returns k columns. A selector that returns nothing is a legitimate outcome on
    a narrow axis -- the shipped FDR gate does exactly that on null data -- and the
    out-of-fold vector must still come back fully populated, or every downstream AUC
    silently becomes NaN.
    """
    from hvantk.algorithms.rerank.reranker import ReRanker

    X, y = _noise(n=200, p=6, seed=11)
    scores = ReRanker(folds=5).score(X, list(X.columns), y, selector=lambda *a: [])
    assert not np.isnan(scores).any()
    assert np.allclose(scores, 0.5, atol=0.05)


def test_evaluator_raw_oof_honours_the_same_contract():
    """The ablation path must be as leakage-free as the headline path.

    `_raw_oof` is what every per-axis delta-AUC is computed from. If the selector were
    threaded into `ReRanker.score` but not here, the headline number would be nested and
    the ablation numbers -- the actual scientific claim -- would not be.
    """
    from hvantk.algorithms.rerank.evaluator import _raw_oof

    X, y = _noise(n=200, p=10, seed=7)
    seen_sizes = []

    def selector(X_tr, y_tr, columns):
        seen_sizes.append(len(X_tr))
        return list(columns)

    _raw_oof(X, list(X.columns), y, selector=selector)
    assert seen_sizes, "selector was never called"
    assert max(seen_sizes) < len(X), "_raw_oof selected on the full matrix"
