import numpy as np
import pandas as pd
import pytest


def test_auc_is_directionless_and_rank_based():
    """Anti-predictive features are informative; |AUC-0.5| must not discard them."""
    from hvantk.algorithms.rerank.selection import univariate_auc

    y = np.array([0, 0, 1, 1])
    up = np.array([1.0, 2.0, 3.0, 4.0])
    down = np.array([4.0, 3.0, 2.0, 1.0])
    assert univariate_auc(up, y)[0] == pytest.approx(1.0)
    assert univariate_auc(down, y)[0] == pytest.approx(0.0)


def test_auc_ignores_monotone_rescaling():
    """Rank-based, so rankscore vs phred vs raw compare fairly."""
    from hvantk.algorithms.rerank.selection import univariate_auc

    y = np.array([0, 0, 1, 1])
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert univariate_auc(x, y)[0] == pytest.approx(univariate_auc(np.exp(x) * 1000, y)[0])


def test_auc_uses_only_rows_where_the_feature_is_defined():
    """NaN rows are excluded and the effective n is reported, so the z-score self-corrects."""
    from hvantk.algorithms.rerank.selection import univariate_auc

    y = np.array([0, 0, 1, 1])
    x = np.array([1.0, np.nan, 3.0, 4.0])
    auc, n_pos, n_neg = univariate_auc(x, y)
    assert (n_pos, n_neg) == (2, 1)
    assert auc == pytest.approx(1.0)


def test_all_nan_feature_is_reported_not_crashed():
    from hvantk.algorithms.rerank.selection import univariate_auc

    auc, n_pos, n_neg = univariate_auc(np.full(4, np.nan), np.array([0, 0, 1, 1]))
    assert np.isnan(auc) and n_pos == 0 and n_neg == 0


def test_fdr_is_stricter_for_a_wider_axis():
    """A 50-column axis is held to a stricter bar than a 1-column axis, automatically."""
    from hvantk.algorithms.rerank.selection import univariate_filter

    rng = np.random.default_rng(42)
    n = 4000
    y = np.repeat([0, 1], n // 2)
    signal = y + rng.normal(0, 3.0, n)          # weak but real
    noise = {f"noise{i}": rng.normal(0, 1, n) for i in range(49)}

    narrow = pd.DataFrame({"signal": signal})
    wide = pd.DataFrame({"signal": signal, **noise})

    r_narrow = univariate_filter(narrow, y, ["signal"], q=0.10)
    r_wide = univariate_filter(wide, y, list(wide.columns), q=0.10)

    assert r_narrow["signal"].passed
    # same feature, same data, more columns -> multiplicity correction is harsher
    assert r_wide["signal"].p == pytest.approx(r_narrow["signal"].p)
    assert sum(s.passed for s in r_wide.values()) < len(wide.columns)


def test_pure_noise_axis_admits_almost_nothing():
    from hvantk.algorithms.rerank.selection import univariate_filter

    rng = np.random.default_rng(0)
    n = 2000
    y = np.repeat([0, 1], n // 2)
    X = pd.DataFrame({f"n{i}": rng.normal(0, 1, n) for i in range(40)})
    res = univariate_filter(X, y, list(X.columns), q=0.10)
    assert sum(s.passed for s in res.values()) <= 2


def test_tied_feature_values_use_average_rank():
    """Real feature columns are full of ties (rankscore buckets, phred caps, integer
    read depths). rankdata's default 'average' handling must be what's used here, or
    AUC would silently depend on how ties happen to break instead of only on the data.
    """
    from hvantk.algorithms.rerank.selection import univariate_auc

    x = np.array([1.0, 2.0, 2.0, 3.0])
    y = np.array([0, 0, 1, 1])
    auc, n_pos, n_neg = univariate_auc(x, y)
    assert auc == pytest.approx(0.875)
    assert (n_pos, n_neg) == (2, 2)


def test_single_class_y_is_reported_not_passing_not_crashed():
    """A training slice can legitimately land on a single class for some fold/axis
    combination (a tiny cohort, a rare-disease arm). AUC is undefined there; the
    filter must report the column as not passing rather than raising, so one
    degenerate fold doesn't take down the whole selection step.
    """
    from hvantk.algorithms.rerank.selection import univariate_filter

    X = pd.DataFrame({"col": [1.0, 2.0, 3.0, 4.0]})
    all_pos = np.array([1, 1, 1, 1])
    all_neg = np.array([0, 0, 0, 0])

    res_pos = univariate_filter(X, all_pos, ["col"], q=0.10)
    res_neg = univariate_filter(X, all_neg, ["col"], q=0.10)

    assert res_pos["col"].passed is False
    assert res_neg["col"].passed is False
    assert res_pos["col"].p == 1.0 and res_neg["col"].p == 1.0


def test_empty_column_list_returns_empty_dict():
    """An axis can legitimately contribute zero columns to a slice (e.g. everything
    upstream already dropped as all-missing); the filter must hand back an empty
    result instead of raising on an empty `columns` (and the internal empty `order`).
    """
    from hvantk.algorithms.rerank.selection import univariate_filter

    X = pd.DataFrame({"col": [1.0, 2.0, 3.0, 4.0]})
    y = np.array([0, 0, 1, 1])
    assert univariate_filter(X, y, [], q=0.10) == {}


def test_bh_admits_more_weak_signals_than_a_fixed_bonferroni_bound():
    """BH's pass bar grows with rank (q*rank/m); a fixed per-column bar (q/m, i.e.
    Bonferroni) would silently reject real-but-modest effects that BH is specifically
    designed to keep. `test_fdr_is_stricter_for_a_wider_axis` and
    `test_pure_noise_axis_admits_almost_nothing` above both still pass unmodified if
    `univariate_filter` were "simplified" to Bonferroni -- one obvious signal, or none
    at all, doesn't distinguish the two procedures. Many weak-but-real signals mixed
    with pure noise does: BH's rank-adaptive threshold must admit strictly more
    columns than a fixed per-column bound would, in this regime.
    """
    from hvantk.algorithms.rerank.selection import univariate_filter

    rng = np.random.default_rng(51)
    n = 400
    y = np.repeat([0, 1], n // 2)
    signal = {f"signal{i}": y + rng.normal(0, 2.5, n) for i in range(20)}   # weak but real
    noise = {f"noise{i}": rng.normal(0, 1, n) for i in range(30)}          # pure noise
    X = pd.DataFrame({**signal, **noise})

    q = 0.10
    stats = univariate_filter(X, y, list(X.columns), q=q)
    bh_passed = sum(s.passed for s in stats.values())

    # Same p-values, scored against Bonferroni's fixed q/m instead of BH's q*rank/m --
    # the exact substitution that the two tests above cannot tell apart from BH.
    testable = [s for s in stats.values() if np.isfinite(s.z)]
    m = len(testable)
    bonferroni_passed = sum(1 for s in testable if s.p <= q / m)

    assert bh_passed > bonferroni_passed
