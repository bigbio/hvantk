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
