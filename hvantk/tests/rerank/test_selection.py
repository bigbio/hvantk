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


def test_redundancy_keeps_the_stronger_member_of_a_correlated_pair():
    """Greedy by univariate strength: the better feature survives, the other is named."""
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(1)
    base = rng.normal(0, 1, 500)
    X = pd.DataFrame({"strong": base, "copy": base * 3.0 + 1.0, "other": rng.normal(0, 1, 500)})
    kept, dropped = redundancy_filter(
        X, ["strong", "copy", "other"], {"strong": 0.30, "copy": 0.10, "other": 0.20}, 0.75
    )
    assert "strong" in kept and "other" in kept
    assert "copy" not in kept
    assert dropped["copy"] == "redundant_with:strong"


def test_redundancy_uses_spearman_so_monotone_rescaling_still_collapses():
    """These scores are monotonically but not linearly related; Pearson would miss it."""
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(2)
    base = rng.uniform(0.1, 5.0, 400)
    X = pd.DataFrame({"a": base, "b": np.exp(base)})
    kept, _ = redundancy_filter(X, ["a", "b"], {"a": 0.30, "b": 0.20}, 0.75)
    assert kept == ["a"]


def test_uncorrelated_columns_all_survive():
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(3)
    X = pd.DataFrame({c: rng.normal(0, 1, 300) for c in "abc"})
    kept, dropped = redundancy_filter(X, list("abc"), {c: 0.1 for c in "abc"}, 0.75)
    assert set(kept) == set("abc") and dropped == {}


def test_redundancy_stronger_column_survives_regardless_of_list_order():
    """Same fixture as `test_redundancy_keeps_the_stronger_member_of_a_correlated_pair`
    above, but with the weaker/duplicate column listed FIRST in `columns`. All three
    shipped fixtures happen to list the highest-scoring column first, so an
    implementation that walks `columns` in caller-supplied order -- instead of sorting
    by descending |score| -- passes them unnoticed. Reordering the input list must not
    change which column survives: the greedy walk is defined over score order, never
    list order.
    """
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(1)
    base = rng.normal(0, 1, 500)
    X = pd.DataFrame({"strong": base, "copy": base * 3.0 + 1.0, "other": rng.normal(0, 1, 500)})
    kept, dropped = redundancy_filter(
        X, ["copy", "strong", "other"], {"strong": 0.30, "copy": 0.10, "other": 0.20}, 0.75
    )
    assert "strong" in kept and "other" in kept
    assert "copy" not in kept
    assert dropped["copy"] == "redundant_with:strong"


def test_redundancy_spearman_keeps_a_pair_pearson_would_drop():
    """`test_redundancy_uses_spearman_so_monotone_rescaling_still_collapses` above uses
    b = exp(a), but for that fixture Pearson(a, b) = 0.868 -- still over the 0.75
    cutoff, so `b` is dropped whether the implementation uses Spearman or Pearson, and
    the test cannot tell the two apart. Here the transform is steep enough that
    Pearson(a, b) = 0.627 (clearly under cutoff, so a Pearson-based filter would keep
    both columns) while Spearman(a, b) is exactly 1.0 (a monotone transform preserves
    rank order perfectly, so it stays over cutoff regardless of shape). Only a
    Spearman-based filter drops `b` here.
    """
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(4)
    base = rng.uniform(0.1, 5.0, 400)
    X = pd.DataFrame({"a": base, "b": np.exp(3 * base)})
    assert X["a"].corr(X["b"], method="pearson") < 0.75
    assert X["a"].corr(X["b"], method="spearman") >= 0.75

    kept, dropped = redundancy_filter(X, ["a", "b"], {"a": 0.30, "b": 0.20}, 0.75)
    assert kept == ["a"]
    assert dropped["b"] == "redundant_with:a"


def test_redundancy_tiny_overlap_is_not_declared_redundant():
    """A sparse dbNSFP-style score can be non-null for only a handful of rows. Two
    columns that happen to agree on their only 2 shared non-null rows are trivially
    "perfectly correlated" -- any 2 distinct points are perfectly monotonic -- so
    without the `len(pair) < 3` guard this reads as rho ~= 1.0 and the sparse column
    would be wrongly discarded on the strength of two observations.
    """
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    X = pd.DataFrame(
        {
            "dense": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            "sparse": [10.0, 20.0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        }
    )
    kept, dropped = redundancy_filter(X, ["dense", "sparse"], {"dense": 0.30, "sparse": 0.20}, 0.75)
    assert set(kept) == {"dense", "sparse"}
    assert dropped == {}


def test_redundancy_checks_survivors_not_everything_walked_so_far():
    """A-B and B-C are each over cutoff, but A-C is not. B is dropped as redundant
    with A; by the time C is considered, B is no longer a *kept* column, so C must be
    compared only against A (still under cutoff) and must survive. Comparing against
    every column visited so far -- including the discarded B -- would wrongly drop C
    too, even though C's only redundancy was with a column that itself got dropped.
    """
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(9)
    corr = np.array(
        [
            [1.00, 0.80, 0.60],
            [0.80, 1.00, 0.86],
            [0.60, 0.86, 1.00],
        ]
    )
    samples = rng.multivariate_normal(mean=[0.0, 0.0, 0.0], cov=corr, size=300)
    X = pd.DataFrame(samples, columns=["A", "B", "C"])
    assert X["A"].corr(X["B"], method="spearman") >= 0.75
    assert X["B"].corr(X["C"], method="spearman") >= 0.75
    assert X["A"].corr(X["C"], method="spearman") < 0.75

    kept, dropped = redundancy_filter(X, ["A", "B", "C"], {"A": 0.30, "B": 0.20, "C": 0.10}, 0.75)
    assert kept == ["A", "C"]
    assert dropped == {"B": "redundant_with:A"}


def test_redundancy_empty_columns_returns_empty():
    """An axis can legitimately contribute zero columns to a slice after upstream
    filtering; the filter must hand back an empty result instead of raising on an
    empty `columns` (and the internal empty `ordered`/`kept`).
    """
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    X = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
    assert redundancy_filter(X, [], {}, 0.75) == ([], {})


def test_redundancy_missing_score_falls_back_to_zero():
    """`scores` comes from the univariate step and may legitimately omit a column
    (e.g. it was untestable there and never got a score). The greedy walk must still
    place such a column instead of raising a KeyError, treating an absent score as the
    weakest possible (0.0) rather than crashing the whole selection step.
    """
    import numpy as np
    import pandas as pd

    from hvantk.algorithms.rerank.selection import redundancy_filter

    rng = np.random.default_rng(20)
    X = pd.DataFrame({"a": rng.normal(0, 1, 300), "b": rng.normal(0, 1, 300)})
    kept, dropped = redundancy_filter(X, ["a", "b"], {"a": 0.5}, 0.75)
    assert set(kept) == {"a", "b"}
    assert dropped == {}
