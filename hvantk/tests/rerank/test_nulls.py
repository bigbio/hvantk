"""Pass A of #247: the permutation null, and the p-value a finite null can license.

`grep -rn "permut\\|selected_max" hvantk/algorithms/rerank/` returned nothing before this.
The selected-maximum null has reversed four results in the motivating project and its
median sits at +0.012 to +0.027 -- larger than most axes' entire measured gain -- so
"axis X adds +0.02" was not interpretable without it.
"""
from __future__ import annotations

import numpy as np
import pytest

from hvantk.algorithms.rerank.nulls import (
    ControlSetting,
    NullConfig,
    oof_scorer,
    p_value,
    permutation_deltas,
)
from hvantk.tests.rerank._synth import cheap_scorer, permuted_labels, planted_signal


def _setting(baseline=("base",), **kw):
    return ControlSetting(
        arm=kw.get("arm", "all"),
        leakage=kw.get("leakage", False),
        baseline=tuple(baseline),
        folds=kw.get("folds", 5),
        blocked=kw.get("blocked", False),
    )


# --- the p-value a finite permutation set can license -----------------------------------


def test_p_value_can_never_be_zero():
    """(1 + #{null >= obs}) / (1 + n_perm). The analysis drivers used the plain mean
    (analysis/crossconfig-multiplicity/run.py:152), which reports p = 0.000 for an observed
    value no permutation reached -- a claim 200 permutations cannot support."""
    null = np.zeros(199)
    assert p_value(null, observed=10.0) == pytest.approx(1 / 200)
    assert p_value(null, observed=10.0) > 0.0


def test_p_value_counts_ties_as_at_least_as_extreme():
    null = np.array([0.0, 0.2, 0.2, 0.2])
    # >=, not >: three of four draws reach 0.2.
    assert p_value(null, observed=0.2) == pytest.approx(4 / 5)


def test_p_value_ignores_nan_draws_rather_than_counting_them():
    """An axis wholly contained in the baseline contributes NaN, not 0.0: a delta of
    exactly zero is a measurement, absence is not."""
    assert p_value(np.array([0.1, np.nan, 0.3]), observed=0.2) == pytest.approx(2 / 3)


def test_p_value_refuses_an_empty_null():
    with pytest.raises(ValueError, match="no finite"):
        p_value(np.array([np.nan, np.nan]), observed=0.1)


# --- chunking ----------------------------------------------------------------------------


def test_chunks_are_contiguous_disjoint_and_cover_everything():
    spans = [NullConfig(n_perm=200, chunk=c, n_chunks=7).span() for c in range(7)]
    assert spans[0][0] == 0 and spans[-1][1] == 200
    assert all(a[1] == b[0] for a, b in zip(spans, spans[1:]))
    assert sum(hi - lo for lo, hi in spans) == 200


def test_a_chunk_reproduces_exactly_the_permutations_the_whole_run_would_have():
    """The seed of permutation i is seed + i, not a draw from a stream, so rerunning chunk
    3 reproduces permutations [lo, hi) and nothing else
    (analysis/rerank-homogenised/perm_null_arm.py:77-86)."""
    matrix, y, baseline, axes = permuted_labels(n=120, n_noise=1)
    axes = {"axis0": axes["axis0"], "axis1": axes["axis1"]}
    scorer = cheap_scorer()
    whole = permutation_deltas(
        matrix, baseline, axes, y,
        config=NullConfig(n_perm=6, seed=5), scorer=scorer,
    )
    part = permutation_deltas(
        matrix, baseline, axes, y,
        config=NullConfig(n_perm=6, chunk=1, n_chunks=3, seed=5), scorer=scorer,
    )
    assert sorted(part.perm.unique().tolist()) == [2, 3]
    merged = whole[whole.perm.isin([2, 3])].reset_index(drop=True)
    assert np.allclose(
        part.sort_values(["perm", "axis"]).delta.to_numpy(),
        merged.sort_values(["perm", "axis"]).delta.to_numpy(),
    )


def test_n_chunks_of_one_is_the_whole_range():
    assert NullConfig(n_perm=13).span() == (0, 13)


@pytest.mark.parametrize(
    "kw,match",
    [
        (dict(n_perm=0), "n_perm"),
        (dict(n_perm=10, n_chunks=0), "n_chunks"),
        (dict(n_perm=10, chunk=3, n_chunks=3), "chunk"),
        (dict(n_perm=10, chunk=-1), "chunk"),
    ],
)
def test_null_config_rejects_impossible_chunking(kw, match):
    with pytest.raises(ValueError, match=match):
        NullConfig(**kw).span()


# --- the permutation loop -----------------------------------------------------------------


def test_every_permutation_refits_the_baseline():
    """The recorded base_auc must MOVE between permutations. Holding the baseline fixed
    measures the permutation rather than the axis
    (analysis/rerank-homogenised/perm_null_arm.py:15-16)."""
    matrix, y, baseline, axes = permuted_labels(n=120, n_noise=1)
    d = permutation_deltas(
        matrix, baseline, {"axis0": axes["axis0"]}, y,
        config=NullConfig(n_perm=6, seed=1), scorer=cheap_scorer(),
    )
    assert d.base_auc.nunique() > 1, "the baseline was not refit per permutation"


def test_the_observed_labels_are_never_scored():
    """A permutation null that included the real label as a draw would be biased toward
    not rejecting. Every yp handed to the scorer must differ from y."""
    matrix, y, baseline, axes = planted_signal(n=120, n_noise=1)
    seen = []

    def spy(m, cols, yy):
        seen.append(np.asarray(yy).copy())
        return cheap_scorer()(m, cols, yy)

    permutation_deltas(
        matrix, baseline, {"axis0": axes["axis0"]}, y,
        config=NullConfig(n_perm=4, seed=1), scorer=spy,
    )
    assert seen and all(not np.array_equal(s, np.asarray(y)) for s in seen)
    assert all(int(s.sum()) == int(np.asarray(y).sum()) for s in seen), "labels not permuted"


def test_an_axis_wholly_inside_the_baseline_is_nan_not_zero():
    matrix, y, baseline, _ = planted_signal(n=120, n_noise=1)
    d = permutation_deltas(
        matrix, baseline, {"same": ["base"]}, y,
        config=NullConfig(n_perm=2, seed=1), scorer=cheap_scorer(),
    )
    assert d.delta.isna().all()


def test_per_axis_p_value_on_a_planted_signal_clears():
    """The control that must pass: an axis with real signal beats its own null."""
    matrix, y, baseline, axes = planted_signal(n=200, n_noise=1)
    scorer = cheap_scorer()
    from sklearn.metrics import roc_auc_score

    a0 = roc_auc_score(y, scorer(matrix, baseline, y))
    a1 = roc_auc_score(y, scorer(matrix, baseline + axes["axis0"], y))
    d = permutation_deltas(
        matrix, baseline, {"axis0": axes["axis0"]}, y,
        config=NullConfig(n_perm=39, seed=1), scorer=scorer,
    )
    p = p_value(d[d.axis == "axis0"].delta.to_numpy(), a1 - a0)
    assert p <= 0.05, p


def test_per_axis_p_values_are_approximately_uniform_under_permuted_labels():
    """Leave-one-out over the null itself: under exchangeability each draw's p against the
    others is uniform. A loose KS bound (0.01), because the point is to catch a broken
    formula or a leaking baseline, not to certify the tail."""
    from scipy.stats import kstest

    matrix, y, baseline, axes = permuted_labels(n=160, n_noise=1)
    d = permutation_deltas(
        matrix, baseline, {"axis0": axes["axis0"]}, y,
        config=NullConfig(n_perm=40, seed=2), scorer=cheap_scorer(),
    )
    draws = d[d.axis == "axis0"].delta.to_numpy()
    ps = [p_value(np.delete(draws, i), draws[i]) for i in range(len(draws))]
    assert kstest(ps, "uniform").pvalue > 0.01, sorted(ps)


# --- the default scorer is the real thing -------------------------------------------------


def test_default_scorer_is_raw_oof():
    """The cheap scorer above exists for the distributional tests. The DEFAULT must be the
    shipped estimator, or the null is computed with a different model from the deltas."""
    from hvantk.algorithms.rerank.evaluator import _raw_oof

    matrix, y, baseline, _ = planted_signal(n=100, n_noise=0)
    assert np.allclose(oof_scorer()(matrix, baseline, y), _raw_oof(matrix, baseline, y))


def test_permutation_deltas_runs_end_to_end_on_the_shipped_estimator():
    """Two permutations only -- the point is that the default path works, not its shape."""
    matrix, y, baseline, axes = planted_signal(n=100, n_noise=0)
    d = permutation_deltas(
        matrix, baseline, {"axis0": axes["axis0"]}, y, config=NullConfig(n_perm=2, seed=1)
    )
    assert len(d) == 2 and d.delta.notna().all()


# --- the control setting is recorded, and hashable ---------------------------------------


def test_control_setting_is_frozen_and_compares_by_value():
    a, b = _setting(), _setting()
    assert a == b and hash(a) == hash(b)
    assert a != _setting(leakage=True)
    assert a != _setting(baseline=("base", "other"))
    with pytest.raises(Exception):
        a.leakage = True


# --- no hail, ever -------------------------------------------------------------------------


def test_nulls_module_imports_without_hail(monkeypatch):
    """The rerank stack is numpy/pandas/sklearn. Importing hail here would make `hvantk
    rerank` unusable on a machine that has no Spark, and would pull the module into the
    `hail`-marked test selection."""
    import sys

    monkeypatch.setitem(sys.modules, "hail", None)
    for name in [m for m in list(sys.modules) if m.startswith("hvantk.algorithms.rerank")]:
        monkeypatch.delitem(sys.modules, name, raising=False)
    import hvantk.algorithms.rerank.nulls  # noqa: F401
