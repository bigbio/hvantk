"""Pass B of #247: a gene family must not straddle a fold.

Paralogues share sequence, constraint, expression pattern and disease status, so a family
split across train and test lets the model recognise a relative rather than generalise --
the standard robustness objection to naive cross-validation in statistical genetics, and
nothing in the library addressed it before this.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hvantk.algorithms.rerank.evaluator import Evaluator, _raw_oof
from hvantk.algorithms.rerank.reranker import ReRanker, _cv
from hvantk.tests.rerank._synth import planted_signal


def _blocked_fixture(n=180, family_size=6, seed=5):
    matrix, y, baseline, axes = planted_signal(n=n, n_noise=1, seed=seed)
    blocks = np.arange(n) // family_size
    return matrix, y, baseline, axes, blocks


def _fold_assignment(matrix, cols, y, groups, folds=5, seed=0):
    cv = _cv(folds, seed, groups)
    assign = np.full(len(y), -1)
    splitter = cv.split(matrix[cols].values, y, groups) if groups is not None else cv.split(
        matrix[cols].values, y
    )
    for k, (_, te) in enumerate(splitter):
        assign[te] = k
    return assign


def test_every_block_sits_in_a_single_test_fold():
    matrix, y, baseline, _, blocks = _blocked_fixture()
    assign = _fold_assignment(matrix, baseline, y, blocks)
    per_block = pd.DataFrame({"block": blocks, "fold": assign}).groupby("block").fold.nunique()
    assert (per_block == 1).all(), per_block[per_block > 1]

    # With every group a singleton, a correctly stratified grouped splitter must balance
    # positives across test folds the same way StratifiedKFold does (differ by at most 1).
    # sklearn's own `StratifiedGroupKFold(shuffle=True)` decides each group's fold from
    # another group's shuffled class counts, not its own, so it does not -- this is the
    # defect the seeded relabelling fixes.
    singleton_groups = np.arange(len(y))
    cv = _cv(5, 0, singleton_groups)
    pos_per_fold = [
        int(y[test_idx].sum())
        for _, test_idx in cv.split(matrix[baseline].values, y, singleton_groups)
    ]
    assert max(pos_per_fold) - min(pos_per_fold) <= 1, pos_per_fold


def test_raw_oof_accepts_groups_and_changes_the_answer():
    """If groups were accepted and ignored, every blocked run would silently be an
    unblocked one -- the exact failure a `groups=` parameter is supposed to close."""
    matrix, y, baseline, axes, blocks = _blocked_fixture()
    cols = baseline + axes["axis0"]
    plain = _raw_oof(matrix, cols, y)
    grouped = _raw_oof(matrix, cols, y, groups=blocks)
    assert grouped.shape == plain.shape
    assert not np.allclose(plain, grouped), "groups= had no effect"
    assert np.isfinite(grouped).all()


def test_raw_oof_honours_groups_on_the_per_fold_selector_path_too():
    """The ablation path runs a selector per fold; it must use the SAME grouped split, or
    the headline is blocked and the deltas are not."""
    matrix, y, baseline, axes, blocks = _blocked_fixture()
    cols = baseline + axes["axis0"]
    assert isinstance(matrix.index, pd.RangeIndex)
    seen = []

    def spy(X_tr, y_tr, columns):
        seen.append(X_tr.index.to_numpy())
        return list(columns)

    _raw_oof(matrix, cols, y, selector=spy, groups=blocks)
    assert len(seen) == 5
    # Each training slice is a union of whole blocks: no block sits on both sides of the
    # train/test boundary within a fold.
    for train in seen:
        test = set(range(len(y))) - set(train)
        assert set(blocks[list(train)]).isdisjoint(blocks[list(test)])


def test_reranker_score_accepts_groups():
    matrix, y, baseline, axes, blocks = _blocked_fixture()
    cols = baseline + axes["axis0"]
    p = ReRanker(folds=5).score(matrix, cols, y, groups=blocks)
    assert p.shape == (len(y),) and np.isfinite(p).all()
    assert not np.allclose(p, ReRanker(folds=5).score(matrix, cols, y))


def test_evaluator_threads_groups_into_the_ablation(monkeypatch):
    """A list-inequality check on the resulting ablation AUCs would still pass if only ONE
    of `Evaluator.evaluate`'s two `_raw_oof` calls received `groups` -- verified with a
    mutant. A spy on `_raw_oof` is the only way to confirm both calls forward the blocks,
    not just that the answer moved."""
    import hvantk.algorithms.rerank.evaluator as evaluator_mod

    matrix, y, baseline, axes, blocks = _blocked_fixture()
    groups_map = {"base": baseline, "axis0": axes["axis0"]}
    scores = _raw_oof(matrix, baseline, y, groups=blocks)

    real_raw_oof = evaluator_mod._raw_oof
    calls = []

    def spy(*args, **kwargs):
        calls.append(kwargs.get("groups"))
        return real_raw_oof(*args, **kwargs)

    monkeypatch.setattr(evaluator_mod, "_raw_oof", spy)
    Evaluator().evaluate(
        matrix, baseline + axes["axis0"], y, scores, groups_map, "base", groups=blocks
    )
    assert len(calls) >= 2
    assert all(np.array_equal(g, blocks) for g in calls)


def test_a_degenerate_blocked_fold_is_refused_with_a_clear_message():
    """Guard the guard: StratifiedGroupKFold cannot split a block, so an oversized block
    unbalances the folds rather than leaking -- but an unbalanced fold can still leave a
    training slice with no example of a class, or a test slice with none at all. Both must
    be refused with a message naming the problem, not left to a silent column of 0.0
    predictions or an opaque estimator error. `gene_blocks`'s dominant-block ceiling
    (``max_block_frac``) is a separate, additional guard for the SILENT unbalancing this
    one does not cover."""
    matrix, y, baseline, _, _ = _blocked_fixture(n=60)
    one_block = np.zeros(60, dtype=int)  # one block for everything
    with pytest.raises(ValueError, match="lower the fold count"):
        _raw_oof(matrix, baseline, y, groups=one_block)

    # All positives confined to ONE block, every other gene a singleton, on the DEFAULT
    # ReRanker.score path -- exactly where a degenerate fold used to collapse into a
    # silent column of 0.0 predictions (a classifier fit on one class alone predicts it
    # with probability 1) instead of raising.
    matrix2, _, baseline2, _ = planted_signal(n=120, n_noise=1, seed=9)
    y2 = np.zeros(120, dtype=int)
    y2[:6] = 1
    blocks2 = np.arange(120)
    blocks2[:6] = 0  # the 6 positives share block 0; genes 6..119 are singletons
    with pytest.raises(ValueError, match="lower the fold count"):
        ReRanker(folds=5).score(matrix2, baseline2, y2, groups=blocks2)
