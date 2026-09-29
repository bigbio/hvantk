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
from hvantk.algorithms.rerank.reranker import ReRanker
from hvantk.tests.rerank._synth import planted_signal


def _blocked_fixture(n=180, family_size=6, seed=5):
    matrix, y, baseline, axes = planted_signal(n=n, n_noise=1, seed=seed)
    blocks = np.arange(n) // family_size
    return matrix, y, baseline, axes, blocks


def _fold_assignment(matrix, cols, y, groups, folds=5, seed=0):
    from hvantk.algorithms.rerank.reranker import _cv

    cv = _cv(folds, seed, groups)
    assign = np.full(len(y), -1)
    splitter = cv.split(matrix[cols].values, y, groups) if groups is not None else cv.split(
        matrix[cols].values, y
    )
    for k, (_, te) in enumerate(splitter):
        assign[te] = k
    return assign


def test_no_fold_contains_two_members_of_one_block():
    matrix, y, baseline, _, blocks = _blocked_fixture()
    assign = _fold_assignment(matrix, baseline, y, blocks)
    per_block = pd.DataFrame({"block": blocks, "fold": assign}).groupby("block").fold.nunique()
    assert (per_block == 1).all(), per_block[per_block > 1]


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


def test_raw_oof_without_groups_is_byte_identical_to_before():
    matrix, y, baseline, _, _ = _blocked_fixture()
    assert np.allclose(_raw_oof(matrix, baseline, y), _raw_oof(matrix, baseline, y, groups=None))


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


def test_evaluator_threads_groups_into_the_ablation():
    matrix, y, baseline, axes, blocks = _blocked_fixture()
    groups_map = {"base": baseline, "axis0": axes["axis0"]}
    scores = _raw_oof(matrix, baseline, y, groups=blocks)
    plain = Evaluator().evaluate(matrix, baseline + axes["axis0"], y, scores, groups_map, "base")
    blocked = Evaluator().evaluate(
        matrix, baseline + axes["axis0"], y, scores, groups_map, "base", groups=blocks
    )
    assert plain.ablation.auc.tolist() != blocked.ablation.auc.tolist()


def test_a_block_larger_than_a_fold_is_a_sklearn_error_not_a_silent_split():
    """Guard the guard: StratifiedGroupKFold cannot split a block, so an oversized block
    unbalances the folds rather than leaking. Task 5's ceiling exists because that
    unbalancing is SILENT -- here we only record that sklearn does not leak."""
    matrix, y, baseline, _, _ = _blocked_fixture(n=60)
    blocks = np.zeros(60, dtype=int)  # one block for everything
    with pytest.raises(ValueError):
        _raw_oof(matrix, baseline, y, groups=blocks)
