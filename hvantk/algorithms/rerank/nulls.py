"""Permutation nulls, and the multiplicity correction the control stack was missing.

``provenance.py`` asks what a predictor was TRAINED on. ``leakage.py`` (#244) asks which
units it was RUN on. Neither asks the question a reader actually has: *the authors offered
N axes and reported the best one -- could that have arisen by chance?* Before this module
`grep -rn "permut\\|selected_max" hvantk/algorithms/rerank/` returned nothing, so a user of
``hvantk rerank`` got a delta-AUC with no way to ask whether it survives selection.

Two distributions, because they answer different questions:

per-axis
    Permute the label, refit baseline and baseline+THIS axis, record the delta. Correct
    when the axis was pre-specified. Its spread scales with axis WIDTH, so a 3-column axis
    and a 25-column axis are not held to one threshold.
selected-maximum
    Permute, refit, add EACH offered axis, record the LARGEST delta. Correct when the axis
    is reported because it came top. **Only this one is a multiplicity correction.** On the
    published CHD matrix the per-axis null sat near zero while the selected-max null had
    median +0.024, which reversed a headline result
    (``analysis/rerank-homogenised/perm_null_arm.py:2-13``,
    ``analysis/crossconfig-multiplicity/run.py:9-14``).

Three contract points that are easy to get wrong, and are enforced here rather than left to
the caller:

1. **Each permutation refits the baseline too.** Scoring a permuted-label augmented model
   against an unpermuted baseline measures the permutation, not the axis.
2. **``p = (1 + #{null >= obs}) / (1 + n_perm)``.** The analysis drivers used the plain
   ``(null >= obs).mean()`` and could therefore report ``p = 0.000``; a finite permutation
   set can never license that (Phipson & Smyth 2010).
3. **A null belongs to ONE control setting.** Leakage on/off, the provenance arm, the exact
   baseline columns, the fold count and whether folds were blocked all change what is being
   permuted. ``NullDistribution`` records the setting it was built under and raises rather
   than answer a question about a different one.

The CV PARTITION is held fixed across permutations on purpose: only the labels move, so
the distribution isolates label information rather than mixing in partition noise. Partition
variance is a separate component and is measured by ``evaluator.seed_sweep``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from hvantk.algorithms.rerank.seeds import DEFAULT_SEED, rng_for


@dataclass(frozen=True)
class ControlSetting:
    """What a null was generated under. Hashable, compared by value.

    ``baseline`` is the exact tuple of baseline columns, not the axis name: the circularity
    arms differ only in WHICH columns survive the provenance filter, so two runs can share
    an arm label and a baseline axis name and still be permuting different matrices.
    """

    arm: str
    leakage: bool
    baseline: tuple
    folds: int
    blocked: bool

    def describe(self) -> str:
        return (
            f"arm={self.arm!r} leakage={self.leakage} folds={self.folds} "
            f"blocked={self.blocked} baseline={list(self.baseline)!r}"
        )


@dataclass(frozen=True)
class NullConfig:
    """How many permutations, and which slice of them this process computes.

    ``chunk``/``n_chunks`` carve a CONTIGUOUS slice of the permutation index rather than
    striding it, so chunks are disjoint, together cover ``[0, n_perm)``, and the seed of
    permutation ``i`` is ``seed + i`` regardless of the division
    (``analysis/rerank-homogenised/perm_null_arm.py:77-80``).
    """

    n_perm: int = 200
    chunk: int = 0
    n_chunks: int = 1
    seed: int = DEFAULT_SEED

    def span(self) -> tuple:
        if self.n_perm < 1:
            raise ValueError(f"n_perm must be >= 1; got {self.n_perm}")
        if self.n_chunks < 1:
            raise ValueError(f"n_chunks must be >= 1; got {self.n_chunks}")
        if not 0 <= self.chunk < self.n_chunks:
            raise ValueError(
                f"chunk must be in [0, n_chunks); got chunk={self.chunk}, "
                f"n_chunks={self.n_chunks}"
            )
        lo = (self.n_perm * self.chunk) // self.n_chunks
        hi = (self.n_perm * (self.chunk + 1)) // self.n_chunks
        return lo, hi


def oof_scorer(*, selector=None, groups=None, seed: int = DEFAULT_SEED) -> Callable:
    """The default scorer: ``evaluator._raw_oof`` with this run's selector, blocks and seed.

    Injectable because the null must be computed with the SAME estimator, selector and fold
    structure the observed delta was -- binding them together in one callable is what makes
    that hard to get wrong -- and because the analysis code does the same
    (``analysis/rerank-homogenised/run_arm.py:197`` wraps ``_raw_oof`` to add grouped folds
    and caching, then hands the wrapper to the permutation driver at ``perm_null_arm.py:88``).
    """
    from hvantk.algorithms.rerank.evaluator import _raw_oof

    def _score(matrix, cols, y):
        return _raw_oof(matrix, list(cols), y, selector=selector, groups=groups, seed=seed)

    return _score


def p_value(null, observed: float) -> float:
    """``(1 + #{null >= obs}) / (1 + n_perm)``.

    The ``+1`` in both places is not a smoothing fudge: the observed statistic is itself one
    draw from the null under the hypothesis being tested, so a permutation p-value that can
    reach 0 is reporting a certainty ``n_perm`` permutations do not buy (Phipson & Smyth,
    SAGMB 2010). ``analysis/crossconfig-multiplicity/run.py:152-153`` and
    ``analysis/group-blocked-cv/run.py:260`` used the plain mean; this is the corrected form
    and the difference is exactly one reported ``0.000``.

    NaN draws are DROPPED, not counted as zero: an axis wholly contained in the baseline has
    no delta to contribute, and a delta of exactly zero is a measurement rather than an
    absence.
    """
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0:
        raise ValueError("no finite draws in the null; nothing to compare against")
    return float((1 + int((null >= observed).sum())) / (1 + null.size))


def permutation_deltas(
    matrix,
    baseline: Sequence[str],
    axes: Mapping[str, Sequence[str]],
    y,
    *,
    config: NullConfig,
    scorer: Callable | None = None,
) -> pd.DataFrame:
    """One row per (permutation, axis): the delta-AUC of baseline+axis over baseline alone.

    Returns columns ``perm, axis, delta, base_auc``. ``base_auc`` is recorded per
    permutation, not per row, so a reader can see that it MOVED -- the tell-tale of the
    mistake this loop exists to avoid.
    """
    from sklearn.metrics import roc_auc_score

    y = np.asarray(y)
    baseline = list(baseline)
    in_baseline = set(baseline)
    score = scorer if scorer is not None else oof_scorer()
    lo, hi = config.span()

    rows = []
    for i in range(lo, hi):
        # seed + i, so this permutation is the same one whichever chunk computes it.
        yp = rng_for(config.seed, i).permutation(y)
        # The baseline is refit HERE, on the permuted label. Reusing an unpermuted baseline
        # AUC would make every delta a measure of the permutation rather than of the axis
        # (analysis/rerank-homogenised/perm_null_arm.py:88).
        a0 = float(roc_auc_score(yp, score(matrix, baseline, yp)))
        for name, cols in axes.items():
            extra = [c for c in cols if c not in in_baseline]
            if not extra:
                rows.append(
                    {"perm": i, "axis": name, "delta": float("nan"), "base_auc": a0}
                )
                continue
            a1 = float(roc_auc_score(yp, score(matrix, baseline + extra, yp)))
            rows.append({"perm": i, "axis": name, "delta": a1 - a0, "base_auc": a0})
    return pd.DataFrame(rows, columns=["perm", "axis", "delta", "base_auc"])
