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
    is reported because it came top. **Only this one is a multiplicity correction**: a
    per-axis null can sit near zero while the selected-maximum null does not, because the
    maximum over several candidate axes is stochastically larger than any single one of
    them -- which is why only the latter corrects best-of-N reporting.

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

The CV SEED is held fixed across permutations, not the partition: the stratified partition
is a function of ``(seed, labels)`` and is recomputed for every permutation exactly as it is
for the observed statistic, which is what makes the observed value and the draws
exchangeable. The null therefore already contains partition noise as one of its components;
a partition-variance sweep (``evaluator.seed_sweep``) must not be added to it as an
independent one, or that source of variation is counted twice.

What null hypothesis this tests: permuting the label tests the GLOBAL null that no feature,
baseline included, carries label information -- not the conditional null that THIS axis adds
nothing given the baseline already in the model. With an informative baseline, the
permutation null of the delta is wider than the delta's true sampling spread under the
conditional null, so the test is valid but conservative, increasingly so as the baseline
strengthens: failing to clear the selected-maximum null is weak evidence of no increment,
not proof of one.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from hvantk.algorithms.rerank.seeds import DEFAULT_SEED, rng_for

if TYPE_CHECKING:
    from hvantk.algorithms.rerank.leakage import LeakagePolicy
    from hvantk.algorithms.rerank.selection import SelectionPolicy


@dataclass(frozen=True)
class ControlSetting:
    """Every choice that changes the statistic a null is generated under. Hashable, compared
    by value: two nulls are comparable only if their settings are equal, and attaching a
    null built under one setting to deltas computed under another is a category error, not a
    warning.

    arm
        Provenance arm ("all" / "clean") -- which columns survived the circularity filter.
    leakage
        The :class:`~hvantk.algorithms.rerank.leakage.LeakagePolicy` in force, or ``None``
        if the leakage guard was off. Stored as the policy itself rather than a bool:
        ``q=0.05`` and ``q=0.2`` bar different columns and so change the statistic, but both
        would compare equal as ``leakage=True``.
    selection
        The :class:`~hvantk.algorithms.rerank.selection.SelectionPolicy` in force, or
        ``None`` if within-axis selection was off. Same reasoning as ``leakage``: selection
        on and selection off can both run under provenance arm "all", so the arm alone
        cannot distinguish them.
    baseline
        The exact baseline columns, normalised to ``tuple[str, ...]`` -- not the axis name.
        The circularity arms differ only in WHICH columns survive the provenance filter, so
        two runs can share an arm label and a baseline axis name and still be permuting
        different matrices.
    candidates
        The OFFERED axes and their columns, as ``((axis, (col, ...)), ...)`` in the order
        given. A selected-maximum null computed over a different candidate set is a
        different statistic even when the baseline and the winning axis are identical.
    folds
        The fold count the scorer ACTUALLY used -- not ``Config.folds``, which governs
        ``ReRanker.score`` only (see ``evaluator.ABLATION_FOLDS``).
    blocked
        Whether folds were group-blocked.

    Both ``leakage`` and ``selection`` are themselves frozen dataclasses with scalar fields,
    so an instance of this record stays hashable.
    """

    arm: str
    leakage: "LeakagePolicy | None"
    selection: "SelectionPolicy | None"
    baseline: tuple
    candidates: tuple
    folds: int
    blocked: bool

    def __post_init__(self) -> None:
        from hvantk.algorithms.rerank.leakage import LeakagePolicy
        from hvantk.algorithms.rerank.selection import SelectionPolicy

        if isinstance(self.baseline, str):
            raise TypeError(
                f"baseline must be a sequence of column names, not a str: {self.baseline!r}"
            )
        object.__setattr__(self, "baseline", tuple(str(c) for c in self.baseline))

        candidates = self.candidates
        if isinstance(candidates, str):
            raise TypeError(
                "candidates must be a mapping or a sequence of (axis, cols) pairs, not a "
                f"str: {candidates!r}"
            )
        if isinstance(candidates, Mapping):
            items = list(candidates.items())
        else:
            items = list(candidates)
        normalised = []
        for axis, cols in items:
            if isinstance(cols, str):
                raise TypeError(
                    f"candidates[{axis!r}] must be a sequence of column names, not a str: "
                    f"{cols!r}"
                )
            normalised.append((str(axis), tuple(str(c) for c in cols)))
        object.__setattr__(self, "candidates", tuple(normalised))

        if self.leakage is not None and not isinstance(self.leakage, LeakagePolicy):
            raise TypeError(
                f"leakage must be None or a LeakagePolicy; got {type(self.leakage).__name__}"
            )
        if self.selection is not None and not isinstance(self.selection, SelectionPolicy):
            raise TypeError(
                "selection must be None or a SelectionPolicy; got "
                f"{type(self.selection).__name__}"
            )

        if isinstance(self.folds, bool) or not isinstance(self.folds, int) or self.folds < 2:
            raise ValueError(f"folds must be an int >= 2; got {self.folds!r}")

    def describe(self) -> str:
        candidates = {axis: len(cols) for axis, cols in self.candidates}
        return (
            f"arm={self.arm!r} leakage={self.leakage!r} selection={self.selection!r} "
            f"folds={self.folds} blocked={self.blocked} baseline={list(self.baseline)!r} "
            f"candidates={candidates!r}"
        )


@dataclass(frozen=True)
class NullConfig:
    """How many permutations, and which slice of them this process computes.

    ``chunk``/``n_chunks`` carve a CONTIGUOUS slice of the permutation index rather than
    striding it, so chunks are disjoint, together cover ``[0, n_perm)``, and the seed of
    permutation ``i`` is ``seed + i`` regardless of the division. Validated here, at
    construction, so an invalid config can never be built and passed around before failing
    later inside ``span()``.
    """

    n_perm: int = 200
    chunk: int = 0
    n_chunks: int = 1
    seed: int = DEFAULT_SEED

    def __post_init__(self) -> None:
        if self.n_perm < 1:
            raise ValueError(f"n_perm must be >= 1; got {self.n_perm}")
        if self.n_chunks < 1:
            raise ValueError(f"n_chunks must be >= 1; got {self.n_chunks}")
        if not 0 <= self.chunk < self.n_chunks:
            raise ValueError(
                f"chunk must be in [0, n_chunks); got chunk={self.chunk}, "
                f"n_chunks={self.n_chunks}"
            )
        if self.n_chunks > self.n_perm:
            raise ValueError(
                f"n_chunks must be <= n_perm, or some chunks would be empty; got "
                f"n_chunks={self.n_chunks}, n_perm={self.n_perm}"
            )
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError(f"seed must be an int >= 0; got {self.seed!r}")

    def span(self) -> tuple:
        lo = (self.n_perm * self.chunk) // self.n_chunks
        hi = (self.n_perm * (self.chunk + 1)) // self.n_chunks
        return lo, hi


def oof_scorer(
    *, selector=None, groups=None, seed: int = DEFAULT_SEED, folds: int | None = None
) -> Callable:
    """The default scorer: ``evaluator._raw_oof`` with this run's selector, blocks, seed and
    fold count.

    Injectable because the null must be computed with the SAME estimator, selector and fold
    structure the observed delta was -- binding them together in one callable is what makes
    that hard to get wrong. A caller wrapping ``_raw_oof`` to add caching, or to swap in a
    different estimator entirely, can still hand the result to ``permutation_deltas`` as
    ``scorer=``, provided it keeps the same ``(matrix, cols, y) -> scores`` shape.
    """
    from hvantk.algorithms.rerank.evaluator import ABLATION_FOLDS, _raw_oof

    resolved_folds = ABLATION_FOLDS if folds is None else folds

    def _score(matrix, cols, y):
        return _raw_oof(
            matrix, list(cols), y, selector=selector, groups=groups, seed=seed,
            folds=resolved_folds,
        )

    return _score


_TIE_TOL = 1e-12
"""Absolute tolerance for counting a null draw as at least as extreme as the observed delta.

sklearn computes AUC by trapezoidal integration, so deltas equal in exact arithmetic often
differ in the last few bits (``0.3 - 0.1 < 0.2``), which undercounts ties and makes the
p-value too small. AUC deltas sit on a grid of spacing ``1 / (2 * n_pos * n_neg)`` --
``>= ~1e-8`` for any gene-level universe -- far above both this tolerance and float noise
(``~1e-16``). Absolute rather than SciPy's relative ``gamma``: these deltas cluster near
zero, where a relative tolerance vanishes.
"""


def p_value(null, observed: float) -> float:
    """``(1 + #{null >= obs}) / (1 + n)``, where ``n`` counts only FINITE draws.

    The ``+1`` in both places is not a smoothing fudge: the observed statistic is itself one
    draw from the null under the hypothesis being tested, so a permutation p-value that can
    reach 0 is reporting a certainty ``n`` permutations do not buy (Phipson & Smyth, SAGMB
    2010).

    ``observed`` must be the unrounded point delta computed with the SAME scorer as the
    null -- never a rounded or bootstrap-median summary; those are a different statistic and
    do not share this null's distribution. It must also be finite: a NaN observed value
    compares False against every draw, which would otherwise silently return the smallest p
    the null can produce instead of signalling that there was nothing to test.

    NaN (and +/-inf) draws are DROPPED, not counted as zero: an axis wholly contained in the
    baseline has no delta to contribute, and a delta of exactly zero is a measurement rather
    than an absence.
    """
    observed = float(observed)
    if not np.isfinite(observed):
        raise ValueError(f"observed delta must be finite; got {observed!r}")
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0:
        raise ValueError("no finite draws in the null; nothing to compare against")
    return float((1 + int((null >= observed - _TIE_TOL).sum())) / (1 + null.size))


def axis_deltas(matrix, baseline, axes, y, *, scorer) -> tuple:
    """``(baseline AUC, {axis: delta-AUC of baseline+axis over baseline})`` for ONE label
    vector.

    The observed statistic is this function called on the real labels; each null draw is
    this function called on a permuted copy of the labels -- the same function, which is
    what makes the observed value and the draws comparable. An axis whose columns all lie
    inside the baseline gets NaN, not 0.0: a delta of exactly zero is a measurement, and an
    axis with nothing left to add has none to report.
    """
    from sklearn.metrics import roc_auc_score

    baseline = list(baseline)
    in_baseline = set(baseline)
    base_auc = float(roc_auc_score(y, scorer(matrix, baseline, y)))
    deltas = {}
    for name, cols in axes.items():
        extra = [c for c in cols if c not in in_baseline]
        if not extra:
            deltas[name] = float("nan")
            continue
        a1 = float(roc_auc_score(y, scorer(matrix, baseline + extra, y)))
        deltas[name] = a1 - base_auc
    return base_auc, deltas


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
    permutation, not per row, so a reader can see that it MOVED -- the tell-tale that the
    baseline was refit on the permuted labels rather than reused from the observed fit,
    which would otherwise measure the permutation instead of the axis.

    ``scorer=None`` means the plain default: ``_raw_oof`` with no selector, no blocks,
    ``DEFAULT_SEED`` and ``ABLATION_FOLDS`` folds. If the observed deltas were computed with
    a selector, with blocked folds, with another seed or with another fold count, the caller
    must pass ``oof_scorer(...)`` built from those SAME inputs -- otherwise the null bounds
    a different model from the one the observed value came from.
    """
    y = np.asarray(y)
    baseline = list(baseline)
    score = scorer if scorer is not None else oof_scorer()
    lo, hi = config.span()

    rows = []
    for i in range(lo, hi):
        # seed + i, so this permutation is the same one whichever chunk computes it.
        yp = rng_for(config.seed, i).permutation(y)
        # axis_deltas refits the baseline HERE, on the permuted label. Reusing an unpermuted
        # baseline AUC would make every delta a measure of the permutation rather than of
        # the axis.
        base_auc, deltas = axis_deltas(matrix, baseline, axes, yp, scorer=score)
        for name, delta in deltas.items():
            rows.append({"perm": i, "axis": name, "delta": delta, "base_auc": base_auc})
    return pd.DataFrame(rows, columns=["perm", "axis", "delta", "base_auc"])
