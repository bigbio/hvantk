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
2. **``p = (1 + #{null >= obs}) / (1 + n_perm)``.** A naive ``p = (null >= obs).mean()``
   reports ``p = 0.000`` for an observed value no permutation reached; a finite permutation
   set can never license that (Phipson & Smyth 2010).
3. **A null belongs to ONE control setting.** The leakage and selection policies, the
   provenance arm, the exact baseline columns, the offered candidate axes and their
   columns, the fold count the scorer actually used, and which blocking (if any) the folds
   used all change what is being permuted. ``NullDistribution`` records the setting it was
   built under and raises rather than answer a question about a different one.

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
strengthens: failing to clear the selected-maximum null is weak evidence that the axis adds
nothing, not proof that it does not -- PROVIDED labels are exchangeable across genes. The
cause of the caveat below is that assumption breaking, not a bookkeeping mismatch with the
folds: paralogues share labels, so genes are not exchangeable with one another when a
universe has them, and a permutation that treats every gene as interchangeable with every
other then understates the null's spread and is ANTI-conservative rather than merely
conservative. When the folds are blocked, the scheme that restores exchangeability is the
block permutation described above (whole blocks permuted among blocks of the same size); a
GLOBAL permutation of the same labels does not restore it, whether or not the folds
themselves happen to be blocked. This is why the DEFAULT (unblocked) null carries the same
risk: its "valid but conservative" statement assumes labels exchangeable across genes, and
the clustering is a property of the DATA -- omitting ``Config.blocks``/``--blocks`` does not
remove it, it only removes the correction for it. Pass blocks (``--blocks``) whenever labels
may cluster by family.
"""

from __future__ import annotations

import logging
import numbers
from dataclasses import dataclass, fields
from typing import TYPE_CHECKING, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from hvantk.algorithms.rerank.seeds import DEFAULT_SEED, rng_for

if TYPE_CHECKING:
    from hvantk.algorithms.rerank.leakage import LeakagePolicy
    from hvantk.algorithms.rerank.selection import SelectionPolicy

logger = logging.getLogger(__name__)


def _coerce_int(value, name: str) -> int:
    """Accept any integral scalar (a bare ``int``, ``np.int64``, ...) but never a ``bool``.

    ``NullConfig(n_perm=1e3)`` must fail HERE, at construction, rather than deep inside
    ``permutation_deltas`` -- a float ``n_perm`` breaks ``span()``'s slicing arithmetic --
    and ``NullConfig(seed=np.int64(7))`` must not fail at all: numpy's own integer types
    are exactly as valid a seed as a bare ``int``.
    """
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an int; got {value!r}")
    return int(value)


def _reject_set(value, what: str) -> None:
    """Refuse a bare ``set``/``frozenset`` wherever column ORDER is part of the setting.

    A set becomes a tuple in hash-randomised order, so two chunk processes constructing
    "the same" setting from a set-valued baseline or candidate column list could disagree
    on ``ControlSetting`` equality for a reason that has nothing to do with the run itself.
    """
    if isinstance(value, (set, frozenset)):
        raise TypeError(
            f"{what} must not be a set: column order is part of the setting; pass a list "
            "or tuple"
        )


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
    seed
        The CV seed the scorer's folds were drawn with (``Config.seed``, as handed to
        ``oof_scorer``). The partition is a function of ``(seed, labels)``, so a null and
        a delta computed under different seeds are different draws of the partition
        component and are not the same statistic: without this field two chunks of one
        ``NullConfig`` scored under CV seeds 0 and 1 merged silently when no ``observed``
        rode with them, and a seed-0 null answered ``p_selected_max`` for a seed-1 delta
        (observed deltas differed by up to 0.049). Not ``NullConfig.seed``, which seeds
        the permutations and is recorded separately as ``perm_seed``.
    block_digest
        A digest of the exact block labels the folds were built from (see
        :func:`~hvantk.algorithms.rerank.blocks.block_digest`), or ``None`` when folds were
        unblocked. Recording WHICH blocking, not only whether one was used, matters because
        nulls computed under two different block tables -- another HGNC release, another
        grouping -- would otherwise compare equal. A non-``None`` digest also means the
        null's labels were permuted BY those same blocks, not globally: the engine builds
        both the scorer's folds and the permutation from the one ``blocks`` array, so this
        digest speaks for both.

    Both ``leakage`` and ``selection`` are themselves frozen dataclasses with scalar fields,
    so an instance of this record stays hashable.
    """

    arm: str
    leakage: "LeakagePolicy | None"
    selection: "SelectionPolicy | None"
    baseline: tuple
    candidates: tuple
    folds: int
    seed: int
    block_digest: "str | None"

    def __post_init__(self) -> None:
        from hvantk.algorithms.rerank.leakage import LeakagePolicy
        from hvantk.algorithms.rerank.selection import SelectionPolicy

        if isinstance(self.baseline, str):
            raise TypeError(
                f"baseline must be a sequence of column names, not a str: {self.baseline!r}"
            )
        _reject_set(self.baseline, "baseline")
        object.__setattr__(self, "baseline", tuple(str(c) for c in self.baseline))

        candidates = self.candidates
        if isinstance(candidates, str):
            raise TypeError(
                "candidates must be a mapping or a sequence of (axis, cols) pairs, not a "
                f"str: {candidates!r}"
            )
        _reject_set(candidates, "candidates")
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
            _reject_set(cols, f"candidates[{axis!r}]")
            normalised.append((str(axis), tuple(str(c) for c in cols)))
        object.__setattr__(self, "candidates", tuple(normalised))

        if self.leakage is not None and not isinstance(self.leakage, LeakagePolicy):
            raise TypeError(
                f"leakage must be None or a LeakagePolicy; got {type(self.leakage).__name__}"
            )
        if self.selection is not None and not isinstance(
            self.selection, SelectionPolicy
        ):
            raise TypeError(
                "selection must be None or a SelectionPolicy; got "
                f"{type(self.selection).__name__}"
            )

        folds = _coerce_int(self.folds, "folds")
        object.__setattr__(self, "folds", folds)
        if folds < 2:
            raise ValueError(f"folds must be an int >= 2; got {folds!r}")

        seed = _coerce_int(self.seed, "seed")
        object.__setattr__(self, "seed", seed)
        if seed < 0:
            raise ValueError(f"seed must be an int >= 0; got {seed!r}")

        if self.block_digest is not None:
            if not isinstance(self.block_digest, str):
                raise TypeError(
                    f"block_digest must be None or a str; got "
                    f"{type(self.block_digest).__name__}"
                )
            if not self.block_digest.strip():
                raise ValueError("block_digest must be None or a non-empty str")

    def describe(self) -> str:
        candidates = {axis: len(cols) for axis, cols in self.candidates}
        digest = None if self.block_digest is None else self.block_digest[:12]
        return (
            f"arm={self.arm!r} leakage={self.leakage!r} selection={self.selection!r} "
            f"folds={self.folds} seed={self.seed} block_digest={digest!r} "
            f"baseline={list(self.baseline)!r} candidates={candidates!r}"
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
        n_perm = _coerce_int(self.n_perm, "n_perm")
        chunk = _coerce_int(self.chunk, "chunk")
        n_chunks = _coerce_int(self.n_chunks, "n_chunks")
        seed = _coerce_int(self.seed, "seed")
        object.__setattr__(self, "n_perm", n_perm)
        object.__setattr__(self, "chunk", chunk)
        object.__setattr__(self, "n_chunks", n_chunks)
        object.__setattr__(self, "seed", seed)

        if n_perm < 1:
            raise ValueError(f"n_perm must be >= 1; got {n_perm}")
        if n_chunks < 1:
            raise ValueError(f"n_chunks must be >= 1; got {n_chunks}")
        if not 0 <= chunk < n_chunks:
            raise ValueError(
                f"chunk must be in [0, n_chunks); got chunk={chunk}, "
                f"n_chunks={n_chunks}"
            )
        if n_chunks > n_perm:
            raise ValueError(
                f"n_chunks must be <= n_perm, or some chunks would be empty; got "
                f"n_chunks={n_chunks}, n_perm={n_perm}"
            )
        if seed < 0:
            raise ValueError(f"seed must be an int >= 0; got {seed!r}")

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

    The returned callable carries the ``groups`` it was built with as a ``.groups``
    attribute. This is what lets ``permutation_deltas`` notice a scorer built with blocked
    folds whose caller forgot to also pass ``blocks=`` -- without it, that mistake would
    silently fall back to a global label permutation while the scorer itself kept scoring
    blocked folds, understating the null's spread for exactly the reason blocking exists.
    """
    from hvantk.algorithms.rerank.evaluator import ABLATION_FOLDS, _raw_oof

    resolved_folds = ABLATION_FOLDS if folds is None else folds

    def _score(matrix, cols, y):
        return _raw_oof(
            matrix,
            list(cols),
            y,
            selector=selector,
            groups=groups,
            seed=seed,
            folds=resolved_folds,
        )

    _score.groups = groups
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


def _block_structure(blocks) -> dict:
    """Per block size, the member-index arrays in canonical order: blocks in order of FIRST
    APPEARANCE, members in ascending gene position.

    Precomputed ONCE per :func:`permutation_deltas` call and reused for every permutation --
    rebuilding it per draw would be pure waste, since it depends only on ``blocks``, never on
    the (permuted) labels. ``blocks`` is not assumed to be the engine's own 0..B-1 codes:
    :func:`permutation_deltas` is public, so any hashable label array is grouped with
    ``pandas.factorize`` -- the same canonicalisation :func:`~hvantk.algorithms.rerank.
    blocks.block_digest` uses, and a label array no longer needs to be sortable.

    First-appearance order, not ``np.unique``'s sorted-by-value order: two arrays that
    encode the SAME partition (which positions share a block) but assign different raw
    label values to each group would otherwise be visited in a different order and so draw
    different permutations for the same ``(seed, i)``, although they share a
    ``block_digest`` and so a ``ControlSetting``. Factorizing by first appearance makes the
    visiting order a function of the partition alone, matching ``block_digest``.
    """
    blocks = np.asarray(blocks)
    codes = pd.factorize(blocks, use_na_sentinel=False)[0]
    n_labels = int(codes.max()) + 1 if codes.size else 0
    members: list = [[] for _ in range(n_labels)]
    for position, code in enumerate(codes):
        members[code].append(position)
    structure: dict = {}
    for member_positions in members:
        arr = np.asarray(member_positions, dtype=int)
        structure.setdefault(arr.size, []).append(arr)
    return structure


def _permute_labels(y, rng, structure=None):
    """Draw ONE permutation of ``y``: global when ``structure`` is ``None``, else a block
    permutation over the exchangeability blocks ``structure`` describes.

    ``structure is None`` reduces to exactly ``rng.permutation(y)`` -- the whole call, no
    extra RNG draws either side of it -- so an unblocked null consumes its random-number
    generator in the identical sequence as before this function existed, and every unblocked
    null (and every test of one) is byte-identical.

    With a ``structure``, permutes whole blocks among blocks of the SAME size, then shuffles
    each block's own members (Winkler et al., NeuroImage 2015, "multi-level block
    permutation"): for each block size in ascending order, one draw
    ``order = rng.permutation(n_blocks_of_that_size)`` decides that source block ``k``'s
    label vector is placed at destination block ``order[k]``, and one further draw per block
    shuffles that vector before it is placed. Singletons (size 1) therefore permute freely
    among singletons, and a size-1 "shuffle" is a no-op, exactly as it should be. The whole
    draw is a function of ``(seed, i)`` alone: the size groups are visited in a fixed
    (ascending) order and nothing here reads global state, so the same ``rng`` in the same
    starting state always produces the same permutation.

    Power caveat: a block whose size occurs only ONCE in ``structure`` has no sibling block
    of the same size to swap with, so ``order = rng.permutation(1)`` is a no-op for it; if
    that block is also label-homogeneous (every member shares one label), the within-block
    shuffle is a no-op too, since shuffling identical values changes nothing. Such a block's
    labels never move, on any draw -- it contributes zero spread to the null, silently
    reducing its effective degrees of freedom below ``n_perm``.
    """
    if structure is None:
        return rng.permutation(y)
    y = np.asarray(y)
    out = np.empty_like(y)
    for size in sorted(structure):
        blocks_of_size = structure[size]
        order = rng.permutation(len(blocks_of_size))
        for source_k, source_members in enumerate(blocks_of_size):
            dest_members = blocks_of_size[order[source_k]]
            out[dest_members] = rng.permutation(y[source_members])
    return out


def permutation_deltas(
    matrix,
    baseline: Sequence[str],
    axes: Mapping[str, Sequence[str]],
    y,
    *,
    config: NullConfig,
    scorer: Callable | None = None,
    blocks=None,
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

    ``blocks=None`` permutes the label globally, which assumes every unit is exchangeable
    with every other. When the scorer's folds are blocked, paralogues share labels and that
    assumption is false, so pass the SAME block labels the scorer used here too: whole blocks
    are then permuted among blocks of the same size, and shuffled within each block (see
    :func:`_permute_labels`), which is the null a blocked scorer needs -- a global permutation
    understates its spread and is anti-conservative. The engine does this for you. Note that a
    permuted labelling can still produce a degenerate blocked fold (e.g. a class confined to
    too few blocks); the scorer's own grouped-split guard raises in that case rather than
    returning a partial null.

    If ``scorer`` declares the blocks it was built with (as ``oof_scorer`` does, via a
    ``.groups`` attribute on the callable it returns), those must describe the SAME
    PARTITION as ``blocks`` (compared by :func:`~hvantk.algorithms.rerank.blocks.
    block_digest`, not by raw label equality -- two arrays that group the same positions
    together but spell the groups differently are the same partition) -- otherwise a caller
    who built a blocked scorer and forgot ``blocks=`` here would silently get the
    anti-conservative global permutation while the scorer kept scoring blocked folds. A
    custom scorer with no ``.groups`` attribute is not checked.
    """
    y = np.asarray(y)
    baseline = list(baseline)
    score = scorer if scorer is not None else oof_scorer()
    # Length checked FIRST, before the declared-blocks comparison below: a wrong-length
    # `blocks` must report the length mismatch, not a confusing "different partition" (or
    # "no blocks passed") message that has nothing to do with the actual problem.
    if blocks is not None:
        blocks = np.asarray(blocks)
        if len(blocks) != len(y):
            raise ValueError(
                f"blocks must have one entry per label; got {len(blocks)} for {len(y)} "
                "labels"
            )
    declared = getattr(score, "groups", None)
    if declared is not None:
        from hvantk.algorithms.rerank.blocks import block_digest

        if blocks is None:
            raise ValueError(
                "the scorer's folds are blocked but no blocks= was passed; pass the same "
                "blocks the scorer was built with"
            )
        if block_digest(blocks) != block_digest(declared):
            raise ValueError(
                "blocks= describes a different partition from the scorer's groups; pass "
                "the same blocks the scorer was built with"
            )
    structure = None if blocks is None else _block_structure(blocks)
    lo, hi = config.span()
    chunk_size = hi - lo
    # ~every 10% of the chunk; for a chunk under 10 draws that floors to 0, so `max(1, ...)`
    # falls back to logging every permutation instead of going silent for the whole chunk.
    log_every = max(1, chunk_size // 10)

    rows = []
    every_draw_unchanged = True
    for i in range(lo, hi):
        # seed + i, so this permutation is the same one whichever chunk computes it.
        yp = _permute_labels(y, rng_for(config.seed, i), structure)
        every_draw_unchanged = every_draw_unchanged and np.array_equal(yp, y)
        # axis_deltas refits the baseline HERE, on the permuted label. Reusing an unpermuted
        # baseline AUC would make every delta a measure of the permutation rather than of
        # the axis.
        base_auc, deltas = axis_deltas(matrix, baseline, axes, yp, scorer=score)
        for name, delta in deltas.items():
            rows.append({"perm": i, "axis": name, "delta": delta, "base_auc": base_auc})
        done = i - lo + 1
        if done % log_every == 0 or done == chunk_size:
            logger.info(
                "permutation null: %d/%d draws (chunk %d/%d)",
                done,
                chunk_size,
                config.chunk,
                config.n_chunks,
            )
    if chunk_size > 0 and every_draw_unchanged:
        # Every block was either a singleton or a size-occurs-once, label-homogeneous block
        # (see _permute_labels's power caveat): nothing in THIS CHUNK ever moved. That is a
        # statement about the chunk, not about the whole null -- a chunk of one draw can land
        # on the identity permutation by chance, so this must not be read as "the null has no
        # spread" until the MERGED null (NullDistribution.merge) confirms it.
        logger.warning(
            "every permutation in this chunk (%d draw(s), chunk %d/%d) reproduced the "
            "observed labels exactly -- this chunk contributed no spread; check the merged "
            "null before concluding the null itself has none",
            chunk_size,
            config.chunk,
            config.n_chunks,
        )
    return pd.DataFrame(rows, columns=["perm", "axis", "delta", "base_auc"])


class ControlSettingMismatch(ValueError):
    """A null was asked about a delta computed under a different control setting."""


def selected_maximum(deltas: pd.DataFrame) -> np.ndarray:
    """Per permutation, the LARGEST delta achieved by any offered axis.

    This is the multiplicity correction and the per-axis null is not: the maximum over
    several candidate axes is stochastically larger than any single one of them, so its
    distribution depends on how many axes were offered. That is why the candidate set must
    be the one actually searched, and why :meth:`NullDistribution.merge` refuses to merge
    chunks that offered different axes.

    Computed over the FULL frame, exactly one entry per permutation present, sorted by
    permutation index: ``groupby("perm").max()`` already skips a NaN axis WITHIN a
    permutation's group (an axis wholly inside the baseline was not a candidate on that
    permutation, so it must not drag the maximum down to it), and it yields NaN -- rather
    than dropping the permutation from the result -- for the one permutation whose every
    offered axis was NaN. Dropping such a permutation instead of recording NaN for it would
    leave the result shorter than the caller's permutation index, which is what
    :meth:`NullDistribution.from_deltas` relies on to stay aligned with ``perms``. Raises
    only when EVERY permutation is like that: then there is no finite delta anywhere for
    this candidate set to report.
    """
    sm = deltas.groupby("perm")["delta"].max().sort_index()
    if sm.isna().all():
        raise ValueError("no finite deltas: every offered axis was inside the baseline")
    return sm.to_numpy(dtype=float)


def _setting_mismatch(generated: ControlSetting, asked: ControlSetting) -> str:
    """Text shared by ``NullDistribution._require`` and ``merge``'s setting check.

    ``describe()`` alone can under-report a difference -- two candidate sets with equal
    per-axis column COUNTS render identically as ``{axis: n_cols}`` -- so this also names
    every differing field (via :func:`dataclasses.fields`) and prints both of its full
    values.
    """
    names = [
        f.name
        for f in fields(ControlSetting)
        if getattr(generated, f.name) != getattr(asked, f.name)
    ]
    detail = "; ".join(
        f"{name} (generated={getattr(generated, name)!r}, asked={getattr(asked, name)!r})"
        for name in names
    )
    return (
        f"generated under {generated.describe()}, but asked about {asked.describe()}. "
        f"Differing field(s): {detail}"
    )


@dataclass(frozen=True, eq=False)
class NullDistribution:
    """Per-axis and selected-maximum nulls, tagged with the setting they were built under.

    ``eq=False``: ``per_axis`` and ``selected_max`` are numpy arrays, whose ``==`` returns
    an array rather than a bool (dataclass-generated equality would therefore raise), and
    arrays are unhashable regardless. Two instances are never compared with ``==``;
    ``merge`` and ``_require`` compare the specific fields that matter (``setting``,
    ``axes``, ``perm_seed``, ``planned_n_perm``, ``observed``) explicitly.

    axes
        Axis names in OFFERED order, derived from ``setting.candidates``.
    per_axis
        ``{axis: np.ndarray of draws}``, each array ordered by permutation index.
    perms
        Sorted permutation indices actually present.
    perm_seed
        The ``NullConfig.seed`` the draws were generated with.
    planned_n_perm
        The ``NullConfig.n_perm`` of the WHOLE run -- chunks of one run share this even
        though each chunk's own ``n_perm`` (draws present) differs.
    observed
        The real-label deltas from the SAME scorer as the null, or ``None`` if this null
        was built without them.
    """

    setting: ControlSetting
    axes: tuple
    per_axis: dict
    selected_max: np.ndarray
    perms: tuple
    perm_seed: int
    planned_n_perm: int
    observed: "dict | None" = None

    @property
    def n_perm(self) -> int:
        return len(self.perms)

    @classmethod
    def from_deltas(
        cls,
        deltas: pd.DataFrame,
        setting: ControlSetting,
        *,
        null_config: NullConfig,
        observed: "dict | None" = None,
    ) -> "NullDistribution":
        """Build from the long-format output of :func:`permutation_deltas`.

        ``setting.candidates`` is the source of truth for which axes this null covers, and
        ``deltas`` must describe exactly the permutations ``null_config.span()`` says it
        should -- checked here, rather than assumed, because a null silently missing a
        permutation or an axis would still produce a plausible-looking p-value. ``axis`` and
        ``perm`` are coerced to ``str``/``int`` before anything else reads them: an integer
        axis name would otherwise compare unequal to the ``str`` names in
        ``setting.candidates`` (an empty, silently-dropped per-axis null), and a ``perm``
        column read back as text sorts ``"10"`` before ``"2"``.

        Every ``(perm, axis)`` pair must appear exactly once. The explicit duplicate check
        just below catches a chunk concatenated into ``deltas`` twice -- including the
        duplicate-plus-missing-cell case where a chunk both doubles one pair and drops
        another, leaving the total row count unchanged and so invisible to a grid-size
        check alone. The grid-completeness check further down
        (``len(deltas) == axes x permutations``) catches the remaining case: cells missing
        without an offsetting duplicate, where the row count falls short.

        ``observed`` is optional; passing it is what lets :meth:`merge` later detect chunks
        computed on different data, a different CV seed or a different estimator. If every
        chunk's ``observed`` is left ``None``, that check is silently off.
        """
        d = deltas.assign(
            axis=deltas["axis"].astype(str), perm=deltas["perm"].astype(int)
        )

        if d.duplicated(["perm", "axis"]).any():
            raise ValueError(
                "duplicate (perm, axis) pair(s): each must appear once; a chunk "
                "concatenated twice double-counts its draws and narrows every p-value"
            )

        axes = tuple(a for a, _ in setting.candidates)
        seen_axes = set(d["axis"].unique())
        if seen_axes != set(axes):
            raise ValueError(
                "the setting must describe the null it is attached to: setting.candidates "
                f"offers {list(axes)!r}, deltas contain {sorted(seen_axes)!r}"
            )
        if observed is not None:
            # Coerce keys to str first, exactly like `deltas["axis"]` above: an integer-keyed
            # `observed` (e.g. {0: ..., 1: ...}) would otherwise compare unequal to the str
            # axis names in `setting.candidates` and fail with a confusing "expected ['0',
            # '1'], got [0, 1]".
            observed = {str(k): v for k, v in observed.items()}
            observed_keys = set(observed)
            if observed_keys != set(axes):
                raise ValueError(
                    "observed must cover exactly the offered axes: expected "
                    f"{list(axes)!r}, got {sorted(observed_keys)!r}"
                )
            observed = {a: float(observed[a]) for a in axes}

        lo, hi = null_config.span()
        expected_perms = set(range(lo, hi))
        seen_perms = set(int(p) for p in d["perm"].unique())
        if seen_perms != expected_perms:
            raise ValueError(
                f"deltas must cover exactly the permutations null_config describes -- "
                f"expected range({lo}, {hi}), got {len(seen_perms)} distinct perm value(s)"
            )
        if len(d) != len(axes) * (hi - lo):
            raise ValueError(
                "deltas must be a complete grid of axes x permutations: expected "
                f"{len(axes)} axes x {hi - lo} permutations = {len(axes) * (hi - lo)} "
                f"rows, got {len(d)}"
            )

        per_axis = {
            a: d[d.axis == a].sort_values("perm")["delta"].to_numpy(dtype=float)
            for a in axes
        }
        return cls(
            setting=setting,
            axes=axes,
            per_axis=per_axis,
            selected_max=selected_maximum(d),
            perms=tuple(sorted(seen_perms)),
            perm_seed=null_config.seed,
            planned_n_perm=null_config.n_perm,
            observed=observed,
        )

    @classmethod
    def merge(cls, parts, *, allow_partial: bool = False) -> "NullDistribution":
        """Combine chunks of ONE run. Every disagreement is an error, not a reconciliation.

        Checked in this order, each an error rather than a reconciliation:

        1. the offered axes -- a selected-maximum over a candidate set nobody searched is
           a different statistic;
        2. the control setting -- ditto;
        3. the base seed and the planned permutation count -- chunks from different base
           seeds can draw the SAME underlying permutations from DISJOINT chunk indices
           (seed 10's chunk ``[20, 40)`` and seed 20's chunk ``[0, 20)`` share the 10 draws
           seeded 30..39), which the plain index-overlap check below cannot see;
        4. the observed deltas riding with the null, WHEN both sides carry one -- chunks
           whose ``observed`` disagree ran on different data, a different CV seed or a
           different estimator. Compared with the same tie tolerance ``p_value`` uses
           (genuinely different data moves a delta by at least the AUC grid spacing, far
           above it), not exact equality: nodes with different CPU instruction sets can
           otherwise differ in the last few bits and refuse a chunk that is actually fine.
           Passing ``observed`` at all is what turns this check on -- if every chunk's
           ``observed`` is ``None``, there is nothing here to disagree, and this check
           cannot catch chunks that were in fact computed on different data;
        5. overlapping permutation indices -- the same permutation counted twice narrows
           every p-value for free;
        6. a gap in the planned permutations -- ``planned_n_perm`` agrees between any two
           chunks of one run, so a chunk that never arrived (a pre-empted or failed array
           task) used to pass every check above and silently shrink the null to the draws
           present, with ``n_perm`` reporting the smaller count as if it were the plan.
           Refused unless ``allow_partial=True``, which keeps that behaviour for a caller
           who has decided a smaller null is acceptable and says so.
        """
        parts = list(parts)
        if not parts:
            raise ValueError("merge() needs at least one chunk")
        head = parts[0]

        for other in parts[1:]:
            if other.axes != head.axes:
                raise ValueError(
                    "cannot merge permutation chunks that offered different axes: "
                    f"{list(head.axes)!r} vs {list(other.axes)!r}; the selected-maximum "
                    "null is defined over the candidate set actually searched"
                )
        for other in parts[1:]:
            if other.setting != head.setting:
                raise ControlSettingMismatch(
                    "cannot merge permutation chunks: one chunk was "
                    + _setting_mismatch(head.setting, other.setting)
                )
        for other in parts[1:]:
            if (
                other.perm_seed != head.perm_seed
                or other.planned_n_perm != head.planned_n_perm
            ):
                raise ValueError(
                    "cannot merge permutation chunks generated from different base seeds "
                    f"or planned permutation counts: seed={head.perm_seed}, "
                    f"n_perm={head.planned_n_perm} vs seed={other.perm_seed}, "
                    f"n_perm={other.planned_n_perm}"
                )
        observed_list = [p.observed for p in parts]
        if any(o is not None for o in observed_list):
            if not all(o is not None for o in observed_list):
                raise ValueError(
                    "cannot merge permutation chunks where some carry observed deltas and "
                    "some do not"
                )
            head_obs = observed_list[0]
            for other_obs in observed_list[1:]:
                same = set(other_obs) == set(head_obs) and all(
                    np.allclose(
                        [head_obs[a]],
                        [other_obs[a]],
                        rtol=0,
                        atol=_TIE_TOL,
                        equal_nan=True,
                    )
                    for a in head_obs
                )
                if not same:
                    raise ValueError(
                        "cannot merge permutation chunks whose observed deltas disagree "
                        f"(the chunks ran on different data): {head_obs!r} vs {other_obs!r}"
                    )
        seen: set = set()
        for part in parts:
            overlap = seen & set(part.perms)
            if overlap:
                raise ValueError(
                    f"permutation indices overlap between chunks: {sorted(overlap)[:5]}...; "
                    "the same permutation counted twice narrows every p-value"
                )
            seen |= set(part.perms)
        planned = set(range(head.planned_n_perm))
        if seen != planned and not allow_partial:
            missing = sorted(planned - seen)
            shown = ", ".join(str(i) for i in missing[:10]) + (
                ", ..." if len(missing) > 10 else ""
            )
            raise ValueError(
                f"merged null covers {len(seen)} of the {head.planned_n_perm} planned "
                f"permutations; missing {len(missing)} index(es): [{shown}]. A chunk is "
                "absent (a pre-empted or failed array task?) -- rerun it, or pass "
                "allow_partial=True to accept the smaller null deliberately."
            )

        order = np.argsort(np.concatenate([np.asarray(p.perms) for p in parts]))
        return cls(
            setting=head.setting,
            axes=head.axes,
            per_axis={
                a: np.concatenate([p.per_axis[a] for p in parts])[order]
                for a in head.axes
            },
            selected_max=np.concatenate([p.selected_max for p in parts])[order],
            perms=tuple(
                int(v)
                for v in np.concatenate([np.asarray(p.perms) for p in parts])[order]
            ),
            perm_seed=head.perm_seed,
            planned_n_perm=head.planned_n_perm,
            observed=head.observed,
        )

    def _require(self, setting: ControlSetting) -> None:
        if setting != self.setting:
            raise ControlSettingMismatch(
                "this null was " + _setting_mismatch(self.setting, setting) + ". These "
                "are different statistics: a null built under one control setting does "
                "not bound a delta measured under another. Regenerate the null under the "
                "setting you are reporting."
            )

    def p_per_axis(
        self, axis: str, observed: float, *, setting: ControlSetting
    ) -> float:
        """Could THIS axis's delta arise by chance? Correct only if the axis was
        pre-specified; if it is reported because it came top, use ``p_selected_max``."""
        self._require(setting)
        if axis not in self.per_axis:
            raise KeyError(
                f"axis {axis!r} is not in this null; offered: {list(self.axes)!r}"
            )
        return p_value(self.per_axis[axis], observed)

    def p_selected_max(self, observed: float, *, setting: ControlSetting) -> float:
        """The multiplicity-corrected p-value: we searched ``len(self.axes)`` axes and
        reported the best -- could that have arisen by chance?"""
        self._require(setting)
        return p_value(self.selected_max, observed)

    def summary(self) -> pd.DataFrame:
        """One row per axis, plus the selected-maximum summary every row is judged against.

        When ``observed`` is set, adds ``observed``, ``p_per_axis`` and ``p_selected_max``:
        the Westfall-Young single-step max-T adjusted p-value for that axis's observed
        delta. Both go to NaN -- ``p_value`` is never called -- when the observed delta
        itself is NaN or this axis's per-axis null has no finite draws: that axis has
        nothing to report, not a real p-value that happens to look favourable.
        """
        sm = self.selected_max[np.isfinite(self.selected_max)]
        rows = []
        for axis in self.axes:
            draws = self.per_axis[axis]
            finite = draws[np.isfinite(draws)]
            row = {
                "axis": axis,
                "n_perm": self.n_perm,
                "n_candidates": len(self.axes),
                "null_mean": float(finite.mean()) if finite.size else float("nan"),
                "null_sd": float(finite.std()) if finite.size else float("nan"),
                "null_p95": float(np.percentile(finite, 95))
                if finite.size
                else float("nan"),
                "selmax_median": float(np.median(sm)) if sm.size else float("nan"),
                "selmax_p95": float(np.percentile(sm, 95)) if sm.size else float("nan"),
            }
            if self.observed is not None:
                obs = self.observed[axis]
                degenerate = not np.isfinite(obs) or finite.size == 0
                row["observed"] = obs
                row["p_per_axis"] = (
                    float("nan")
                    if degenerate
                    else self.p_per_axis(axis, obs, setting=self.setting)
                )
                row["p_selected_max"] = (
                    float("nan")
                    if degenerate
                    else self.p_selected_max(obs, setting=self.setting)
                )
            rows.append(row)
        return pd.DataFrame(rows)
