"""Headline metrics and the per-axis ablation table for one rerank run.

Three numbers describe the model as a whole (ROC AUC, average precision, Brier on
min-max-rescaled scores, plus a quantile-binned calibration curve), and then one row per
feature axis describes what that axis ADDS: the out-of-fold AUC of baseline+axis, with a
paired bootstrap interval on the difference from the baseline alone.

``_raw_oof`` is uncalibrated on purpose. The ablation compares discrimination between two
feature sets, and isotonic calibration is a monotone transform -- it cannot change AUC, and
fitting it 5x per fold would multiply the cost of the table by the calibration folds for no
change in the answer. ``ReRanker.score`` still calibrates, because the SCORES it produces
are read as probabilities; these deltas are not.

The per-axis interval here resamples GENES only. It is one of at least three variance
components; the CV partition and the best-of-N axis selection are the other two, and they
live in ``nulls.py`` and in ``seed_sweep`` below.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import cross_val_predict

from hvantk.algorithms.rerank.reranker import _cv, _gbm, _grouped_splits
from hvantk.algorithms.rerank.seeds import DEFAULT_SEED


ABLATION_FOLDS = 5
"""The ablation's fold count when no caller threads one through.

The engine always passes ``Config.folds`` explicitly -- to ``Evaluator.evaluate``,
to ``seed_sweep`` and, via ``engine._run_nulls``, to ``nulls.oof_scorer`` -- so the
headline score, the ablation and the permutation null all use the SAME fold count on
every path `hvantk rerank`/`rerank()` take. This constant is only still read as a
default by a direct call to ``_raw_oof``/``seed_sweep``/``Evaluator.evaluate`` that
passes no ``folds=`` of its own (e.g. a script calling this module directly)."""


def _raw_oof(
    matrix, cols, y, selector=None, groups=None, seed=DEFAULT_SEED, folds=ABLATION_FOLDS
):
    """Uncalibrated out-of-fold probabilities for one feature subset.

    Mirrors ReRanker.score's nesting contract: when a selector is given it runs per fold on
    the training slice only, so every ablation delta-AUC is as leakage-free as the headline.

    ``groups`` switches to paralogue-blocked folds. It is threaded to the ABLATION path as
    well as the headline deliberately: a blocked headline beside unblocked deltas would
    report a corrected level and an uncorrected claim.
    """
    y = np.asarray(y)
    cv = _cv(folds, seed, groups)
    X = matrix[cols].values
    if selector is None:
        cv_arg = cv if groups is None else _grouped_splits(cv, X, y, groups)
        return cross_val_predict(
            _gbm(seed), X, y, cv=cv_arg, groups=groups, method="predict_proba"
        )[:, 1]

    oof = np.full(len(y), np.nan)
    splits = cv.split(X, y) if groups is None else _grouped_splits(cv, X, y, groups)
    for train_idx, test_idx in splits:
        X_tr = matrix.iloc[train_idx]
        sel = list(selector(X_tr, y[train_idx], list(cols)))
        if not sel:
            oof[test_idx] = y[train_idx].mean()
            continue
        m = _gbm(seed).fit(X_tr[sel].values, y[train_idx])
        oof[test_idx] = m.predict_proba(matrix.iloc[test_idx][sel].values)[:, 1]
    return oof


@dataclass
class EvalResult:
    auc: float
    pr_auc: float
    brier: float
    ablation: pd.DataFrame
    calibration: tuple
    seed_spread: dict = None
    """axis -> SeedSpread, or ``None`` when the run was not swept (``n_seeds == 1``,
    the default). Populated by ``Evaluator.evaluate``'s ``n_seeds`` argument."""


@dataclass(frozen=True)
class SeedSpread:
    """One axis's delta-AUC across several CV partitions.

    Reported BESIDE the bootstrap interval, never instead of it: they measure different
    things. The bootstrap asks how the delta would move if the gene SAMPLE moved; this asks
    how it moves when only the fold ASSIGNMENT does -- the CV partition is a second variance
    component the bootstrap cannot see at all. The gap between a single partition's delta and
    the across-seed mean can be larger than several axes' entire measured gain, which is why
    it is reported as a separate number rather than assumed away.

    "Only the fold assignment" holds exactly for fits of up to 10,000 rows. Above that,
    ``_gbm`` leaves ``HistGradientBoostingClassifier``'s ``early_stopping="auto"`` and
    scikit-learn (1.7.2) switches it on, so each fit also holds out a random, unblocked
    validation split drawn with the same seed -- the sweep then varies that split as well
    as the partition (see ``ReRanker.score`` for the limit and why it is kept).
    """

    seeds: tuple
    deltas: tuple
    lo: float
    hi: float
    mean: float
    sd: float

    def __post_init__(self) -> None:
        """The two relationships ``seed_sweep`` already guarantees: one delta per seed,
        and a minimum that is at most the maximum it was computed alongside."""
        if len(self.seeds) != len(self.deltas):
            raise ValueError(
                f"seeds has {len(self.seeds)} entrie(s) but deltas has "
                f"{len(self.deltas)}; seed_sweep records exactly one delta per seed"
            )
        # `not (lo <= hi)` rather than `lo > hi`: every comparison with NaN is False,
        # so the latter would let a NaN bound through.
        if not (self.lo <= self.hi):
            raise ValueError(f"lo ({self.lo!r}) must be <= hi ({self.hi!r})")


def seed_sweep(
    matrix,
    base_cols,
    axis_cols,
    y,
    *,
    seeds,
    selector=None,
    groups=None,
    base_oof=None,
    aug_oof=None,
    folds=ABLATION_FOLDS,
):
    """Delta-AUC of baseline+axis over baseline, recomputed under each CV seed.

    ``base_oof``/``aug_oof``, each given as ``{seed: OOF array}`` (baseline-only /
    baseline+axis respectively), are used instead of refitting for a seed already present in
    them (a seed missing from either mapping is still computed). This is what lets
    ``Evaluator.evaluate`` fit the baseline once per sweep seed and share it across every
    axis, and reuse the headline seed's already-fitted augmented model here instead of
    refitting it a second time, rather than refitting both once per (axis, seed) pair. The
    standalone call with neither mapping computes every seed's baseline and augmented fit
    itself.

    ``folds`` must be the SAME fold count any cached ``base_oof``/``aug_oof`` arrays were
    fitted with -- ``Evaluator.evaluate`` passes its own ``folds`` through here for exactly
    that reason. Left at ``ABLATION_FOLDS`` for a standalone call that threads none.
    """
    from sklearn.metrics import roc_auc_score

    y = np.asarray(y)
    cc = list(dict.fromkeys(list(base_cols) + list(axis_cols)))
    deltas = []
    for seed in seeds:
        if base_oof is not None and seed in base_oof:
            p0 = base_oof[seed]
        else:
            p0 = _raw_oof(
                matrix,
                list(base_cols),
                y,
                selector,
                groups=groups,
                seed=seed,
                folds=folds,
            )
        if aug_oof is not None and seed in aug_oof:
            p1 = aug_oof[seed]
        else:
            p1 = _raw_oof(
                matrix, cc, y, selector, groups=groups, seed=seed, folds=folds
            )
        deltas.append(float(roc_auc_score(y, p1) - roc_auc_score(y, p0)))
    arr = np.asarray(deltas, dtype=float)
    return SeedSpread(
        seeds=tuple(int(s) for s in seeds),
        deltas=tuple(deltas),
        lo=float(arr.min()),
        hi=float(arr.max()),
        mean=float(arr.mean()),
        sd=float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
    )


def _envelope(boot, spread):
    """The union of the gene-resampling interval and the across-seed range.

    A union, NOT a calibrated interval over both components: the two are not independent
    draws from one distribution and combining them properly needs a nested design nobody has
    run. It is reported as an ENVELOPE, is guaranteed never to be narrower than the bootstrap
    interval alone, and exists so that a reader cannot mistake a lucky partition for a
    measured effect.

    Only these two components go into it. The permutation null (``nulls.py``) already
    recomputes the stratified partition for every permutation, so it already contains
    partition noise as one of its components; this envelope is a separate report on the
    OBSERVED delta and must never be folded into the null on top of that, or partition
    variance would be counted twice.
    """
    lo, hi = boot
    if spread is None:
        return lo, hi
    return min(lo, spread.lo), max(hi, spread.hi)


def _boot_ci(y, p1, p0, n=1000, seed=DEFAULT_SEED):
    """Paired bootstrap over GENES for the delta-AUC of ``p1`` over ``p0``.

    Draws that end up single-class are skipped rather than counted: AUC is undefined there,
    and substituting 0.5 would drag the interval toward no-difference.
    """
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y))
    d = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), True)
        if 0 < y[b].sum() < len(b):
            d.append(roc_auc_score(y[b], p1[b]) - roc_auc_score(y[b], p0[b]))
    return [round(v, 3) for v in np.percentile(d, [2.5, 50, 97.5])]


class Evaluator:
    def evaluate(
        self,
        matrix,
        feat_cols,
        y,
        scores,
        axis_groups,
        baseline_axis,
        selector=None,
        groups=None,
        seed=DEFAULT_SEED,
        n_seeds=1,
        folds=ABLATION_FOLDS,
    ):
        """Headline metrics plus the per-axis ablation table (see the module docstring).

        ``folds`` governs every out-of-fold fit THIS method makes (the baseline, each
        axis, and -- via ``seed_sweep`` -- every swept seed's refit): one fold count for
        the whole ablation. The engine passes ``Config.folds`` here so the ablation
        matches the headline score and the permutation null (``nulls.oof_scorer`` is
        built from the same ``Config.folds`` in ``engine._run_nulls``); a direct caller
        that passes none gets ``ABLATION_FOLDS``, the historic fixed value.
        """
        y = np.asarray(y)
        scores = np.asarray(scores)
        auc = roc_auc_score(y, scores)
        pr = average_precision_score(y, scores)
        cs = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
        brier = brier_score_loss(y, cs)
        cal = calibration_curve(y, cs, n_bins=8, strategy="quantile")

        base_cols = axis_groups[baseline_axis]
        p_base = _raw_oof(
            matrix, base_cols, y, selector, groups=groups, seed=seed, folds=folds
        )
        base_auc = roc_auc_score(y, p_base)
        base_row = {
            "family": baseline_axis,
            "auc": round(base_auc, 3),
            "d_lo": 0.0,
            "d_md": 0.0,
            "d_hi": 0.0,
        }
        # The baseline OOF is fit once per sweep seed, reusing the headline fit for the first seed.
        # With A axes and S seeds that is S baseline fits where A x S would otherwise suffice.
        base_oof_by_seed = None
        if n_seeds > 1:
            base_row["d_lo_env"] = 0.0
            base_row["d_hi_env"] = 0.0
            base_row["n_seeds"] = n_seeds
            base_oof_by_seed = {seed: p_base}
            base_oof_by_seed.update(
                {
                    s: _raw_oof(
                        matrix,
                        base_cols,
                        y,
                        selector,
                        groups=groups,
                        seed=s,
                        folds=folds,
                    )
                    for s in (seed + k for k in range(1, n_seeds))
                }
            )
        rows = [base_row]
        spreads: dict = {}
        for fam, cols in axis_groups.items():
            if fam == baseline_axis:
                continue
            cc = list(dict.fromkeys(base_cols + cols))
            p = _raw_oof(matrix, cc, y, selector, groups=groups, seed=seed, folds=folds)
            lo, md, hi = _boot_ci(y, p, p_base, seed=seed)
            row = {
                "family": fam,
                "auc": round(roc_auc_score(y, p), 3),
                "d_lo": lo,
                "d_md": md,
                "d_hi": hi,
            }
            if n_seeds > 1:
                # seed, seed+1, ... rather than a spawned stream: the same reproducibility
                # argument as rng_for (see seeds.py).
                spread = seed_sweep(
                    matrix,
                    base_cols,
                    cols,
                    y,
                    seeds=tuple(seed + k for k in range(n_seeds)),
                    selector=selector,
                    groups=groups,
                    base_oof=base_oof_by_seed,
                    aug_oof={seed: p},
                    folds=folds,
                )
                spreads[fam] = spread
                env_lo, env_hi = _envelope((lo, hi), spread)
                row["d_lo_env"] = round(env_lo, 3)
                row["d_hi_env"] = round(env_hi, 3)
                row["n_seeds"] = n_seeds
            rows.append(row)
        return EvalResult(
            round(auc, 3),
            round(pr, 3),
            round(brier, 4),
            pd.DataFrame(rows),
            cal,
            spreads or None,
        )
