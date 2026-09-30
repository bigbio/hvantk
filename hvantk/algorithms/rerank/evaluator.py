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
"""The ablation's fold count. ``Config.folds`` governs ``ReRanker.score`` only."""


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
    """

    seeds: tuple
    deltas: tuple
    lo: float
    hi: float
    mean: float
    sd: float


def seed_sweep(matrix, base_cols, axis_cols, y, *, seeds, selector=None, groups=None,
                base_oof=None):
    """Delta-AUC of baseline+axis over baseline, recomputed under each CV seed.

    ``base_oof``, given as ``{seed: baseline OOF array}``, is used instead of refitting the
    baseline for a seed already present in it (a seed missing from the mapping is still
    computed). This is what lets ``Evaluator.evaluate`` fit the baseline once per sweep seed
    and share it across every axis, rather than refitting it once per (axis, seed) pair. The
    standalone call with no ``base_oof`` computes every seed's baseline itself and behaves
    exactly as the brief specifies.
    """
    from sklearn.metrics import roc_auc_score

    y = np.asarray(y)
    cc = list(dict.fromkeys(list(base_cols) + list(axis_cols)))
    deltas = []
    for seed in seeds:
        if base_oof is not None and seed in base_oof:
            p0 = base_oof[seed]
        else:
            p0 = _raw_oof(matrix, list(base_cols), y, selector, groups=groups, seed=seed)
        p1 = _raw_oof(matrix, cc, y, selector, groups=groups, seed=seed)
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
        self, matrix, feat_cols, y, scores, axis_groups, baseline_axis, selector=None,
        groups=None, seed=DEFAULT_SEED, n_seeds=1,
    ):
        y = np.asarray(y)
        scores = np.asarray(scores)
        auc = roc_auc_score(y, scores)
        pr = average_precision_score(y, scores)
        cs = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
        brier = brier_score_loss(y, cs)
        cal = calibration_curve(y, cs, n_bins=8, strategy="quantile")

        base_cols = axis_groups[baseline_axis]
        p_base = _raw_oof(matrix, base_cols, y, selector, groups=groups, seed=seed)
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
            base_oof_by_seed.update({
                s: _raw_oof(matrix, base_cols, y, selector, groups=groups, seed=s)
                for s in (seed + k for k in range(1, n_seeds))
            })
        rows = [base_row]
        spreads: dict = {}
        for fam, cols in axis_groups.items():
            if fam == baseline_axis:
                continue
            cc = list(dict.fromkeys(base_cols + cols))
            p = _raw_oof(matrix, cc, y, selector, groups=groups, seed=seed)
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
                    matrix, base_cols, cols, y,
                    seeds=tuple(seed + k for k in range(n_seeds)),
                    selector=selector, groups=groups, base_oof=base_oof_by_seed,
                )
                spreads[fam] = spread
                env_lo, env_hi = _envelope((lo, hi), spread)
                row["d_lo_env"] = round(env_lo, 3)
                row["d_hi_env"] = round(env_hi, 3)
                row["n_seeds"] = n_seeds
            rows.append(row)
        return EvalResult(
            round(auc, 3), round(pr, 3), round(brier, 4), pd.DataFrame(rows), cal,
            spreads or None,
        )
