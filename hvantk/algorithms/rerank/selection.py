"""Within-axis feature selection, run on a TRAINING SLICE ONLY.

Three steps, in the order benchmarked by Perez-Riverol et al. (PLOS ONE 2017,
doi:10.1371/journal.pone.0189875): univariate filter -> redundancy filter -> wrapper RFE.
Their Table 1 is the reason for that order: prefiltering cut RFE runtime ~3x (35 min -> 11)
at equal accuracy, and inner CV3/CV7/CV10 all gave the same RMSE, so 3 inner folds suffice.

Only the first TWO run by default here; the wrapper is opt-in (see ``SelectionPolicy``).
That is a departure from the reference workflow, taken on measurements from four cohorts
rather than on principle, and it is reversible per run.

EVERY function here takes a training slice. Nothing in this module may see held-out rows.
Selection that has seen the test labels inflates the reported metric (Ambroise & McLachlan
2002), and here the reported delta-AUC IS the scientific claim. `test_selection_nesting.py`
is the guardrail that proves the caller respects this.

Departure from the reference workflow: it used a fixed correlation cutoff, appropriate to a
continuous outcome. With a binary label and a few hundred positives, a fixed |AUC-0.5|
threshold is roughly a 1-SE gate -- for 219 positives against 19,229 negatives the null SE
is about 0.020 -- so noise passes. The univariate step uses within-axis BH-FDR instead,
which adapts to both cohort size and axis width with no per-cohort tuning.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class UnivariateStat:
    """Per-column univariate screening result.

    For untestable columns (all-NaN, or one arm empty) ``p`` defaults to ``1.0``
    while ``z`` is ``nan`` -- they are not interchangeable stand-ins for "no result".
    Consumers must gate on ``.passed``, never on ``np.isfinite(p)`` (always true,
    even for untestable columns) or on ``p`` alone.
    """

    auc: float
    z: float
    p: float
    passed: bool
    n_pos: int
    n_neg: int


def univariate_auc(x, y):
    """Rank-based AUC of a single feature, using only rows where it is defined.

    Returns ``(auc, n_pos_eff, n_neg_eff)``. Directionless use is intended: an
    anti-predictive feature (AUC near 0) is informative, and the caller scores on
    ``|auc - 0.5|``. Rank-based so that rankscore, phred and raw scales compare fairly.
    """
    from scipy.stats import rankdata

    x = np.asarray(x, dtype=float)
    y = np.asarray(y)
    ok = ~np.isnan(x)
    xs, ys = x[ok], y[ok]
    n_pos = int((ys == 1).sum())
    n_neg = int((ys == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), n_pos, n_neg
    ranks = rankdata(xs)
    auc = (ranks[ys == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return float(auc), n_pos, n_neg


def _null_se(n_pos: int, n_neg: int) -> float:
    """SE of AUC under the null. Bamber's bound; adapts the gate to cohort size."""
    return float(np.sqrt((1.0 / 12.0) * (1.0 / n_pos + 1.0 / n_neg)))


def univariate_filter(X, y, columns, q: float = 0.10):
    """Score each column by |AUC-0.5| and keep those surviving BH-FDR at ``q``.

    FDR is computed WITHIN the axis, so a 50-column dbNSFP axis is automatically held to a
    stricter bar than a 1-column expression axis. Columns that are entirely missing, or
    whose effective arm is empty, are reported as not passing rather than raising.
    """
    from scipy.stats import norm

    stats = {}
    for col in columns:
        auc, n_pos, n_neg = univariate_auc(X[col].to_numpy(), y)
        if np.isnan(auc) or n_pos == 0 or n_neg == 0:
            stats[col] = UnivariateStat(auc, float("nan"), 1.0, False, n_pos, n_neg)
            continue
        z = (auc - 0.5) / _null_se(n_pos, n_neg)
        p = 2.0 * norm.sf(abs(z))
        stats[col] = UnivariateStat(auc, float(z), float(p), False, n_pos, n_neg)

    testable = [c for c in columns if np.isfinite(stats[c].z)]
    if not testable:
        return stats
    order = sorted(testable, key=lambda c: stats[c].p)
    m = len(order)
    passing = set()
    for rank, col in enumerate(order, start=1):
        if stats[col].p <= q * rank / m:
            passing.update(order[:rank])   # BH: everything up to the largest passing rank
    for col in passing:
        s = stats[col]
        stats[col] = UnivariateStat(s.auc, s.z, s.p, True, s.n_pos, s.n_neg)
    return stats


def redundancy_filter(X, columns, scores, max_rho: float = 0.75):
    """Drop columns that duplicate a stronger sibling WITHIN the same axis.

    Greedy in descending univariate strength: walk the columns best-first and drop any
    whose |Spearman rho| against an already-kept column reaches ``max_rho``. The 0.75
    default is the correlation-matrix cutoff from the reference workflow.

    Spearman, not Pearson: these predictors are monotonically but not linearly related
    (rankscore vs phred vs raw), so a linear coefficient would understate the redundancy.

    Deliberately WITHIN-axis only. An axis that duplicates the baseline axis is not pruned
    here -- the per-axis delta-AUC exists to reveal exactly that, and pruning it first
    would hide the finding.
    """
    import pandas as pd

    ordered = sorted(columns, key=lambda c: -abs(scores.get(c, 0.0)))
    kept: list[str] = []
    dropped: dict[str, str] = {}
    for col in ordered:
        redundant_with = None
        for k in kept:
            pair = X[[col, k]].dropna()
            if len(pair) < 3:
                continue
            rho = pair[col].corr(pair[k], method="spearman")
            if pd.notna(rho) and abs(rho) >= max_rho:
                redundant_with = k
                break
        if redundant_with is None:
            kept.append(col)
        else:
            dropped[col] = f"redundant_with:{redundant_with}"
    return kept, dropped


@dataclass(frozen=True)
class SelectionPolicy:
    """How to filter one axis. An all-defaults instance is the recommended policy.

    ``wrapper`` defaults to "none" -- the two filters only -- on measured evidence, not on
    principle. Across four real cohorts (CHD, CHD-NDD, epilepsy-DEE, SCHEMA-SCZ; 21 axes)
    RFECV eliminated a column on 4 of the 11 axes wide enough for it to run, and 3 of
    those 4 were in the cohort with the FEWEST positives (54). Its aggression tracks label
    scarcity inversely, which is the signature of a wrapper fitting inner-CV noise: at 54
    positives an inner CV3 fold holds ~18, and "the feature count that maximised inner
    AUC" is barely distinguishable from chance.

    Worse, it prunes the axis whose composition defines the headline metric. On CHD-NDD it
    cut the constraint axis to a single column, and constraint is the ABLATION BASELINE --
    a thinner baseline silently inflates every other axis's delta-AUC. On CHD it removed
    ``n_case_var``, one of the two columns that condition the gene universe.

    Set ``wrapper="rfecv"`` deliberately, on an axis wide enough to need it (the ~45-column
    dbNSFP predictor axis is the motivating case) and with enough positives to trust the
    inner CV. It is kept, not deleted, because it demonstrably has behaviour -- it has
    simply not yet been shown to have BENEFICIAL behaviour, which needs an out-of-fold
    outcome comparison rather than an elimination count.
    """

    univariate: str = "auc"          # "auc" | "none"
    q: float = 0.10                  # BH-FDR level, within axis
    redundancy: str = "spearman"     # "spearman" | "none"
    redundancy_max: float = 0.75
    wrapper: str = "none"            # "none" | "rfecv" -- see the class docstring
    wrapper_estimator: str = "random_forest"
    inner_folds: int = 3             # Table 1: CV3 == CV7 == CV10; more is wasted compute
    seed: int = 42


@dataclass(frozen=True)
class SelectionReport:
    kept: tuple[str, ...]
    dropped: dict[str, str]
    stats: dict
    wrapper_ran: bool


def _rfecv(X, y, columns, policy):
    """Recursive feature elimination with internal CV, choosing the count itself.

    RandomForest rather than the downstream HistGradientBoostingClassifier: HistGBM
    exposes neither ``coef_`` nor ``feature_importances_``, so sklearn's RFE cannot wrap
    it. The reference workflow used RandomForest and SVM, so this follows it -- but the
    selector and the final model are then DIFFERENT estimators, and features RF ranks as
    important are not guaranteed to be the ones HistGBM would favour. Recorded as a known
    caveat rather than presented as free.

    Requires scikit-learn >= 1.4, where tree estimators accept NaN natively; these matrices
    are 20-50% missing by design and imputing here would leak column statistics.
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.feature_selection import RFECV
    from sklearn.model_selection import StratifiedKFold

    est = RandomForestClassifier(
        n_estimators=200, min_samples_leaf=20, class_weight="balanced",
        random_state=policy.seed, n_jobs=-1,
    )
    sel = RFECV(
        estimator=est,
        step=1,
        min_features_to_select=1,
        cv=StratifiedKFold(policy.inner_folds, shuffle=True, random_state=policy.seed),
        scoring="roc_auc",
    )
    sel.fit(X[list(columns)].to_numpy(), y)
    return [c for c, keep in zip(columns, sel.support_) if keep]


def select_axis(X, y, columns, policy) -> SelectionReport:
    """Filter one axis down to its informative, non-redundant subset.

    ``X`` must already be restricted to the TRAINING slice. Order is univariate ->
    redundancy -> wrapper, which is the benchmarked order: the filters make the wrapper
    roughly 3x cheaper at equal accuracy.
    """
    columns = list(columns)
    dropped: dict[str, str] = {}
    stats: dict = {}

    kept = columns
    if policy.univariate == "auc":
        stats = univariate_filter(X, y, kept, q=policy.q)
        survivors = [c for c in kept if stats[c].passed]
        for c in kept:
            if c not in survivors:
                dropped[c] = "univariate_fdr"
        kept = survivors

    if policy.redundancy == "spearman" and len(kept) > 1:
        strength = {c: abs(stats[c].auc - 0.5) if c in stats else 0.0 for c in kept}
        kept, red_dropped = redundancy_filter(X, kept, strength, policy.redundancy_max)
        dropped.update(red_dropped)

    wrapper_ran = False
    n_pos = int((np.asarray(y) == 1).sum())
    n_neg = int((np.asarray(y) == 0).sum())
    if policy.wrapper == "rfecv" and len(kept) > 1:
        if min(n_pos, n_neg) < policy.inner_folds:
            # Degrade to filter-only rather than raise: a small arm is a property of the
            # cohort, not an error, and the report records that the wrapper was skipped.
            wrapper_ran = False
        else:
            survivors = _rfecv(X, y, kept, policy)
            for c in kept:
                if c not in survivors:
                    dropped[c] = "rfecv"
            kept = survivors
            wrapper_ran = True

    return SelectionReport(tuple(kept), dropped, stats, wrapper_ran)


_POLICY_SCHEMA_PATH = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "schemas"
    / "selection_policy.schema.json"
)


def load_policy(path):
    """Read and validate a selection policy. Returns ``(SelectionPolicy, equivalence)``.

    Every field is optional: an empty document is a valid, working policy that yields the
    defaults above. That matters because it makes the feature-selection stage adoptable
    without anyone first having to understand the knobs.

    ``equivalence`` is provenance vocabulary rather than a filter knob, so it is returned
    alongside the policy instead of living on it -- the resolver consumes it, the filters
    never see it.
    """
    import json

    import jsonschema
    import yaml

    from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE

    doc = yaml.safe_load(Path(path).read_text()) or {}
    jsonschema.validate(doc, json.loads(_POLICY_SCHEMA_PATH.read_text()))
    equivalence = doc.pop("equivalence", None) or DEFAULT_EQUIVALENCE
    return SelectionPolicy(**doc), equivalence
