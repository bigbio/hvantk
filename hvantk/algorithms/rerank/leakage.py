"""Bar feature columns whose MISSINGNESS predicts the label.

Circularity, handled in ``provenance.py``, asks what a predictor was TRAINED on. This asks
which genes it was RUN on. They are independent, and only the first had a control.

The case that motivated this module is EVE. It is computed from deep multiple-sequence
alignments for a curated subset of proteins; that subset is enriched for well-studied
genes; well-studied genes are disease genes. So the bare flag "was EVE computed for this
gene" can score a real AUC against a ClinGen/GenCC-derived label -- on its own, without
any EVE score. Gradient-boosted trees learn a default branch direction for NaN, so "not
computed" reaches the model as a usable feature. EVE is unsupervised on alignments and so
passes every circularity check correctly.

Note what is NOT the defect. Sparsity is not: a column can be sparser than EVE and still
have an inert presence indicator (presence AUC near 0.5), when whatever decides its
coverage is unrelated to the label. A coverage threshold cannot separate the two, and one
strict enough to catch EVE would delete sparse, non-leaking axes for the wrong reason. What
matters is missingness that CORRELATES WITH THE LABEL -- which, exactly like circularity,
makes this a property of the feature/LABEL pair rather than of the feature. The same column
can leak strongly against one disease label and weakly against another.

Statistical selection cannot find this on its own; it rewards it. A univariate filter ranks
EVE highly BECAUSE its missingness tracks the label, which is the argument
``provenance.py`` already makes about circular predictors, transposed.

Two conditions are required to bar a column, because either alone fails at one end of the
sample-size range:

* an FDR-significant presence/label association -- effect size alone flags noise in a small
  cohort;
* a minimum effect -- significance alone flags a 1.5-point coverage difference once n is
  large enough, and near-fully-covered reference axes sit just above presence AUC 0.5;
  they must never be barred, or the control deletes the baseline it exists to protect.

Dropping is the bluntest available remedy and deliberately not the only one. A column whose
presence predicts the label can also be handled by restricting the universe to where it is
measured, or by admitting the presence indicator as an explicit feature and requiring the
score to beat it -- which separates "the score is informative" from "having the score is
informative". This module implements the blunt one because it is what a headline result can
be defended with.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hvantk.algorithms.rerank.selection import _null_se, univariate_auc

# Presence AUC at or beyond 0.5 +/- this is "material". An empirical default, not a derived
# bound: near-fully-covered reference axes (constraint, expression) land within a few
# hundredths of 0.5, so 0.05 sits just above what a non-leaking column reaches, while a
# presence flag that genuinely tracks the label (EVE) lands well outside it. Re-check it
# for a different spine or label.
DEFAULT_MIN_EFFECT = 0.05

# Presence AUC threshold in the same units callers think in.
DEFAULT_MIN_AUC = 0.5 + DEFAULT_MIN_EFFECT


@dataclass(frozen=True)
class LeakagePolicy:
    """How aggressively to bar columns whose missingness predicts the label.

    An all-defaults instance is the recommended policy. Both knobs must be cleared for a
    column to be barred -- see the module docstring for why either alone is wrong at one
    end of the sample-size range.
    """

    q: float = 0.10
    """BH-FDR level for the presence/label association, across the columns supplied."""

    min_auc: float = DEFAULT_MIN_AUC
    """Presence AUC beyond which the association is material. Two-sided about 0.5."""


@dataclass(frozen=True)
class LeakageStat:
    """Presence/label association for one column.

    ``auc`` is ``nan`` when the test is undefined -- a fully observed column, a fully
    missing one, or a single-class label. Undefined is not the same as null: consumers must
    gate on ``.leaking``, never on ``auc`` alone.
    """

    coverage: float
    auc: float
    z: float
    p: float
    leaking: bool
    n_pos: int
    n_neg: int


@dataclass(frozen=True)
class LeakageReport:
    """A partition of the supplied columns. ``clean`` and ``leaking`` are disjoint and
    together account for every column asked about."""

    clean: tuple[str, ...]
    leaking: dict[str, float]
    stats: dict[str, LeakageStat]


def presence_leakage(
    X, y, columns, q: float = 0.10, min_auc: float = DEFAULT_MIN_AUC
) -> dict[str, LeakageStat]:
    """Score every column by how well its presence indicator alone predicts ``y``.

    No score VALUES are used -- only whether each cell is observed.
    """
    from scipy.stats import norm

    if not 0.0 < q <= 1.0:
        raise ValueError(f"q must be in (0, 1]; got {q}")
    # NaN would pass silently and disable the effect floor entirely: max(0.0, nan) is 0.0,
    # so every FDR-significant column would be barred with no effect-size protection --
    # exactly the over-correction the floor exists to prevent, and invisible in the output.
    _m = float(min_auc)
    if not np.isfinite(_m) or not 0.0 <= _m <= 1.0:
        raise ValueError(f"min_auc must be finite and in [0, 1]; got {min_auc!r}")

    y = np.asarray(y)
    min_effect = max(0.0, _m - 0.5)
    stats: dict[str, LeakageStat] = {}

    for col in columns:
        present = X[col].notna().to_numpy()
        coverage = float(present.mean()) if len(present) else 0.0
        # A column with no missingness has no indicator, and one with nothing but
        # missingness has no contrast. Neither can leak; both are undefined, not null.
        if coverage in (0.0, 1.0):
            stats[col] = LeakageStat(
                coverage,
                float("nan"),
                float("nan"),
                1.0,
                False,
                int((y == 1).sum()),
                int((y == 0).sum()),
            )
            continue
        auc, n_pos, n_neg = univariate_auc(present.astype(float), y)
        if np.isnan(auc) or n_pos == 0 or n_neg == 0:
            stats[col] = LeakageStat(
                coverage, float("nan"), float("nan"), 1.0, False, n_pos, n_neg
            )
            continue
        z = (auc - 0.5) / _null_se(n_pos, n_neg)
        # Two-sided: presence predicting the NEGATIVE class is informative missingness too.
        # A column measured only on controls leaks exactly as hard as one measured only on
        # cases, and a one-sided test would wave it through.
        p = 2.0 * norm.sf(abs(z))
        stats[col] = LeakageStat(
            coverage, float(auc), float(z), float(p), False, n_pos, n_neg
        )

    # BH-FDR across the columns actually supplied, mirroring univariate_filter: a 115-column
    # axis gets more chances at a spurious hit than a 1-column axis, so the correction has
    # to see how many were tried.
    testable = [c for c in columns if np.isfinite(stats[c].z)]
    if not testable:
        return stats
    order = sorted(testable, key=lambda c: stats[c].p)
    m = len(order)
    significant: set[str] = set()
    for rank, col in enumerate(order, start=1):
        if stats[col].p <= q * rank / m:
            significant.update(order[:rank])

    for col in significant:
        s = stats[col]
        if abs(s.auc - 0.5) >= min_effect:
            stats[col] = LeakageStat(
                s.coverage, s.auc, s.z, s.p, True, s.n_pos, s.n_neg
            )
    return stats


def resolve_leakage(
    X, y, columns, q: float = 0.10, min_auc: float = DEFAULT_MIN_AUC
) -> LeakageReport:
    """Split columns into those safe to use and those whose missingness carries the label."""
    stats = presence_leakage(X, y, columns, q=q, min_auc=min_auc)
    clean = tuple(c for c in columns if not stats[c].leaking)
    leaking = {c: stats[c].auc for c in columns if stats[c].leaking}
    return LeakageReport(clean=clean, leaking=leaking, stats=stats)


def leakage_selector(q: float = 0.10, min_auc: float = DEFAULT_MIN_AUC):
    """A selector for ``_raw_oof(..., selector=...)``, so the control runs PER FOLD.

    Computing leakage once over all rows and then filtering would let held-out labels choose
    the feature set -- the error ``selection.py``'s module docstring exists to prevent, and
    the reported delta-AUC is the scientific claim. Returned columns keep the caller's
    order, since downstream code indexes matrices positionally.
    """

    def _select(X_tr, y_tr, cols):
        return resolve_leakage(X_tr, y_tr, list(cols), q=q, min_auc=min_auc).clean

    return _select
