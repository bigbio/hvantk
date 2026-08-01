"""Presence-leakage control: a column whose MISSINGNESS predicts the label.

Motivating measurement (CHD cohort, 1,362 genes): the bare flag "was EVE computed for this
gene" scores AUC 0.716 against a ClinGen/GenCC-derived label -- higher than the entire
nine-feature gnomAD constraint axis. EVE is computed from deep alignments for a curated
subset of proteins, that subset is enriched for well-studied genes, and well-studied genes
are disease genes. A gradient-boosted tree learns a default direction for NaN, so
"not computed" is available to it as a feature.

Critically, EVE is unsupervised on multiple-sequence alignments and therefore PASSES the
existing circularity check in provenance.py: that machinery asks what a predictor was
TRAINED on. This asks which genes it was RUN on. The two are independent, and the second
had no control.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hvantk.algorithms.rerank.leakage import (
    LeakageStat,
    presence_leakage,
    resolve_leakage,
)


def _frame(n=600, seed=0):
    rng = np.random.default_rng(seed)
    y = np.zeros(n, dtype=int)
    y[: n // 6] = 1  # ~17% prevalence
    rng.shuffle(y)
    return rng, y


def test_column_present_only_for_positives_is_flagged_as_leaking():
    """The EVE case: a column measured mostly on labelled genes."""
    rng, y = _frame()
    x = rng.normal(size=len(y))
    # present for 90% of positives but only 20% of negatives
    keep = np.where(y == 1, rng.random(len(y)) < 0.9, rng.random(len(y)) < 0.2)
    x[~keep] = np.nan
    X = pd.DataFrame({"eve_mean": x})

    report = resolve_leakage(X, y, ["eve_mean"])

    assert report.leaking, "a column present mainly on positives must be flagged"
    assert "eve_mean" in report.leaking
    assert "eve_mean" not in report.clean
    assert report.stats["eve_mean"].auc > 0.6


def test_column_missing_at_random_is_clean():
    """Sparsity alone is not the defect. GTEx eQTL is sparser than EVE and inert."""
    rng, y = _frame()
    x = rng.normal(size=len(y))
    x[rng.random(len(y)) < 0.8] = np.nan  # 20% coverage, unrelated to y
    X = pd.DataFrame({"eqtl_beta": x})

    report = resolve_leakage(X, y, ["eqtl_beta"])

    assert "eqtl_beta" in report.clean
    assert not report.leaking


def test_fully_observed_column_cannot_leak():
    """With no missingness there is no indicator, so the test is undefined, not failed."""
    rng, y = _frame()
    X = pd.DataFrame({"pLI": rng.normal(size=len(y))})

    report = resolve_leakage(X, y, ["pLI"])

    assert "pLI" in report.clean
    assert report.stats["pLI"].coverage == 1.0
    assert np.isnan(report.stats["pLI"].auc)


def test_fully_missing_column_cannot_leak_either():
    rng, y = _frame()
    X = pd.DataFrame({"empty": np.full(len(y), np.nan)})

    report = resolve_leakage(X, y, ["empty"])

    assert "empty" in report.clean
    assert report.stats["empty"].coverage == 0.0


def test_leakage_is_detected_in_both_directions():
    """Presence predicting the NEGATIVE class is equally informative missingness.

    A one-sided test would wave through a column measured only on controls, which leaks
    just as hard as one measured only on cases.
    """
    rng, y = _frame()
    x = rng.normal(size=len(y))
    keep = np.where(y == 0, rng.random(len(y)) < 0.9, rng.random(len(y)) < 0.2)
    x[~keep] = np.nan
    X = pd.DataFrame({"inverted": x})

    report = resolve_leakage(X, y, ["inverted"])

    assert "inverted" in report.leaking
    assert report.stats["inverted"].auc < 0.4


def test_small_effect_is_not_flagged_even_when_significant():
    """Significance is not sufficient. At large n a trivial imbalance reaches p<0.05.

    The reference axes in the motivating cohort (constraint, expression) sit at presence
    AUC 0.503-0.530 and must never be barred, or the control would delete the baseline it
    exists to protect.
    """
    rng = np.random.default_rng(7)
    n = 200_000
    y = (rng.random(n) < 0.2).astype(int)
    x = rng.normal(size=n)
    # a 1.5-point coverage difference: overwhelming p, negligible effect
    keep = np.where(y == 1, rng.random(n) < 0.815, rng.random(n) < 0.80)
    x[~keep] = np.nan
    X = pd.DataFrame({"barely": x})

    report = resolve_leakage(X, y, ["barely"])

    assert abs(report.stats["barely"].auc - 0.5) < 0.05
    assert "barely" in report.clean, "effect-size floor must protect near-null columns"


def test_large_effect_is_not_flagged_when_it_could_be_chance():
    """Effect size is not sufficient either. At tiny n a big AUC is unremarkable.

    Requiring BOTH an FDR-significant test and an effect floor is the point: either alone
    is wrong at one end of the sample-size range.
    """
    rng = np.random.default_rng(3)
    n, y = 24, np.array([1] * 6 + [0] * 18)
    x = rng.normal(size=n)
    x[[0, 1, 7, 8, 9, 10]] = np.nan
    X = pd.DataFrame({"tiny": x})

    report = resolve_leakage(X, y, ["tiny"], min_auc=0.0)  # effect floor disabled

    # With 6 positives nothing should clear BH-FDR on its own.
    assert "tiny" in report.clean


def test_fdr_is_computed_across_the_supplied_columns():
    """Testing 50 columns and testing 1 must not use the same bar.

    Same argument univariate_filter already makes for within-axis FDR: a wide axis gets
    more chances to throw a spurious hit, so the correction has to see how many were tried.
    """
    rng, y = _frame(n=800)
    cols, data = [], {}
    for i in range(50):
        x = rng.normal(size=len(y))
        x[rng.random(len(y)) < 0.5] = np.nan  # all missing-at-random
        data[f"noise{i}"] = x
        cols.append(f"noise{i}")
    X = pd.DataFrame(data)

    report = resolve_leakage(X, y, cols)

    assert not report.leaking, f"MAR columns must not be flagged: {report.leaking}"


def test_presence_leakage_returns_a_stat_for_every_column():
    rng, y = _frame()
    X = pd.DataFrame({"a": rng.normal(size=len(y)), "b": rng.normal(size=len(y))})

    stats = presence_leakage(X, y, ["a", "b"])

    assert set(stats) == {"a", "b"}
    assert all(isinstance(v, LeakageStat) for v in stats.values())


def test_report_partitions_the_columns_exactly():
    """clean and leaking must together account for every column, with no overlap."""
    rng, y = _frame()
    x_leak = rng.normal(size=len(y))
    x_leak[
        np.where(y == 1, rng.random(len(y)) < 0.1, rng.random(len(y)) < 0.9)
    ] = np.nan
    X = pd.DataFrame({"leaky": x_leak, "fine": rng.normal(size=len(y))})

    report = resolve_leakage(X, y, ["leaky", "fine"])

    assert set(report.clean) | set(report.leaking) == {"leaky", "fine"}
    assert not (set(report.clean) & set(report.leaking))


def test_selector_is_usable_as_a_nested_selection_step():
    """The control must be runnable per fold on a training slice.

    Computing it once on all rows would let held-out labels choose the feature set, the
    error selection.py's module docstring exists to prevent. The selector signature matches
    what _raw_oof already accepts.
    """
    from hvantk.algorithms.rerank.leakage import leakage_selector

    rng, y = _frame()
    x = rng.normal(size=len(y))
    x[np.where(y == 1, rng.random(len(y)) < 0.1, rng.random(len(y)) < 0.9)] = np.nan
    X = pd.DataFrame({"leaky": x, "fine": rng.normal(size=len(y))})

    sel = leakage_selector()
    kept = list(sel(X, y, ["leaky", "fine"]))

    assert kept == ["fine"]


def test_selector_preserves_input_order_of_kept_columns():
    rng, y = _frame()
    X = pd.DataFrame({c: rng.normal(size=len(y)) for c in ("c", "a", "b")})

    from hvantk.algorithms.rerank.leakage import leakage_selector

    assert list(leakage_selector()(X, y, ["c", "a", "b"])) == ["c", "a", "b"]


def test_empty_column_list_is_not_an_error():
    rng, y = _frame()
    report = resolve_leakage(pd.DataFrame({"a": rng.normal(size=len(y))}), y, [])
    assert report.clean == ()
    assert report.leaking == {}


def test_single_class_label_is_untestable_not_leaking():
    """A degenerate label cannot be predicted by anything, including missingness."""
    rng = np.random.default_rng(1)
    y = np.ones(300, dtype=int)
    x = rng.normal(size=300)
    x[rng.random(300) < 0.5] = np.nan
    report = resolve_leakage(pd.DataFrame({"a": x}), y, ["a"])
    assert "a" in report.clean


@pytest.mark.parametrize("bad_q", [-0.1, 0.0, 1.5])
def test_invalid_fdr_level_is_rejected(bad_q):
    """A q outside (0, 1] silently disables or inverts the filter; refuse it."""
    rng, y = _frame()
    with pytest.raises(ValueError, match="q"):
        resolve_leakage(pd.DataFrame({"a": rng.normal(size=len(y))}), y, ["a"], q=bad_q)
