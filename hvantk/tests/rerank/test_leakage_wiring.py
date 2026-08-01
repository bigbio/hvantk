"""Wiring the presence-leakage control into Config and the engine's selector chain.

The unit behaviour lives in test_leakage.py. This covers the part that makes it reachable:
a control nobody can switch on does not close the gap it was written for.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hvantk.algorithms.rerank.config import Config
from hvantk.algorithms.rerank.engine import _compose_selector
from hvantk.algorithms.rerank.leakage import LeakagePolicy


def _leaky_frame(n=600, seed=0):
    """A matrix with one leaking column, one clean column, and a baseline axis."""
    rng = np.random.default_rng(seed)
    y = np.zeros(n, dtype=int)
    y[: n // 5] = 1
    rng.shuffle(y)
    leaky = rng.normal(size=n)
    keep = np.where(y == 1, rng.random(n) < 0.9, rng.random(n) < 0.15)
    leaky[~keep] = np.nan
    X = pd.DataFrame(
        {
            "pLI": rng.normal(size=n),
            "eve_mean": leaky,
            "expr": rng.normal(size=n),
        }
    )
    return X, y


def test_config_defaults_to_no_leakage_control():
    """Default must reproduce the pre-existing code path exactly, like `selection`."""
    assert Config.__dataclass_fields__["leakage"].default is None


def test_compose_selector_returns_none_when_nothing_is_configured():
    """No selection and no leakage means no selector at all -- not an identity wrapper.

    ReRanker.score takes a different, cheaper path when selector is None, and an identity
    selector would silently change which code path every existing run takes.
    """
    assert _compose_selector(selection_selector=None, leakage_policy=None) is None


def test_leakage_alone_produces_a_working_selector():
    """Leakage control must not require a SelectionPolicy.

    They are independent controls: one asks whether a column's VALUES are worth keeping,
    the other whether its MISSINGNESS carries the label.
    """
    X, y = _leaky_frame()
    sel = _compose_selector(selection_selector=None, leakage_policy=LeakagePolicy())

    kept = list(sel(X, y, ["pLI", "eve_mean", "expr"]))

    assert "eve_mean" not in kept
    assert kept == ["pLI", "expr"]


def test_leakage_runs_before_axis_selection():
    """Order matters: a leaking column must never reach the univariate filter.

    That filter scores |AUC-0.5| on the column's VALUES and would rank a leaking column
    highly -- it rewards the defect, which is the whole reason the control cannot be left
    to statistics.
    """
    X, y = _leaky_frame()
    seen: list[list[str]] = []

    def spy(X_tr, y_tr, columns):
        seen.append(list(columns))
        return list(columns)

    sel = _compose_selector(selection_selector=spy, leakage_policy=LeakagePolicy())
    sel(X, y, ["pLI", "eve_mean", "expr"])

    assert seen, "the downstream selector must still be called"
    assert "eve_mean" not in seen[0], "leaking column reached axis selection"


def test_selection_alone_is_unchanged_by_this_feature():
    """With no leakage policy the composed selector must be the original object.

    Wrapping it would change behaviour for every existing caller for no reason.
    """

    def spy(X_tr, y_tr, columns):
        return list(columns)

    assert _compose_selector(selection_selector=spy, leakage_policy=None) is spy


def test_composed_selector_keeps_caller_column_order():
    X, y = _leaky_frame()
    sel = _compose_selector(selection_selector=None, leakage_policy=LeakagePolicy())
    assert list(sel(X, y, ["expr", "pLI"])) == ["expr", "pLI"]


def test_policy_thresholds_are_honoured():
    """A permissive effect floor must let a mildly-leaking column through, so the knob is
    demonstrably wired rather than merely present."""
    X, y = _leaky_frame()
    strict = _compose_selector(None, LeakagePolicy(min_auc=0.55))
    permissive = _compose_selector(None, LeakagePolicy(min_auc=0.99))

    assert "eve_mean" not in list(strict(X, y, ["pLI", "eve_mean", "expr"]))
    assert "eve_mean" in list(permissive(X, y, ["pLI", "eve_mean", "expr"]))


def test_config_rejects_an_invalid_policy_type():
    """A bare float or dict here would be silently ignored by the engine."""
    with pytest.raises((TypeError, ValueError)):
        Config(name="x", features=[], labels=None, leakage=0.55).__post_init__()


def test_global_selection_summary_excludes_leaking_columns():
    """The global summary must not report a column the nested folds barred.

    `_selection_summary` runs its own selection over ALL the data to produce the
    human-readable "these are the features" list. That pass is separate from the per-fold
    selector, so without explicit filtering a leaking column can appear in
    `global_features` -- and be read as endorsed -- while no fold ever used it.
    """
    from hvantk.algorithms.rerank.engine import _leakage_filtered_groups

    X, y = _leaky_frame()
    groups = {"constraint": ["pLI"], "dbnsfp": ["eve_mean"], "expression": ["expr"]}

    filtered = _leakage_filtered_groups(X, y, groups, LeakagePolicy())

    assert "eve_mean" not in filtered.get("dbnsfp", [])
    assert filtered["constraint"] == ["pLI"]
    assert filtered["expression"] == ["expr"]


def test_leakage_filtered_groups_drops_an_axis_left_empty():
    """An axis whose every column leaks must disappear, not survive as an empty list.

    A zero-column axis reaching select_axis is reported as present-but-contributing-
    nothing, which reads as a measured null rather than a barred axis.
    """
    from hvantk.algorithms.rerank.engine import _leakage_filtered_groups

    X, y = _leaky_frame()
    filtered = _leakage_filtered_groups(
        X, y, {"only_leaky": ["eve_mean"], "fine": ["expr"]}, LeakagePolicy()
    )

    assert "only_leaky" not in filtered
    assert filtered["fine"] == ["expr"]


def test_leakage_filtered_groups_is_identity_without_a_policy():
    from hvantk.algorithms.rerank.engine import _leakage_filtered_groups

    X, y = _leaky_frame()
    groups = {"a": ["pLI"], "b": ["eve_mean"]}
    assert _leakage_filtered_groups(X, y, groups, None) is groups


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -0.2, 1.4])
def test_policy_rejects_non_finite_or_out_of_range_min_auc(bad):
    """NaN silently disables the effect floor: max(0.0, nan) is 0.0, so every
    FDR-significant column would be barred with no effect-size protection at all."""
    with pytest.raises(ValueError, match="min_auc"):
        Config(
            name="x", features=[], labels=None, leakage=LeakagePolicy(min_auc=bad)
        ).__post_init__()


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_config_rejects_invalid_policy_q(bad):
    """presence_leakage raises for these, but only once a per-fold selector reaches it."""
    with pytest.raises(ValueError, match="q"):
        Config(
            name="x", features=[], labels=None, leakage=LeakagePolicy(q=bad)
        ).__post_init__()
