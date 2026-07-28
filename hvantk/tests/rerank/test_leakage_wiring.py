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
