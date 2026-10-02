"""The multiplicity correction proper, plus the two ways it is silently misused.

A null is generated under ONE control setting and covers ONE candidate set. Attaching it to
deltas computed under a different feature set, or merging chunks that disagree, produces a
plausible number that answers a question nobody asked.
"""

from __future__ import annotations

import dataclasses
import functools

import numpy as np
import pytest

from hvantk.algorithms.rerank.evaluator import ABLATION_FOLDS
from hvantk.algorithms.rerank.leakage import LeakagePolicy
from hvantk.algorithms.rerank.nulls import (
    ControlSetting,
    ControlSettingMismatch,
    NullConfig,
    NullDistribution,
    axis_deltas,
    permutation_deltas,
    selected_maximum,
)
from hvantk.algorithms.rerank.selection import SelectionPolicy
from hvantk.tests.rerank._synth import cheap_scorer, permuted_labels, planted_signal


def _setting(candidates, **kw):
    base = dict(
        arm="all",
        leakage=None,
        selection=None,
        baseline=("base",),
        candidates=candidates,
        folds=ABLATION_FOLDS,
        block_digest=None,
    )
    base.update(kw)
    return ControlSetting(**base)


@functools.lru_cache(maxsize=None)
def _null(n_axes, n_perm=30, seed=4, n=160):
    """Built repeatedly across the guard/merge tests with the SAME arguments -- a
    `NullDistribution` is frozen, so caching by argument tuple is safe, and every caller
    that needs a different shape just gets a fresh (also cached) build under its own key."""
    matrix, y, baseline, axes = permuted_labels(n=n, n_noise=5)
    offered = {k: axes[k] for k in sorted(axes)[:n_axes]}
    cfg = NullConfig(n_perm=n_perm, seed=seed)
    d = permutation_deltas(
        matrix, baseline, offered, y, config=cfg, scorer=cheap_scorer()
    )
    return NullDistribution.from_deltas(d, _setting(offered), null_config=cfg)


# --- the selected maximum -----------------------------------------------------------------


def test_selected_maximum_is_the_per_permutation_max_over_axes():
    import pandas as pd

    d = pd.DataFrame(
        {
            "perm": [0, 0, 0, 1, 1, 1],
            "axis": ["a", "b", "c"] * 2,
            "delta": [0.01, 0.05, np.nan, -0.02, 0.00, 0.03],
            "base_auc": [0.6] * 6,
        }
    )
    assert selected_maximum(d).tolist() == [0.05, 0.03]


def test_a_planted_signal_clears_the_selected_maximum_null():
    from sklearn.metrics import roc_auc_score

    matrix, y, baseline, axes = planted_signal(n=200, n_noise=3)
    scorer = cheap_scorer()
    a0 = roc_auc_score(y, scorer(matrix, baseline, y))
    a1 = roc_auc_score(y, scorer(matrix, baseline + axes["axis0"], y))
    cfg = NullConfig(n_perm=39, seed=1)
    d = permutation_deltas(matrix, baseline, axes, y, config=cfg, scorer=scorer)
    setting = _setting(axes)
    nd = NullDistribution.from_deltas(d, setting, null_config=cfg)
    p = nd.p_selected_max(a1 - a0, setting=setting)
    assert p <= 0.05, (p, a1 - a0, float(np.median(nd.selected_max)))


def test_selected_max_is_never_below_the_per_axis_draw_it_contains():
    nd = _null(n_axes=3)
    for axis, draws in nd.per_axis.items():
        assert np.all(nd.selected_max >= np.nan_to_num(draws, nan=-np.inf)), axis


# --- the control-setting guard ------------------------------------------------------------


@pytest.mark.parametrize(
    "different",
    [
        dict(leakage=LeakagePolicy()),
        dict(selection=SelectionPolicy()),
        dict(arm="clean"),
        dict(baseline=("base", "burden")),
        dict(block_digest="0" * 40),
        dict(folds=10),
        dict(candidates={"axis0": ("elsewhere",)}),
    ],
)
def test_a_null_refuses_a_delta_from_another_control_setting(different):
    """Attaching a null generated under one control setting to a delta computed under
    another is a category error: the two quantities are not the same statistic. It must
    raise, not silently mismatch."""
    nd = _null(n_axes=2)
    other = dataclasses.replace(nd.setting, **different)
    with pytest.raises(ControlSettingMismatch) as info:
        nd.p_selected_max(0.02, setting=other)
    msg = str(info.value)
    assert "generated under" in msg and "asked about" in msg
    # `describe()` prints every field name, so asserting the changed field's name is
    # anywhere in `msg` can never fail -- it must show up specifically in the differing-
    # field list, and an unrelated field must NOT. A bare substring search on just the
    # field name is not enough either: e.g. "arm" is itself a substring of "spearman",
    # which appears inside SelectionPolicy's own repr. Anchor on the exact
    # "<field> (generated=" marker `_setting_mismatch` emits for each differing field.
    changed = next(iter(different))
    differing = msg.split("Differing field(s): ", 1)[1]
    assert f"{changed} (generated=" in differing
    unchanged = "arm" if changed != "arm" else "folds"
    assert f"{unchanged} (generated=" not in differing


def test_the_guard_also_covers_the_per_axis_question():
    nd = _null(n_axes=2)
    other = dataclasses.replace(nd.setting, leakage=LeakagePolicy())
    with pytest.raises(ControlSettingMismatch):
        nd.p_per_axis("axis0", 0.02, setting=other)


def test_a_matching_setting_is_accepted():
    nd = _null(n_axes=2)
    assert 0.0 < nd.p_selected_max(0.02, setting=nd.setting) <= 1.0


def test_an_unknown_axis_is_a_key_error_naming_what_was_offered():
    nd = _null(n_axes=2)
    with pytest.raises(KeyError, match="axis0"):
        nd.p_per_axis("never_offered", 0.02, setting=nd.setting)


# --- chunk merging --------------------------------------------------------------------------


def test_merged_chunks_equal_the_whole_run():
    matrix, y, baseline, axes = permuted_labels(n=140, n_noise=2)
    offered = {k: axes[k] for k in ("axis0", "axis1", "axis2")}
    scorer = cheap_scorer()
    setting = _setting(offered)

    whole_cfg = NullConfig(n_perm=6, seed=8)
    whole = NullDistribution.from_deltas(
        permutation_deltas(
            matrix, baseline, offered, y, config=whole_cfg, scorer=scorer
        ),
        setting,
        null_config=whole_cfg,
    )

    parts = []
    for c in range(3):
        cfg = NullConfig(n_perm=6, chunk=c, n_chunks=3, seed=8)
        deltas = permutation_deltas(
            matrix, baseline, offered, y, config=cfg, scorer=scorer
        )
        parts.append(NullDistribution.from_deltas(deltas, setting, null_config=cfg))

    merged = NullDistribution.merge(parts)
    assert merged.n_perm == whole.n_perm == 6
    assert np.allclose(np.sort(merged.selected_max), np.sort(whole.selected_max))


def test_merge_refuses_chunks_from_different_control_settings():
    a = _null(n_axes=2, seed=1)
    b = _null(n_axes=2, seed=1)
    b = dataclasses.replace(
        b, setting=dataclasses.replace(b.setting, leakage=LeakagePolicy())
    )
    with pytest.raises(ControlSettingMismatch):
        NullDistribution.merge([a, b])


def test_merge_refuses_overlapping_permutation_indices():
    """Two runs of chunk 0 are the SAME permutations. Concatenating them would double-count
    them and narrow every p-value for free."""
    a = _null(n_axes=2, n_perm=4, seed=1)
    b = _null(n_axes=2, n_perm=4, seed=1)
    with pytest.raises(ValueError, match="overlap"):
        NullDistribution.merge([a, b])

    # A chunk concatenated into its own deltas frame leaves the perm SET and axis SET
    # unchanged (both checks above would miss it), yet double-counts every draw.
    import pandas as pd

    matrix, y, baseline, axes = permuted_labels(n=60, n_noise=1)
    offered = {"axis0": axes["axis0"]}
    cfg = NullConfig(n_perm=2, seed=1)
    d = permutation_deltas(
        matrix, baseline, offered, y, config=cfg, scorer=cheap_scorer()
    )
    with pytest.raises(ValueError, match="duplicate"):
        NullDistribution.from_deltas(
            pd.concat([d, d]), _setting(offered), null_config=cfg
        )


def test_merge_refuses_chunks_that_offered_different_axes():
    with pytest.raises(ValueError, match="axes"):
        NullDistribution.merge(
            [_null(n_axes=2, n_perm=4), _null(n_axes=3, n_perm=4, seed=9)]
        )


def test_merge_of_one_is_that_one():
    a = _null(n_axes=2, n_perm=4)
    assert NullDistribution.merge([a]).n_perm == a.n_perm


def test_merge_refuses_chunks_from_different_permutation_seeds():
    """Two chunks whose BASE seeds differ can draw the same underlying permutations from
    disjoint chunk indices, which the plain perm-index overlap check cannot see."""
    matrix, y, baseline, axes = permuted_labels(n=60, n_noise=1)
    offered = {"axis0": axes["axis0"]}
    scorer = cheap_scorer()
    setting = _setting(offered)

    cfg_a = NullConfig(40, chunk=1, n_chunks=2, seed=42)
    cfg_b = NullConfig(40, chunk=0, n_chunks=2, seed=52)
    a = NullDistribution.from_deltas(
        permutation_deltas(matrix, baseline, offered, y, config=cfg_a, scorer=scorer),
        setting,
        null_config=cfg_a,
    )
    b = NullDistribution.from_deltas(
        permutation_deltas(matrix, baseline, offered, y, config=cfg_b, scorer=scorer),
        setting,
        null_config=cfg_b,
    )
    with pytest.raises(ValueError, match="seed"):
        NullDistribution.merge([a, b])


def test_observed_deltas_ride_with_the_null_and_must_agree_across_chunks():
    matrix, y, baseline, axes = permuted_labels(n=140, n_noise=2)
    offered = {k: axes[k] for k in ("axis0", "axis1", "axis2")}
    scorer = cheap_scorer()
    setting = _setting(offered)
    _, observed = axis_deltas(matrix, baseline, offered, y, scorer=scorer)

    parts = []
    for c in range(2):
        cfg = NullConfig(n_perm=6, chunk=c, n_chunks=2, seed=8)
        deltas = permutation_deltas(
            matrix, baseline, offered, y, config=cfg, scorer=scorer
        )
        parts.append(
            NullDistribution.from_deltas(
                deltas, setting, null_config=cfg, observed=observed
            )
        )

    merged = NullDistribution.merge(parts)
    assert merged.observed == observed

    summary = merged.summary().set_index("axis")
    for axis in merged.axes:
        row = summary.loc[axis]
        assert row["observed"] == pytest.approx(observed[axis])
        assert row["p_per_axis"] == pytest.approx(
            merged.p_per_axis(axis, observed[axis], setting=setting)
        )
        assert row["p_selected_max"] == pytest.approx(
            merged.p_selected_max(observed[axis], setting=setting)
        )

    other_cfg = NullConfig(n_perm=6, chunk=0, n_chunks=2, seed=8)
    other_deltas = permutation_deltas(
        matrix, baseline, offered, y, config=other_cfg, scorer=scorer
    )
    first_axis = next(iter(observed))

    different_observed = dict(observed)
    different_observed[first_axis] += 1.0
    bad = NullDistribution.from_deltas(
        other_deltas, setting, null_config=other_cfg, observed=different_observed
    )
    with pytest.raises(ValueError, match="observed"):
        NullDistribution.merge([parts[0], bad])

    no_obs = NullDistribution.from_deltas(other_deltas, setting, null_config=other_cfg)
    with pytest.raises(ValueError, match="observed"):
        NullDistribution.merge([parts[0], no_obs])

    nan_observed = dict(observed)
    nan_observed[first_axis] = float("nan")
    nd_nan = NullDistribution.from_deltas(
        other_deltas, setting, null_config=other_cfg, observed=nan_observed
    )
    nan_row = nd_nan.summary().set_index("axis").loc[first_axis]
    assert np.isnan(nan_row["observed"])
    assert np.isnan(nan_row["p_per_axis"])
    assert np.isnan(nan_row["p_selected_max"])


# --- config wiring ---------------------------------------------------------------------------


def test_config_defaults_to_no_null():
    from hvantk.algorithms.rerank.config import Config

    assert Config.__dataclass_fields__["nulls"].default is None


def test_config_rejects_a_non_nullconfig():
    """A bare int here would be silently ignored and the correction would be off while the
    caller believed it was on -- the failure mode Config.leakage's type check exists for."""
    from hvantk.algorithms.rerank.config import Config

    with pytest.raises(TypeError, match="NullConfig"):
        Config(name="x", features=[], labels=None, nulls=200).__post_init__()


def test_a_requested_null_with_no_candidate_axis_is_an_error_not_a_none():
    """`RerankResult.nulls is None` means "no correction was requested". The engine used to
    log a warning and return exactly that when a correction WAS requested but no axis
    besides the baseline had columns left -- the state Config.nulls's type check calls
    worse than never offering a correction, now reachable through the API with nothing
    to show for it. It must raise, naming the configured axes, the one that was dropped
    and why, and how to run without a correction."""
    import pandas as pd

    from hvantk.algorithms.rerank.config import Config, FeatureAxis, LabelSpec
    from hvantk.algorithms.rerank.engine import _run_nulls

    frame = pd.DataFrame({"gene": ["A", "B"], "x": [0.1, 0.2]})
    empty = pd.DataFrame({"gene": ["A", "B"]})  # a table with no feature column
    cfg = Config(
        name="x",
        features=[
            FeatureAxis("constraint", lambda: frame),
            FeatureAxis("expression", lambda: empty),
        ],
        labels=LabelSpec(lambda: {"A"}),
        nulls=NullConfig(n_perm=2),
    )

    with pytest.raises(ValueError) as info:
        _run_nulls(
            cfg, frame, "constraint", {"constraint": ["x"]}, [1, 0], None, "all", None
        )
    msg = str(info.value)
    assert "Config.nulls" in msg and "'constraint'" in msg
    assert "'expression'" in msg and "no feature column" in msg
    assert "unset" in msg
