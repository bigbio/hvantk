"""The CV partition as a measured variance component, not a fixed choice.

`_boot_ci` resamples GENES. Which genes land in which fold is a second source of variation,
and the shipped interval could not see it at all.
"""
from __future__ import annotations

import numpy as np
import pytest

from hvantk.algorithms.rerank.evaluator import (
    Evaluator,
    SeedSpread,
    _envelope,
    _raw_oof,
    seed_sweep,
)
from hvantk.algorithms.rerank.seeds import DEFAULT_SEED
from hvantk.tests.rerank._synth import planted_signal


def _fixture(n=180):
    matrix, y, baseline, axes = planted_signal(n=n, n_noise=1)
    return matrix, y, baseline, axes["axis0"]


def test_a_sweep_reports_one_delta_per_seed():
    matrix, y, base, axis = _fixture()
    s = seed_sweep(matrix, base, axis, y, seeds=(1, 2, 3))
    assert isinstance(s, SeedSpread)
    assert s.seeds == (1, 2, 3) and len(s.deltas) == 3
    assert s.lo <= s.hi
    assert s.mean == pytest.approx(float(np.mean(s.deltas)))


def test_different_seeds_give_different_deltas():
    """If they did not, the sweep would be reporting a constant and the widening below
    would be an artefact of the code rather than a measurement."""
    matrix, y, base, axis = _fixture()
    s = seed_sweep(matrix, base, axis, y, seeds=(1, 2, 3, 4))
    assert len(set(round(d, 9) for d in s.deltas)) > 1


def test_seed_zero_offset_matches_the_headline():
    """Guard the guard: a sweep that disagreed with the shipped estimator at the default
    seed would invalidate every comparison built on it."""
    from sklearn.metrics import roc_auc_score

    matrix, y, base, axis = _fixture()
    s = seed_sweep(matrix, base, axis, y, seeds=(DEFAULT_SEED,))
    expected = roc_auc_score(y, _raw_oof(matrix, base + axis, y)) - roc_auc_score(
        y, _raw_oof(matrix, base, y)
    )
    assert s.deltas[0] == pytest.approx(expected)


def test_a_sweep_is_reproducible():
    matrix, y, base, axis = _fixture()
    a = seed_sweep(matrix, base, axis, y, seeds=(1, 2))
    b = seed_sweep(matrix, base, axis, y, seeds=(1, 2))
    assert a.deltas == b.deltas


def test_the_envelope_can_only_widen():
    boot = (-0.01, 0.05)
    assert _envelope(boot, SeedSpread((1, 2), (0.0, 0.02), 0.0, 0.02, 0.01, 0.01)) == boot
    wider = _envelope(boot, SeedSpread((1, 2), (-0.03, 0.09), -0.03, 0.09, 0.03, 0.06))
    assert wider[0] <= boot[0] and wider[1] >= boot[1]
    assert wider == (-0.03, 0.09)


def test_the_envelope_of_no_sweep_is_the_bootstrap_interval():
    assert _envelope((-0.01, 0.05), None) == (-0.01, 0.05)


# --- the reported interval widens -----------------------------------------------------------


def test_the_sweep_envelope_is_wired_through_evaluate():
    """The ablation's envelope columns are exactly the rounded union of the bootstrap interval
    and the across-seed range, and never narrower than the bootstrap interval (a strict
    widening depends on the data and is not guaranteed)."""
    matrix, y, base, axis = _fixture()
    groups = {"base": base, "axis0": axis}
    scores = _raw_oof(matrix, base, y)
    single = Evaluator().evaluate(matrix, base + axis, y, scores, groups, "base")
    swept = Evaluator().evaluate(matrix, base + axis, y, scores, groups, "base", n_seeds=5)

    row_single = single.ablation.set_index("family").loc["axis0"]
    row_swept = swept.ablation.set_index("family").loc["axis0"]
    spread = swept.seed_spread["axis0"]
    env_lo, env_hi = _envelope((row_single.d_lo, row_single.d_hi), spread)
    assert (row_swept.d_lo_env, row_swept.d_hi_env) == (round(env_lo, 3), round(env_hi, 3))
    assert row_swept.d_lo_env <= row_single.d_lo and row_swept.d_hi_env >= row_single.d_hi


def test_the_bootstrap_columns_are_untouched_by_the_sweep():
    """d_lo/d_md/d_hi keep meaning exactly what they meant: the gene-resampling interval.
    The sweep adds columns beside them rather than redefining them."""
    matrix, y, base, axis = _fixture()
    groups = {"base": base, "axis0": axis}
    scores = _raw_oof(matrix, base, y)
    single = Evaluator().evaluate(matrix, base + axis, y, scores, groups, "base")
    swept = Evaluator().evaluate(matrix, base + axis, y, scores, groups, "base", n_seeds=5)
    for col in ("d_lo", "d_md", "d_hi"):
        assert single.ablation[col].tolist() == swept.ablation[col].tolist()


def test_the_default_run_adds_no_sweep_columns():
    """`n_seeds=1` (the default) must reproduce today's ablation table exactly: no
    `d_lo_env`/`d_hi_env`/`n_seeds` columns, and no seed_spread. `test_evaluator_pin.py`
    pins the exact records this default run produces; this test pins the shape."""
    matrix, y, base, axis = _fixture()
    groups = {"base": base, "axis0": axis}
    ev = Evaluator().evaluate(matrix, base + axis, y, _raw_oof(matrix, base, y), groups, "base")
    assert not {"d_lo_env", "d_hi_env", "n_seeds"} & set(ev.ablation.columns)
    assert ev.seed_spread is None


def test_the_spread_is_attached_per_axis():
    matrix, y, base, axis = _fixture()
    groups = {"base": base, "axis0": axis}
    ev = Evaluator().evaluate(
        matrix, base + axis, y, _raw_oof(matrix, base, y), groups, "base", n_seeds=3
    )
    assert set(ev.seed_spread) == {"axis0"}
    assert ev.seed_spread["axis0"].seeds == (DEFAULT_SEED, DEFAULT_SEED + 1, DEFAULT_SEED + 2)


# --- config wiring ---------------------------------------------------------------------------


def test_config_seed_sweep_defaults_to_one():
    from hvantk.algorithms.rerank.config import Config

    assert Config.__dataclass_fields__["seed_sweep"].default == 1


@pytest.mark.parametrize("bad", [0, -3, 2.5, "5"])
def test_config_rejects_an_impossible_sweep(bad):
    from hvantk.algorithms.rerank.config import Config

    with pytest.raises((TypeError, ValueError), match="seed_sweep"):
        Config(name="x", features=[], labels=None, seed_sweep=bad).__post_init__()
