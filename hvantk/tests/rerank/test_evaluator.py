import numpy as np, pandas as pd
from hvantk.algorithms.rerank.evaluator import Evaluator, _boot_ci, _raw_oof
from hvantk.tests.rerank._synth import planted_signal


def test_metrics_and_ablation_shapes():
    rng = np.random.default_rng(0)
    n = 300
    y = (rng.random(n) < 0.3).astype(int)
    base = y + rng.normal(0, 0.5, n)
    extra = rng.normal(0, 1, n)  # extra is noise
    m = pd.DataFrame(
        {"gene": [f"g{i}" for i in range(n)], "base": base, "extra": extra}
    )
    scores = (base - base.min()) / (base.max() - base.min())
    ev = Evaluator().evaluate(
        m,
        ["base", "extra"],
        y,
        scores,
        axis_groups={"baseline": ["base"], "extra": ["extra"]},
        baseline_axis="baseline",
    )
    assert 0.5 < ev.auc <= 1.0 and 0.0 < ev.brier < 0.5
    assert "extra" in set(ev.ablation.family) and {"d_lo", "d_md", "d_hi"} <= set(
        ev.ablation.columns
    )


def test_evaluate_uses_the_block_bootstrap_when_groups_is_given():
    """Guards the one line that forwards `groups` from `evaluate` into `_boot_ci`'s
    per-axis call: dropping it there (reverting to the gene-only bootstrap even when the
    run is blocked) passes every other rerank test and would only fail here."""
    # 30 families of 6, matching test_blocked_cv.py's own fixture.
    matrix, y, baseline, axes = planted_signal(n=180, n_noise=1, seed=5)
    blocks = np.arange(180) // 6
    cols = baseline + axes["axis0"]

    p_base = _raw_oof(matrix, baseline, y, groups=blocks)
    p_aug = _raw_oof(matrix, cols, y, groups=blocks)
    row = (
        Evaluator()
        .evaluate(
            matrix,
            cols,
            y,
            p_base,
            {"base": baseline, "axis0": axes["axis0"]},
            "base",
            groups=blocks,
        )
        .ablation.set_index("family")
        .loc["axis0"]
    )
    blocked = _boot_ci(y, p_aug, p_base, groups=blocks)
    unblocked = _boot_ci(y, p_aug, p_base)
    assert [row.d_lo, row.d_md, row.d_hi] == blocked
    assert [row.d_lo, row.d_md, row.d_hi] != unblocked
