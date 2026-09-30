"""A byte-for-byte pin on evaluator.py, captured BEFORE the #247 reformat.

`evaluator.py` opened with the stale line `# local/rerank_engine/evaluator.py` and packed
several statements onto most lines. It is the module #247 edits most, so it was reformatted
first and separately -- but "no behaviour change" is a claim, and a claim about a stochastic
estimator needs a fixed seed and recorded numbers rather than a reading of the diff. These
values come from the unmodified file at dev@693c490c; every later change to this module must
leave them untouched too, which is why this file stays in the rerank test selection every
time evaluator.py is touched.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from hvantk.algorithms.rerank.evaluator import Evaluator, _boot_ci, _raw_oof

# --- captured from dev@693c490c, via a fixed-seed run of the fixture below --------------
OOF_HEAD = [
    0.294651741677,
    0.164421721731,
    0.444986827035,
    0.930760857286,
    0.89923815832,
    0.485374629101,
    0.93073060035,
    0.134549842236,
]  # list[float], 8 values
OOF_SUM = 70.790331372377  # float
BOOT_CI = [-0.027, 0.009, 0.04]  # list[float], 3 values
ABLATION = [
    {"family": "b", "auc": 0.819, "d_lo": 0.0, "d_md": 0.0, "d_hi": 0.0},
    {"family": "e", "auc": 0.829, "d_lo": -0.027, "d_md": 0.01, "d_hi": 0.041},
]  # list[dict]


def _fixture():
    rng = np.random.default_rng(11)
    n = 160
    y = (rng.random(n) < 0.3).astype(int)
    m = pd.DataFrame(
        {
            "gene": [f"g{i}" for i in range(n)],
            "base": y + rng.normal(0, 0.8, n),
            "extra": rng.normal(0, 1, n),
        }
    )
    return m, y


def test_raw_oof_is_unchanged():
    m, y = _fixture()
    p = _raw_oof(m, ["base"], y)
    assert [round(float(v), 12) for v in p[:8]] == OOF_HEAD
    assert round(float(p.sum()), 12) == OOF_SUM


def test_boot_ci_is_unchanged():
    m, y = _fixture()
    p0 = _raw_oof(m, ["base"], y)
    p1 = _raw_oof(m, ["base", "extra"], y)
    assert _boot_ci(y, p1, p0, n=200) == BOOT_CI


def test_ablation_table_is_unchanged():
    m, y = _fixture()
    p = _raw_oof(m, ["base"], y)
    ev = Evaluator().evaluate(
        m, ["base", "extra"], y, p, {"b": ["base"], "e": ["extra"]}, "b"
    )
    assert ev.ablation.to_dict("records") == ABLATION


def test_the_prototype_header_is_gone():
    """The file was moved out of a prototype tree and the path comment survived the move."""
    from pathlib import Path

    import hvantk.algorithms.rerank.evaluator as mod

    first = Path(mod.__file__).read_text().splitlines()[0]
    assert "local/rerank_engine" not in first, first
    assert mod.__doc__, "the module needs a docstring, not a path comment"
