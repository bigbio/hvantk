# local/rerank_engine/tests/test_evaluator.py
import numpy as np, pandas as pd
from hvantk.algorithms.rerank.evaluator import Evaluator

def test_metrics_and_ablation_shapes():
    rng = np.random.default_rng(0); n=300
    y = (rng.random(n)<0.3).astype(int)
    base = y + rng.normal(0,0.5,n); extra = rng.normal(0,1,n)   # extra is noise
    m = pd.DataFrame({"gene":[f"g{i}" for i in range(n)], "base":base, "extra":extra})
    scores = (base - base.min())/(base.max()-base.min())
    ev = Evaluator().evaluate(m, ["base","extra"], y, scores,
            axis_groups={"baseline":["base"], "extra":["extra"]}, baseline_axis="baseline")
    assert 0.5 < ev.auc <= 1.0 and 0.0 < ev.brier < 0.5
    assert "extra" in set(ev.ablation.family) and {"d_lo","d_md","d_hi"} <= set(ev.ablation.columns)
