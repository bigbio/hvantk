# local/rerank_engine/tests/test_tiers.py
import numpy as np, pandas as pd
from hvantk.algorithms.rerank.tiers import TierAssigner

def test_tiers_and_veto_override():
    scores = np.linspace(0, 1, 10)
    veto = pd.Series([False]*9 + [True])      # highest-scoring unit is vetoed
    out = TierAssigner(tiers=5).assign(scores, veto)
    assert out.loc[9, "tier"] == "T0_FRAGILE" and out.loc[9, "verdict"] == "FRAGILE"
    assert out.loc[0, "verdict"] == "INTERMEDIATE"   # lowest non-vetoed -> low tier, not ROBUST
    assert set(out.verdict) <= {"FRAGILE","INTERMEDIATE","ROBUST"}
    assert "ROBUST" in set(out.verdict)
