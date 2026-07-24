# hvantk/tests/rerank/test_tiers.py
import numpy as np, pandas as pd
from hvantk.algorithms.rerank.tiers import TierAssigner


def test_tiers_are_pure_credibility():
    scores = np.linspace(0, 1, 10)
    out = TierAssigner(tiers=5).assign(scores)      # note: no veto argument
    assert set(out.verdict) <= {"ROBUST", "INTERMEDIATE"}
    assert out.loc[9, "verdict"] == "ROBUST"        # highest score -> top tier
    assert out.loc[0, "verdict"] == "INTERMEDIATE"  # lowest score -> low tier
    assert set(out.tier) <= {f"T{i}" for i in range(1, 6)}
    assert "T0_FRAGILE" not in set(out.tier) and "FRAGILE" not in set(out.verdict)
