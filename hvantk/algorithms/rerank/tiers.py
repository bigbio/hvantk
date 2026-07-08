# hvantk/algorithms/rerank/tiers.py
import numpy as np, pandas as pd


class TierAssigner:
    """Pure credibility tiers over scored genes. The audit flag is a SEPARATE
    axis (see engine.py) and is deliberately NOT folded into the tier ladder."""
    def __init__(self, tiers=5):
        self.tiers = tiers

    def assign(self, scores) -> pd.DataFrame:
        s = pd.Series(np.asarray(scores))
        labels = [f"T{i + 1}" for i in range(self.tiers)]
        tier = pd.qcut(s.rank(method="first"), self.tiers, labels=labels).astype(str)
        top2 = set(labels[-2:])
        verdict = np.where(pd.Series(tier).isin(top2), "ROBUST", "INTERMEDIATE")
        return pd.DataFrame({"tier": tier, "verdict": verdict})
