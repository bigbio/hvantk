# local/rerank_engine/tiers.py
import numpy as np, pandas as pd

class TierAssigner:
    def __init__(self, tiers=5): self.tiers = tiers
    def assign(self, scores, veto_mask):
        s = pd.Series(np.asarray(scores)); veto = pd.Series(np.asarray(veto_mask)).reset_index(drop=True)
        labels = [f"T{i+1}" for i in range(self.tiers)]
        tier = pd.qcut(s.rank(method="first"), self.tiers, labels=labels).astype(str)
        final_tier = np.where(veto, "T0_FRAGILE", tier)
        top2 = set(labels[-2:])
        verdict = np.where(veto, "FRAGILE", np.where(pd.Series(final_tier).isin(top2), "ROBUST", "INTERMEDIATE"))
        return pd.DataFrame({"tier": final_tier, "verdict": verdict})
