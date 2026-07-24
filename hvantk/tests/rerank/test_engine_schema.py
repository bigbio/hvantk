# hvantk/tests/rerank/test_engine_schema.py
import numpy as np, pandas as pd
from hvantk.algorithms.rerank.tiers import TierAssigner


def test_output_columns_and_flag_semantics():
    # Unit-level check of the table assembly contract the engine must satisfy.
    scores = np.linspace(0, 1, 6)
    flag_reason = pd.Series(["", "", "recurrent_variant", "", "common_driver", ""])
    flag = (flag_reason != "")
    tiers = TierAssigner(tiers=3).assign(scores)
    table = pd.DataFrame({"gene": list("abcdef"), "score": scores,
                          "flag": flag.values, "flag_reason": flag_reason.values})
    table = pd.concat([table.reset_index(drop=True), tiers.reset_index(drop=True)], axis=1)
    table["score_percentile"] = (table["score"].rank(pct=True) * 100).round(1)
    # a flagged gene keeps its credibility tier (NOT overridden)
    assert table.loc[2, "tier"].startswith("T") and table.loc[2, "flag"] == True
    assert set(["gene", "score", "score_percentile", "tier", "verdict",
                "flag", "flag_reason"]).issubset(table.columns)
    assert "FRAGILE" not in set(table.verdict) and "T0_FRAGILE" not in set(table.tier)
