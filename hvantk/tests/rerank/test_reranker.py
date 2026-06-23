# local/rerank_engine/tests/test_reranker.py
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from hvantk.algorithms.rerank.reranker import ReRanker

def test_oof_scores_are_calibrated_and_discriminate():
    rng = np.random.default_rng(0); n=400
    y = (rng.random(n) < 0.3).astype(int)
    x = y + rng.normal(0, 0.5, n)            # informative feature
    m = pd.DataFrame({"gene":[f"g{i}" for i in range(n)], "x": x})
    p = ReRanker(folds=5).score(m, ["x"], y)
    assert p.shape == (n,) and ((p>=0)&(p<=1)).all()
    assert roc_auc_score(y, p) > 0.8         # informative feature -> good OOF AUC
