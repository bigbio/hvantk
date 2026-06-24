# hvantk/tests/rerank/test_smoke.py
import numpy as np, pandas as pd
from hvantk.algorithms.rerank import rerank, Config, PriorSpec, FeatureAxis, LabelSpec


def test_engine_end_to_end(tmp_path):
    rng = np.random.default_rng(0); n = 200
    genes = [f"g{i}" for i in range(n)]; y = (rng.random(n) < 0.3).astype(int)
    feat = pd.DataFrame({"gene": genes, "x": y + rng.normal(0, 0.5, n)})
    prior = pd.DataFrame({"gene": genes, "p": rng.random(n)}); pf = tmp_path/"p.tsv"
    prior.to_csv(pf, sep="\t", index=False)
    pos = {g for g, yy in zip(genes, y) if yy == 1}
    cfg = Config(name="t", prior=PriorSpec(str(pf), "gene", "p"),
                 features=[FeatureAxis("x", lambda: feat)], labels=LabelSpec(lambda: pos),
                 min_label_coverage=0.0)
    res = rerank(cfg)
    assert len(res.table) == n and res.metrics.auc > 0.7
    assert {"gene", "score", "tier", "verdict"} <= set(res.table.columns)
