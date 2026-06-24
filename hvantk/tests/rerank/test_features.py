# local/rerank_engine/tests/test_features.py
import pandas as pd
from hvantk.algorithms.rerank.config import FeatureAxis
from hvantk.algorithms.rerank.features import FeatureAssembler

class _Cfg:
    def __init__(self):
        self.features=[
          FeatureAxis("a", lambda: pd.DataFrame({"gene":["A","B"],"x":[1.0,2.0]})),
          FeatureAxis("b", lambda: pd.DataFrame({"gene":["A"],"y":[9.0]})),
        ]
        self.cohort=None

def test_assemble_outer_join_and_coverage():
    m, cov = FeatureAssembler().assemble(_Cfg())
    assert set(m.gene) == {"A","B"}
    assert list(m.columns) == ["gene","x","y"]
    assert m.loc[m.gene=="B","y"].isna().all()         # B missing from axis b -> NaN
    assert cov["b"] == 0.5                              # 1 of 2 genes covered
