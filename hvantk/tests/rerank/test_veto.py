# local/rerank_engine/tests/test_veto.py
import pandas as pd
from hvantk.algorithms.rerank.veto import CaseControlArchitectureVeto, NoOpVeto

def test_noop_vetoes_nothing():
    t = pd.DataFrame({"gene":["A","B"]})
    assert NoOpVeto().apply(t).tolist() == [False, False]

def test_architecture_veto_rules():
    t = pd.DataFrame({
        "gene":     ["few","recurrent","qc","clean"],
        "n_case_var":[2,      9,          24,   8],
        "conc":      [0.3,    0.9,        0.5,  0.2],
        "driver_af": [0.0,    1e-4,       6e-3, 1e-6],
    })
    v = CaseControlArchitectureVeto().apply(t).tolist()
    # few(n<=2)=T ; recurrent(conc>=.6 & af>5e-5)=T ; qc(af>1e-3)=T ; clean=F
    assert v == [True, True, True, False]
