# hvantk/tests/rerank/test_audit.py
import pandas as pd
from hvantk.algorithms.rerank.audit import CaseControlArchitectureAudit, NoAudit


def test_noaudit_flags_nothing():
    t = pd.DataFrame({"gene": ["A", "B"]})
    assert NoAudit().apply(t).tolist() == ["", ""]


def test_architecture_audit_reasons():
    t = pd.DataFrame({
        "gene":       ["few", "recurrent", "qc",  "clean"],
        "n_case_var": [2,      9,           24,    8],
        "conc":       [0.3,    0.9,         0.5,   0.2],
        "driver_af":  [0.0,    1e-4,        6e-3,  1e-6],
    })
    r = CaseControlArchitectureAudit().apply(t).tolist()
    # few(n<=2)=insufficient_data ; recurrent(conc>=.6 & af>5e-5)=recurrent_variant ;
    # qc(af>1e-3)=common_driver ; clean="" ; flagged set == old veto union
    assert r == ["insufficient_data", "recurrent_variant", "common_driver", ""]
    flagged = [x != "" for x in r]
    assert flagged == [True, True, True, False]


def test_architecture_audit_precedence_on_overlap():
    # A gene matching multiple rules gets the highest-precedence reason:
    # common_driver > recurrent_variant > insufficient_data.
    t = pd.DataFrame({
        "gene":       ["few_and_common", "recurrent_and_common", "few_and_recurrent"],
        "n_case_var": [2,                 9,                      2],
        "conc":       [0.9,               0.9,                    0.9],
        "driver_af":  [6e-3,              6e-3,                   1e-4],
    })
    r = CaseControlArchitectureAudit().apply(t).tolist()
    assert r == ["common_driver", "common_driver", "recurrent_variant"]
