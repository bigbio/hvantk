# hvantk/tests/rerank/test_audit.py
import pandas as pd
import pytest
from hvantk.algorithms.rerank.audit import (
    ARCHITECTURE_AUDIT_COLUMNS,
    CaseControlArchitectureAudit,
    NoAudit,
)


def test_noaudit_flags_nothing():
    t = pd.DataFrame({"gene": ["A", "B"]})
    assert NoAudit().apply(t).tolist() == ["", ""]


def test_architecture_audit_reasons():
    t = pd.DataFrame(
        {
            "gene": ["few", "recurrent", "qc", "clean"],
            "n_case_var": [2, 9, 24, 8],
            "conc": [0.3, 0.9, 0.5, 0.2],
            "driver_af": [0.0, 1e-4, 6e-3, 1e-6],
        }
    )
    r = CaseControlArchitectureAudit().apply(t).tolist()
    # few(n<=2)=insufficient_data ; recurrent(conc>=.6 & af>5e-5)=recurrent_variant ;
    # qc(af>1e-3)=common_driver ; clean="" ; flagged set == old veto union
    assert r == ["insufficient_data", "recurrent_variant", "common_driver", ""]
    flagged = [x != "" for x in r]
    assert flagged == [True, True, True, False]


def test_architecture_audit_missing_columns_error_names_every_required_column():
    # Finding 5 (re-review): apply() must derive its own requirement check from
    # ARCHITECTURE_AUDIT_COLUMNS (the constant it and has_architecture_columns() are
    # both supposed to agree on) rather than re-listing ("n_case_var", "conc")
    # separately -- otherwise the two can drift apart silently if the constant is
    # ever extended.
    t = pd.DataFrame({"gene": ["A", "B"]})
    with pytest.raises(ValueError) as excinfo:
        CaseControlArchitectureAudit().apply(t)
    msg = str(excinfo.value)
    for col in ARCHITECTURE_AUDIT_COLUMNS:
        assert col in msg
    # Finding 3 (re-review): the remedy text must not claim there's a way to declare
    # a required column "directly" outside cohort_axes (there isn't -- the manifest
    # schema allows only 'prior' and 'cohort_axes'), and must not blindly point CLI
    # users at leaving Config.audit unset (rerank_cli.py always sets it explicitly).
    assert "directly or via a cohort_axes entry" not in msg
    assert "hvantk rerank" in msg


def test_architecture_audit_precedence_on_overlap():
    # A gene matching multiple rules gets the highest-precedence reason:
    # common_driver > recurrent_variant > insufficient_data.
    t = pd.DataFrame(
        {
            "gene": ["few_and_common", "recurrent_and_common", "few_and_recurrent"],
            "n_case_var": [2, 9, 2],
            "conc": [0.9, 0.9, 0.9],
            "driver_af": [6e-3, 6e-3, 1e-4],
        }
    )
    r = CaseControlArchitectureAudit().apply(t).tolist()
    assert r == ["common_driver", "common_driver", "recurrent_variant"]
