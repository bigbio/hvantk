# local/rerank_engine/tests/test_config.py
import pytest, pandas as pd
from hvantk.algorithms.rerank.config import Config, PriorSpec, FeatureAxis, LabelSpec, validate
from hvantk.algorithms.rerank.veto import NoOpVeto

def _mk():
    return Config(name="t",
        prior=PriorSpec(path="x.tsv", unit_col="gene", stat_col="minp"),
        features=[FeatureAxis("constraint", lambda: pd.DataFrame({"gene":["A"],"mis_z":[1.0]}))],
        labels=LabelSpec(lambda: {"A"}))

def test_defaults():
    c = _mk()
    assert c.units == "gene" and c.cohort is None and isinstance(c.veto, NoOpVeto) and c.tiers == 5
    assert c.calibration == "isotonic" and c.folds == 5

def test_validate_rejects_variant_units():
    c = _mk(); c.units = "variant"
    with pytest.raises(NotImplementedError):
        validate(c)

def test_validate_requires_noop_veto_when_no_cohort():
    from hvantk.algorithms.rerank.veto import CaseControlArchitectureVeto
    c = _mk(); c.veto = CaseControlArchitectureVeto()   # cohort is None -> invalid
    with pytest.raises(ValueError):
        validate(c)

def test_validate_passes_minimal():
    validate(_mk())   # no raise

def test_validate_rejects_empty_labels():
    import pytest, pandas as pd
    from hvantk.algorithms.rerank.config import Config, PriorSpec, FeatureAxis, LabelSpec, validate
    c = Config(name="t", prior=PriorSpec("x.tsv","gene","minp"),
               features=[FeatureAxis("a", lambda: pd.DataFrame({"gene":["A"],"x":[1.0]}))],
               labels=LabelSpec(lambda: set()))
    with pytest.raises(ValueError):
        validate(c)
