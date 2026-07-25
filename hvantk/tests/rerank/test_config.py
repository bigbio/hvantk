# local/rerank_engine/tests/test_config.py
import pytest, pandas as pd
from hvantk.algorithms.cohort.spec import CohortManifest, CohortPrior
from hvantk.algorithms.rerank.config import (
    Config,
    PriorSpec,
    FeatureAxis,
    LabelSpec,
    validate,
)
from hvantk.algorithms.rerank.audit import NoAudit


def _mk(cohort=None):
    return Config(
        name="t",
        prior=PriorSpec(path="x.tsv", unit_col="gene", stat_col="minp"),
        features=[
            FeatureAxis(
                "constraint", lambda: pd.DataFrame({"gene": ["A"], "mis_z": [1.0]})
            )
        ],
        labels=LabelSpec(lambda: {"A"}),
        cohort=cohort,
    )


def test_defaults():
    c = _mk()
    assert (
        c.units == "gene"
        and c.cohort is None
        and isinstance(c.audit, NoAudit)
        and c.tiers == 5
    )
    assert c.calibration == "isotonic" and c.folds == 5


def test_validate_rejects_variant_units():
    c = _mk()
    c.units = "variant"
    with pytest.raises(NotImplementedError):
        validate(c)


def test_validate_requires_a_cohort_manifest():
    # M3: rerank always needs a prior, and a CohortManifest is now its only supported
    # source -- a PriorSpec-only config (cohort=None) is rejected even with NoAudit.
    with pytest.raises(ValueError, match="cohort"):
        validate(_mk())


def test_validate_passes_minimal():
    manifest = CohortManifest(
        name="t",
        key="gene",
        table="x.tsv",
        prior=CohortPrior(column="minp", direction="lower_is_better"),
    )
    validate(_mk(cohort=manifest))  # no raise


def test_validate_rejects_empty_labels():
    import pytest, pandas as pd
    from hvantk.algorithms.rerank.config import (
        Config,
        PriorSpec,
        FeatureAxis,
        LabelSpec,
        validate,
    )

    manifest = CohortManifest(
        name="t",
        key="gene",
        table="x.tsv",
        prior=CohortPrior(column="minp", direction="lower_is_better"),
    )
    c = Config(
        name="t",
        prior=PriorSpec("x.tsv", "gene", "minp"),
        features=[FeatureAxis("a", lambda: pd.DataFrame({"gene": ["A"], "x": [1.0]}))],
        labels=LabelSpec(lambda: set()),
        cohort=manifest,
    )
    with pytest.raises(ValueError):
        validate(c)
