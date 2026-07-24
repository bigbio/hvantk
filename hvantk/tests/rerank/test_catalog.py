# hvantk/tests/rerank/test_catalog.py
# Trimmed for hvantk CI: domain-specific axis builders (constraint/burden/subgenic/labels/ptm)
# stay in local/rerank_engine/catalog/axes.py (workspace recipes).
# Only registry mechanics + generic builders are tested here.
from hvantk.algorithms.rerank.catalog.profile import DiseaseProfile


def test_profile_defaults():
    p = DiseaseProfile(name="t")
    assert p.cohort is None and p.precomputed is None and p.tissue is None
    assert p.min_label_coverage == 0.5
    assert p.label_classifications == ["Definitive", "Strong", "Moderate"]
    assert p.cell_types == [] and p.disease_terms == [] and p.extra_flagged_genes == []
    assert p.extra_positive_genes == set()


# Task 2: Registry mechanics tests
import pandas as pd
from hvantk.algorithms.rerank.config import FeatureAxis
from hvantk.algorithms.rerank.catalog.profile import DiseaseProfile
from hvantk.algorithms.rerank.catalog.registry import axis, AXES, build_config, _constraint_first, _default_audit
from hvantk.algorithms.rerank.audit import NoAudit, CaseControlArchitectureAudit


def test_constraint_first_reorders():
    feats = [FeatureAxis("expression", lambda: pd.DataFrame({"gene": ["A"], "x": [1.0]})),
             FeatureAxis("constraint", lambda: pd.DataFrame({"gene": ["A"], "z": [1.0]}))]
    out = _constraint_first(feats)
    assert out[0].name == "constraint"


def test_default_audit_depends_on_cohort():
    from hvantk.algorithms.rerank.config import CohortSpec
    assert isinstance(_default_audit(DiseaseProfile(name="t")), NoAudit)
    p = DiseaseProfile(name="t", cohort=CohortSpec("x.tsv", params={}))
    assert isinstance(_default_audit(p), CaseControlArchitectureAudit)


def test_build_config_string_and_override():
    import pandas as pd
    from hvantk.algorithms.rerank.config import FeatureAxis, LabelSpec
    from hvantk.algorithms.rerank.catalog.registry import axis, AXES, build_config
    from hvantk.algorithms.rerank.catalog.profile import DiseaseProfile
    from hvantk.algorithms.rerank.audit import NoAudit
    _saved = dict(AXES)
    try:
        @axis("labels")
        def _labels(profile):
            return LabelSpec(lambda: {"A"})

        @axis("dummy_constraint")
        def _dummy_constraint(profile):
            return FeatureAxis("constraint", lambda: pd.DataFrame({"gene": ["A", "B"], "z": [1.0, 2.0]}))

        from hvantk.algorithms.rerank.config import PriorSpec
        p = DiseaseProfile(name="t", disease_terms=[],
                           prior=PriorSpec(path="/dev/null", unit_col="gene", stat_col="score"))
        override = FeatureAxis("expression", lambda: pd.DataFrame({"gene": ["A", "B"], "x": [3.0, 4.0]}))
        cfg = build_config(p, axes=["dummy_constraint", override])
        assert cfg.name == "t"
        assert [f.name for f in cfg.features][0] == "constraint"
        assert "expression" in [f.name for f in cfg.features]
        assert isinstance(cfg.audit, NoAudit)
    finally:
        AXES.clear()
        AXES.update(_saved)


# Task 3: Generic builders (table_axis / genelist_labels) - path-free synthetic tests
def test_table_axis_parquet(tmp_path):
    import pandas as pd
    from hvantk.algorithms.rerank.catalog.builders import table_axis
    df = pd.DataFrame({"gene": ["A", "B"], "x": [1.0, 2.0]})
    p = tmp_path / "feat.parquet"
    df.to_parquet(p)
    ax = table_axis("myaxis", str(p))
    loaded = ax.load()
    assert list(loaded.columns) == ["gene", "x"]
    assert ax.name == "myaxis"


def test_table_axis_tsv(tmp_path):
    from hvantk.algorithms.rerank.catalog.builders import table_axis
    p = tmp_path / "feat.tsv"
    p.write_text("gene\tx\nA\t1.0\nB\t2.0\n")
    ax = table_axis("tsv_axis", str(p))
    loaded = ax.load()
    assert list(loaded["gene"]) == ["A", "B"]


def test_genelist_labels(tmp_path):
    from hvantk.algorithms.rerank.catalog.builders import genelist_labels
    p = tmp_path / "genes.txt"
    p.write_text("GENE1\nGENE2\n# comment\n\n")
    ls = genelist_labels(str(p))
    s = ls.load()
    assert s == {"GENE1", "GENE2"}


def test_axis_registry_decorator():
    """axis() decorator registers a builder in AXES and it can be looked up."""
    from hvantk.algorithms.rerank.catalog.registry import axis, AXES
    _saved = dict(AXES)
    try:
        @axis("test_synthetic")
        def _build(profile):
            return FeatureAxis("test_synthetic", lambda: pd.DataFrame({"gene": ["X"], "v": [1.0]}))

        assert "test_synthetic" in AXES
        p = DiseaseProfile(name="t")
        ax = AXES["test_synthetic"](p)
        assert ax.name == "test_synthetic"
    finally:
        AXES.clear()
        AXES.update(_saved)
