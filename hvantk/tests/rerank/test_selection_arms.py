"""Two provenance arms, the global pass, and the selection summary.

The `clean` arm is the headline; `all` exists only so the circularity channel is a
measured number. Both arms run over identical folds, so their difference is paired.
"""
import numpy as np
import pandas as pd


def _cfg(tmp_path, **kw):
    """Minimal Config over a tiny in-memory matrix with one real signal column.

    A CohortManifest is mandatory for `rerank` (the cohort contract), so the prior is
    written out and declared here even though nothing in these tests reads it -- the
    subject under test is column selection, not the prior.
    """
    from hvantk.algorithms.cohort.spec import CohortManifest, CohortPrior
    from hvantk.algorithms.rerank.config import Config, FeatureAxis, LabelSpec

    rng = np.random.default_rng(11)
    n = 800
    genes = [f"G{i}" for i in range(n)]
    y = np.repeat([0, 1], n // 2)
    base = pd.DataFrame({"gene": genes, "c1": y + rng.normal(0, 1.0, n),
                         "c2": rng.normal(0, 1, n)})
    extra = pd.DataFrame({"gene": genes, "REVEL_rankscore": y * 1.0})   # perfect + circular
    pos = {g for g, v in zip(genes, y) if v == 1}

    prior_path = tmp_path / "prior.tsv"
    pd.DataFrame({"gene": genes, "p": rng.random(n)}).to_csv(
        prior_path, sep="\t", index=False
    )
    cohort = CohortManifest(
        name="t",
        key="gene",
        table=str(prior_path),
        prior=CohortPrior(column="p", direction="lower_is_better"),
    )
    return Config(
        name="t",
        cohort=cohort,
        features=[FeatureAxis("constraint", lambda: base),
                  FeatureAxis("trained", lambda: extra)],
        labels=LabelSpec(loader=lambda: pos),
        **kw,
    )


def test_arms_split_on_provenance_and_delta_is_reported(tmp_path):
    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(
        tmp_path,
        selection=SelectionPolicy(wrapper="none"),
        feature_provenance={"c1": frozenset(), "c2": frozenset(),
                            "REVEL_rankscore": frozenset({"ClinVar"})},
        label_provenance=frozenset({"GenCC"}),
        min_label_coverage=0.0,
    )

    arms = rerank_arms(cfg)

    assert set(arms) == {"clean", "all"}
    assert arms["all"].metrics.auc > arms["clean"].metrics.auc  # circular column helps
    assert arms["clean"].selection.n_conflicted == 1
    assert "REVEL_rankscore" not in arms["clean"].selection.global_features.get("trained", ())


def test_undeclared_column_is_unknown_and_excluded_from_clean(tmp_path):
    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(
        tmp_path,
        selection=SelectionPolicy(wrapper="none"),
        feature_provenance={"c1": frozenset(), "c2": None,
                            "REVEL_rankscore": frozenset({"ClinVar"})},
        label_provenance=frozenset({"GenCC"}),
        min_label_coverage=0.0,
    )

    arms = rerank_arms(cfg)

    assert arms["clean"].selection.n_unknown == 1


def test_selection_summary_carries_frequency_and_global_list(tmp_path):
    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(tmp_path, selection=SelectionPolicy(wrapper="none"),
               feature_provenance=None, label_provenance=frozenset(),
               min_label_coverage=0.0)

    s = rerank_arms(cfg)["all"].selection

    assert s.frequency["constraint"]["c1"] >= 1
    assert "c1" in s.global_features["constraint"]
    assert np.isfinite(s.auc_nested) and np.isfinite(s.auc_global)


def test_no_selection_policy_returns_a_single_arm_unchanged(tmp_path):
    from hvantk.algorithms.rerank.engine import rerank, rerank_arms

    cfg = _cfg(tmp_path, min_label_coverage=0.0)

    arms = rerank_arms(cfg)

    assert list(arms) == ["all"]
    assert arms["all"].selection is None
    assert rerank(cfg).metrics.auc == arms["all"].metrics.auc


def test_noise_column_is_dropped_from_the_global_list(tmp_path):
    """The readable feature list must reflect selection, not just echo the inputs."""
    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(tmp_path, selection=SelectionPolicy(wrapper="none"),
               feature_provenance=None, label_provenance=frozenset(),
               min_label_coverage=0.0)

    s = rerank_arms(cfg)["all"].selection

    assert "c2" not in s.global_features["constraint"]


def test_clean_arm_refuses_to_run_when_every_column_conflicts(tmp_path):
    """An empty clean arm is a configuration error, not a silently degraded run.

    Scoring nothing would otherwise surface as an opaque failure deep inside the
    assembler, or worse, as a plausible-looking chance-level AUC.
    """
    import pytest

    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(
        tmp_path,
        selection=SelectionPolicy(wrapper="none"),
        feature_provenance={c: frozenset({"ClinVar"})
                            for c in ("c1", "c2", "REVEL_rankscore")},
        label_provenance=frozenset({"GenCC"}),
        min_label_coverage=0.0,
    )

    with pytest.raises(ValueError, match="clean"):
        rerank_arms(cfg)


def test_column_missing_from_the_provenance_map_is_undeclared_not_deleted(tmp_path):
    """An axis nobody declared must land in `all` as unknown -- never silently vanish.

    Omission and an explicit None mean the same thing per the plugin contract: usable,
    never clean. If an undeclared column were dropped from both arms instead, adding an
    axis and forgetting to declare it would quietly change the headline with nothing in
    the output to show a column had gone missing.
    """
    from hvantk.algorithms.rerank.engine import rerank_arms
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    cfg = _cfg(
        tmp_path,
        selection=SelectionPolicy(wrapper="none"),
        # c2 and REVEL_rankscore are simply absent from the map
        feature_provenance={"c1": frozenset()},
        label_provenance=frozenset({"GenCC"}),
        min_label_coverage=0.0,
    )

    arms = rerank_arms(cfg)

    assert arms["all"].selection.n_unknown == 2
    assert "REVEL_rankscore" in arms["all"].selection.global_features.get("trained", ())
    assert "trained" not in arms["clean"].selection.global_features
