# hvantk/tests/rerank/test_config_cohort.py
"""Config/engine consumption of CohortManifest -- the legacy per-rerank cohort loader's
replacement.

The legacy loader is deleted (M1). Config.cohort is now a CohortManifest, and Config
derives its prior from the manifest when the caller doesn't supply one directly
(M2). rerank always needs a prior, and a cohort manifest is now the only supported
source of one, so validate()/rerank() require Config.cohort unconditionally (M3).
The engine's cohort merge fails loud on a column collision instead of silently
dropping the cohort's column (M5).
"""
import logging

import numpy as np
import pandas as pd
import pytest

from hvantk.algorithms.cohort.spec import CohortAxis, CohortManifest, CohortPrior
from hvantk.algorithms.rerank import rerank
from hvantk.algorithms.rerank.audit import CaseControlArchitectureAudit
from hvantk.algorithms.rerank.config import (
    Config,
    FeatureAxis,
    LabelSpec,
    PriorSpec,
    validate,
)


def _write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(str(v) for v in row) for row in rows)
    path.write_text("\n".join(lines) + "\n")


def _manifest(table, prior_col="minp", axes=()):
    return CohortManifest(
        name="demo",
        key="gene",
        table=str(table),
        prior=CohortPrior(column=prior_col, direction="lower_is_better"),
        cohort_axes=axes,
    )


def _mk(cohort=None, prior=None, audit=None, **kw):
    return Config(
        name="t",
        features=[
            FeatureAxis(
                "constraint",
                lambda: pd.DataFrame({"gene": ["A", "B"], "z": [1.0, 2.0]}),
            )
        ],
        labels=LabelSpec(lambda: {"A"}),
        cohort=cohort,
        prior=prior,
        audit=audit,
        **kw,
    )


# ---------------------------------------------------------------------------
# Config.cohort / Config.prior derivation (M1, M2)
# ---------------------------------------------------------------------------


def test_cohort_defaults_to_none():
    assert _mk().cohort is None
    assert _mk().prior is None


def test_config_accepts_a_cohort_manifest(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.2)])
    m = _manifest(p)
    cfg = _mk(cohort=m)
    assert cfg.cohort is m


def test_config_derives_prior_from_cohort_when_none_given(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.2)])
    cfg = _mk(cohort=_manifest(p))

    assert cfg.prior is not None
    loaded = cfg.prior.load()
    assert list(loaded.columns) == ["unit", "prior_stat"]
    assert dict(zip(loaded["unit"], loaded["prior_stat"])) == {"A": 0.01, "B": 0.2}


def test_config_keeps_an_explicitly_given_prior_even_with_a_cohort(tmp_path, caplog):
    # M2 derives the prior only when the caller didn't already supply one -- an
    # explicit PriorSpec is never silently clobbered by the cohort's own prior
    # column (this is exactly registry.build_config's existing calling convention:
    # DiseaseProfile.prior is passed straight through alongside DiseaseProfile.cohort).
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.2)])
    explicit = PriorSpec(path=str(p), unit_col="gene", stat_col="minp")

    with caplog.at_level(logging.WARNING, logger="hvantk.algorithms.rerank.config"):
        cfg = _mk(cohort=_manifest(p), prior=explicit)

    assert cfg.prior is explicit
    # Same (path, key/unit_col, stat/prior column) on both sides -> the sources agree,
    # so no disagreement warning should fire (a warning here would be noise: this is
    # registry.build_config's normal, correct calling convention).
    assert caplog.records == []


def test_config_warns_when_an_explicit_prior_disagrees_with_the_cohorts_own_prior(
    tmp_path, caplog
):
    # The manifest's own prior source and the explicitly-supplied PriorSpec point at
    # genuinely different data (different file, different values), unlike the test
    # above. This lets the assertions below tell "correctly used the explicit prior"
    # apart from "silently used a stale/wrong value" -- a test built on identical
    # values on both sides cannot distinguish those two outcomes.
    cohort_path = tmp_path / "cohort.tsv"
    _write_tsv(cohort_path, ["gene", "minp"], [("A", 0.01), ("B", 0.2)])
    explicit_path = tmp_path / "explicit_prior.tsv"
    _write_tsv(explicit_path, ["gene", "minp"], [("A", 0.9), ("B", 0.4)])
    explicit = PriorSpec(path=str(explicit_path), unit_col="gene", stat_col="minp")

    with caplog.at_level(logging.WARNING, logger="hvantk.algorithms.rerank.config"):
        cfg = _mk(cohort=_manifest(cohort_path), prior=explicit)

    # The warning fires and names both sources.
    assert len(caplog.records) == 1
    message = caplog.records[0].getMessage()
    assert "demo" in message  # names the cohort (Config.name and cohort.name both
    assert str(explicit_path) in message  # the explicit source
    assert str(cohort_path) in message  # the cohort's own source
    assert "disagree" in message

    # The explicit prior still wins -- it must be the *values from explicit_path*
    # that actually reach the output, not the cohort's own (different) values. If the
    # override direction were ever flipped (cohort wins on disagreement instead of
    # the explicit prior), this assertion -- not just the warning above -- would fail.
    loaded = cfg.prior.load()
    assert dict(zip(loaded["unit"], loaded["prior_stat"])) == {"A": 0.9, "B": 0.4}


def test_config_derived_prior_matches_load_prior_frame(tmp_path):
    from hvantk.algorithms.cohort.frame import load_prior_frame

    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.2), ("C", 0.9)])
    manifest = _manifest(p)

    cfg = _mk(cohort=manifest)

    pd.testing.assert_frame_equal(cfg.prior.load(), load_prior_frame(manifest))


# ---------------------------------------------------------------------------
# validate() requires a cohort manifest (M3)
# ---------------------------------------------------------------------------


def test_validate_requires_a_cohort_manifest():
    cfg = _mk(prior=PriorSpec(path="x.tsv", unit_col="gene", stat_col="minp"))
    with pytest.raises(ValueError, match="cohort"):
        validate(cfg)


def test_validate_passes_with_a_cohort_manifest(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.2)])
    validate(_mk(cohort=_manifest(p)))  # no raise


# ---------------------------------------------------------------------------
# Engine cohort merge: loud failure on collision (M5)
# ---------------------------------------------------------------------------


def _synthetic(n=60, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n)]
    y = (rng.random(n) < 0.3).astype(int)
    feat = pd.DataFrame({"gene": genes, "x": y + rng.normal(0, 0.5, n)})
    pos = {g for g, yy in zip(genes, y) if yy == 1}
    return genes, feat, pos


def test_engine_raises_loud_on_cohort_feature_collision(tmp_path):
    genes, feat, pos = _synthetic()
    # The cohort declares a "x" column -- the same name as the "x" feature axis.
    cohort_path = tmp_path / "cohort.tsv"
    _write_tsv(
        cohort_path,
        ["gene", "minp", "x"],
        [(g, 0.1, 99.0) for g in genes],
    )
    cohort = _manifest(cohort_path, axes=(CohortAxis(axis="clash", columns=("x",)),))
    cfg = Config(
        name="t",
        features=[FeatureAxis("x", lambda: feat)],
        labels=LabelSpec(lambda: pos),
        cohort=cohort,
        min_label_coverage=0.0,
    )
    with pytest.raises(ValueError) as excinfo:
        rerank(cfg)
    msg = str(excinfo.value)
    assert "'x'" in msg
    assert "demo" in msg  # names the cohort
    assert "feature axis" in msg  # names the other source


def test_engine_raises_loud_on_cohort_prior_collision(tmp_path):
    genes, feat, pos = _synthetic()
    # The cohort declares a "prior_stat" column -- the exact name engine.py assigns
    # to the merged-in prior, before the cohort merge ever runs.
    cohort_path = tmp_path / "cohort.tsv"
    _write_tsv(
        cohort_path,
        ["gene", "minp", "prior_stat"],
        [(g, 0.1, 7.0) for g in genes],
    )
    cohort = _manifest(
        cohort_path, axes=(CohortAxis(axis="clash", columns=("prior_stat",)),)
    )
    cfg = Config(
        name="t",
        features=[FeatureAxis("x", lambda: feat)],
        labels=LabelSpec(lambda: pos),
        cohort=cohort,
        min_label_coverage=0.0,
    )
    with pytest.raises(ValueError) as excinfo:
        rerank(cfg)
    msg = str(excinfo.value)
    assert "'prior_stat'" in msg
    assert "the prior" in msg


def test_engine_does_not_collide_when_prior_column_reused_as_a_feature(tmp_path):
    # Findings 1+2 (whole-branch review): a manifest honestly declaring its prior
    # column under the same name a feature axis also carries as a model feature (the
    # CHD shape: prior.column='minp', and a 'burden' FeatureAxis whose own column is
    # also 'minp') must run clean. The cohort's prior column was already consumed into
    # 'prior_stat' before the audit merge runs, so re-merging it under its raw name
    # must never happen -- it is not a genuine collision, just the same statistic
    # reachable under two names.
    genes, feat, pos = _synthetic()
    cohort_path = tmp_path / "cohort.tsv"
    _write_tsv(cohort_path, ["gene", "minp"], [(g, 0.1) for g in genes])
    cohort = _manifest(cohort_path, prior_col="minp")
    burden = pd.DataFrame({"gene": genes, "minp": [0.2] * len(genes)})
    cfg = Config(
        name="t",
        features=[FeatureAxis("burden", lambda: burden)],
        labels=LabelSpec(lambda: pos),
        cohort=cohort,
        min_label_coverage=0.0,
    )
    res = rerank(cfg)  # must not raise
    assert "prior_stat" in res.table.columns
    assert len(res.table) == len(genes)
    # 'minp' really was used as the 'burden' axis's model feature (not silently
    # dropped by the collision guard).
    assert "burden" in res.metrics.ablation["family"].tolist()


# ---------------------------------------------------------------------------
# Audit contract: cohort columns feed config.audit.apply() but never become
# model features.
# ---------------------------------------------------------------------------


def test_cohort_columns_feed_audit_but_never_become_model_features(tmp_path):
    genes, feat, pos = _synthetic(n=200)
    prior_vals = {g: (i % 100) / 100.0 for i, g in enumerate(genes)}

    minimal_path = tmp_path / "cohort_minimal.tsv"
    _write_tsv(minimal_path, ["gene", "minp"], [(g, prior_vals[g]) for g in genes])
    cfg_minimal = Config(
        name="t",
        features=[FeatureAxis("x", lambda: feat)],
        labels=LabelSpec(lambda: pos),
        cohort=_manifest(minimal_path),
        min_label_coverage=0.0,
    )
    res_minimal = rerank(cfg_minimal)

    rich_path = tmp_path / "cohort_rich.tsv"
    rows = [
        (
            g,
            prior_vals[g],
            2 if i == 0 else 10,  # n_case_var: gene 0 -> insufficient_data
            0.9 if i == 1 else 0.1,  # conc
            1e-2
            if i == 1
            else 0.0,  # driver_af: gene 1 -> common_driver (conc high too)
        )
        for i, g in enumerate(genes)
    ]
    _write_tsv(rich_path, ["gene", "minp", "n_case_var", "conc", "driver_af"], rows)
    cfg_rich = Config(
        name="t",
        features=[FeatureAxis("x", lambda: feat)],
        labels=LabelSpec(lambda: pos),
        cohort=_manifest(
            rich_path,
            axes=(
                CohortAxis(
                    axis="architecture", columns=("n_case_var", "conc", "driver_af")
                ),
            ),
        ),
        audit=CaseControlArchitectureAudit(),
        min_label_coverage=0.0,
    )
    res_rich = rerank(cfg_rich)

    # Same feature matrix, same labels, same prior values in both runs -> scores and
    # tiers (which depend only on the fitted model) must be bit-identical. If a
    # cohort-only column ("n_case_var" et al.) had leaked into the feature matrix, the
    # fitted GBM -- and therefore these columns -- would differ between the two runs.
    pd.testing.assert_series_equal(res_minimal.table["score"], res_rich.table["score"])
    pd.testing.assert_series_equal(res_minimal.table["tier"], res_rich.table["tier"])

    # The cohort's audit-only columns DID feed the audit.
    assert not res_minimal.table["flag"].any()
    assert res_rich.table["flag"].any()

    # And no cohort-only column ever became (or looks like) a model feature: the
    # per-axis ablation table is keyed strictly by config.features' axis names.
    assert set(res_rich.metrics.ablation["family"]) == {"x"}
    assert "n_case_var" not in set(res_rich.metrics.ablation["family"])
