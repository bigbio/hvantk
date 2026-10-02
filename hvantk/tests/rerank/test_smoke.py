# hvantk/tests/rerank/test_smoke.py
import dataclasses

import numpy as np, pandas as pd
from hvantk.algorithms.cohort.spec import CohortManifest, CohortPrior
from hvantk.algorithms.rerank import (
    rerank,
    Config,
    PriorSpec,
    FeatureAxis,
    LabelSpec,
    NullConfig,
    BlockPolicy,
)
from hvantk.algorithms.rerank.evaluator import ABLATION_FOLDS


def test_engine_end_to_end(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    n = 200
    genes = [f"g{i}" for i in range(n)]
    y = (rng.random(n) < 0.3).astype(int)
    feat = pd.DataFrame({"gene": genes, "x": y + rng.normal(0, 0.5, n)})
    prior = pd.DataFrame({"gene": genes, "p": rng.random(n)})
    pf = tmp_path / "p.tsv"
    prior.to_csv(pf, sep="\t", index=False)
    pos = {g for g, yy in zip(genes, y) if yy == 1}
    # rerank() requires a cohort manifest; reuse pf as the cohort table too, so
    # Config.prior (given explicitly here, matching the old shape) and the cohort's
    # own prior column agree on the exact same values.
    cohort = CohortManifest(
        name="t",
        key="gene",
        table=str(pf),
        prior=CohortPrior(column="p", direction="lower_is_better"),
    )
    cfg = Config(
        name="t",
        prior=PriorSpec(str(pf), "gene", "p"),
        cohort=cohort,
        features=[FeatureAxis("x", lambda: feat)],
        labels=LabelSpec(lambda: pos),
        min_label_coverage=0.0,
    )
    res = rerank(cfg)
    assert len(res.table) == n and res.metrics.auc > 0.7
    assert {"gene", "score", "tier", "verdict"} <= set(res.table.columns)

    # Guard the engine's null wiring: `_run_nulls` computing the observed deltas with a
    # different scorer from the one `permutation_deltas` used would silently produce wrong
    # p-values with nothing here to catch it. `z` is pure noise, offered alongside the
    # existing `x` axis (which stays the baseline, offered first), so the null must cover
    # only `z`, and the null's own observed delta for `z` must agree with the ablation
    # table's independently-computed AUC gap.
    noise = pd.DataFrame({"gene": genes, "z": rng.normal(0, 1, n)})
    # Engine-wiring guard: a real block table (families of 4, status Approved) so
    # `_run_nulls` is exercised on an actual blocking, not just blocks=None.
    groups_path = tmp_path / "hgnc_groups.tsv"
    pd.DataFrame(
        {
            "symbol": genes,
            "gene_group": [f"Family {i // 4}" for i in range(n)],
            "status": ["Approved"] * n,
        }
    ).to_csv(groups_path, sep="\t", index=False)
    cfg2 = dataclasses.replace(
        cfg,
        features=[cfg.features[0], FeatureAxis("z", lambda: noise)],
        nulls=NullConfig(n_perm=2, seed=1),
        folds=3,
        blocks=BlockPolicy(table=str(groups_path)),
        # A non-default Config.seed: a _run_nulls that dropped config.seed on the floor
        # (e.g. reverted to the historic DEFAULT_SEED) would then compute its own `observed`
        # delta under a DIFFERENT CV partition from the one the ablation table below used,
        # and the two independently-computed deltas would no longer agree to 1e-3.
        seed=7,
    )

    # Guard the headline block wiring: spy on ReRanker.score at the class level so every
    # call this run makes -- however it is reached -- is checked against the SAME block
    # array the engine resolved, not merely that the run as a whole "looks blocked".
    from hvantk.algorithms.rerank.reranker import ReRanker

    calls = []
    real_score = ReRanker.score

    def _spy_score(self, *args, **kwargs):
        calls.append(kwargs.get("groups"))
        return real_score(self, *args, **kwargs)

    monkeypatch.setattr(ReRanker, "score", _spy_score)
    res2 = rerank(cfg2)
    assert calls, "ReRanker.score was never called for the blocked run"
    assert all(np.array_equal(g, res2.blocks.blocks) for g in calls)
    assert res2.nulls is not None
    assert res2.nulls.axes == ("z",)
    assert res2.nulls.setting.folds == ABLATION_FOLDS
    assert res2.nulls.setting.seed == 7  # the CV seed the null's scorer ran under
    assert res2.nulls.setting.block_digest is not None
    assert res2.blocks.digest == res2.nulls.setting.block_digest
    abl = res2.metrics.ablation.set_index("family")["auc"]
    assert abs((abl["z"] - abl["x"]) - res2.nulls.observed["z"]) <= 1e-3
