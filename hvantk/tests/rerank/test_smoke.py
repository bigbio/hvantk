# hvantk/tests/rerank/test_smoke.py
import dataclasses

import numpy as np, pandas as pd
from hvantk.algorithms.cohort.spec import CohortManifest, CohortPrior
from hvantk.algorithms.rerank import (
    rerank, Config, PriorSpec, FeatureAxis, LabelSpec, NullConfig,
)
from hvantk.algorithms.rerank.evaluator import ABLATION_FOLDS


def test_engine_end_to_end(tmp_path):
    rng = np.random.default_rng(0)
    n = 200
    genes = [f"g{i}" for i in range(n)]
    y = (rng.random(n) < 0.3).astype(int)
    feat = pd.DataFrame({"gene": genes, "x": y + rng.normal(0, 0.5, n)})
    prior = pd.DataFrame({"gene": genes, "p": rng.random(n)})
    pf = tmp_path / "p.tsv"
    prior.to_csv(pf, sep="\t", index=False)
    pos = {g for g, yy in zip(genes, y) if yy == 1}
    # rerank() requires a cohort manifest (M3); reuse pf as the cohort table too, so
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

    # Guard the engine's null wiring: Tasks 4 and 6 both edit `_run_nulls`, and a
    # regression that computed the observed deltas with a different scorer from the one
    # `permutation_deltas` used would silently produce wrong p-values with nothing here to
    # catch it. `z` is pure noise, offered alongside the existing `x` axis (which stays the
    # baseline, offered first), so the null must cover only `z`, and the null's own observed
    # delta for `z` must agree with the ablation table's independently-computed AUC gap.
    noise = pd.DataFrame({"gene": genes, "z": rng.normal(0, 1, n)})
    cfg2 = dataclasses.replace(
        cfg,
        features=[cfg.features[0], FeatureAxis("z", lambda: noise)],
        nulls=NullConfig(n_perm=2, seed=1),
    )
    res2 = rerank(cfg2)
    assert res2.nulls is not None
    assert res2.nulls.axes == ("z",)
    assert res2.nulls.setting.folds == ABLATION_FOLDS
    abl = res2.metrics.ablation.set_index("family")["auc"]
    assert abs((abl["z"] - abl["x"]) - res2.nulls.observed["z"]) <= 1e-3
