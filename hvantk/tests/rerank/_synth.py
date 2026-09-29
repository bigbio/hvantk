"""Three tiny synthetic datasets, built offline in milliseconds.

Deliberately not fixtures in conftest.py: Tasks 2-7 each run their own test file in
isolation on the login node, and an importable module keeps that selection honest while
still guaranteeing all of them see the same three datasets.

  planted_signal  -- one axis genuinely predicts the label; the control that must CLEAR.
  permuted_labels -- the same matrix with the label shuffled; nothing predicts anything,
                     so every p-value computed on it must be approximately uniform.
  paralogue_groups-- genes in families, families of a controllable size; the fixture the
                     blocked-CV tests and the dominant-block abort run on.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def planted_signal(n=200, n_noise=4, seed=7):
    """`base` carries a real but modest signal; `axis0` carries a planted one; the rest are
    pure noise. Returned as (matrix, y, baseline_cols, axis_groups)."""
    rng = np.random.default_rng(seed)
    y = (rng.random(n) < 0.3).astype(int)
    cols = {"gene": [f"g{i}" for i in range(n)]}
    cols["base"] = y * 0.6 + rng.normal(0, 1.0, n)
    cols["axis0_a"] = y * 1.2 + rng.normal(0, 1.0, n)
    cols["axis0_b"] = y * 1.0 + rng.normal(0, 1.0, n)
    for k in range(n_noise):
        cols[f"axis{k + 1}_a"] = rng.normal(0, 1, n)
        cols[f"axis{k + 1}_b"] = rng.normal(0, 1, n)
    matrix = pd.DataFrame(cols)
    axes = {f"axis{k}": [f"axis{k}_a", f"axis{k}_b"] for k in range(n_noise + 1)}
    return matrix, y, ["base"], axes


def permuted_labels(n=200, n_noise=4, seed=7, shuffle_seed=99):
    """`planted_signal` with the label shuffled: the matrix is identical, so any difference
    in the resulting nulls is attributable to the label and to nothing else."""
    matrix, y, baseline, axes = planted_signal(n=n, n_noise=n_noise, seed=seed)
    return matrix, np.random.default_rng(shuffle_seed).permutation(y), baseline, axes


def paralogue_groups(n=180, family_size=6, dominant=0, seed=3):
    """Genes assigned to HGNC-style families.

    `family_size` genes share each `gene_group`; `dominant` genes (if any) are forced into
    one oversized family, which is how the dominant-block abort is exercised. Returns
    (genes, gene_group mapping) where the mapping is pipe-separated exactly as HGNC's
    `gene_group` column is, with a secondary group on every member so that a naive
    connected-components blocker would collapse the whole matrix.
    """
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n)]
    mapping = {}
    for i, g in enumerate(genes):
        if i < dominant:
            mapping[g] = "Huge family|Shared secondary"
        elif rng.random() < 0.15:
            mapping[g] = ""  # ungrouped: a gene with no family has no paralogue to leak
        else:
            mapping[g] = f"Family {(i - dominant) // family_size}|Shared secondary"
    return genes, mapping


def cheap_scorer(folds=3, seed=0):
    """A fast, honest out-of-fold scorer for the STATISTICAL tests.

    The distributional claims (`p` uniform under a permuted label, a wider candidate set
    raising the null median) need tens of permutations to mean anything, and tens of
    permutations of the real HistGBM is a multi-minute job -- not something that belongs on
    a login node. `permutation_deltas` takes an injectable `scorer` for exactly this reason
    -- swapping in a cheap model keeps the statistics honest without paying for the shipped
    estimator -- so the statistics are exercised here with logistic regression, and the
    DEFAULT scorer -- the real GBM -- is pinned separately by test_default_scorer_is_raw_oof.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    def _score(matrix, cols, y):
        X = np.nan_to_num(matrix[list(cols)].to_numpy(dtype=float))
        cv = StratifiedKFold(folds, shuffle=True, random_state=seed)
        return cross_val_predict(
            LogisticRegression(max_iter=200), X, np.asarray(y), cv=cv,
            method="predict_proba",
        )[:, 1]

    return _score
