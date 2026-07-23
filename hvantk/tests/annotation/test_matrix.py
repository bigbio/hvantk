"""matrix.py: EWCE specificity + summary-AnnData -> per-gene table."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_ewce_specificity_is_row_normalized_mean():
    from hvantk.algorithms.annotation.matrix import ewce_specificity

    # gene X: mean 3 in group A, 1 in B -> spec 0.75/0.25; gene Y: 0 everywhere -> 0.
    mean_gg = pd.DataFrame({"A": [3.0, 0.0], "B": [1.0, 0.0]}, index=["X", "Y"])
    spec = ewce_specificity(mean_gg)
    assert spec.loc["X", "A"] == pytest.approx(0.75)
    assert spec.loc["X", "B"] == pytest.approx(0.25)
    assert spec.loc["Y", "A"] == pytest.approx(0.0)  # unexpressed gene -> 0, not NaN


def _summary_adata():
    import anndata as ad

    # 2 groups x 3 genes. group ' CM ' (whitespace) is the specificity target.
    genes = ["TNNT2", "ACTB", "DUP"]
    groups = [" CM ", "Other"]
    mean = np.array([[10.0, 5.0, 4.0], [0.0, 5.0, 4.0]])  # groups x genes
    frac = np.array([[0.9, 0.8, 0.7], [0.1, 0.8, 0.7]])
    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": groups, "n_cells": [100, 200]}, index=groups),
        var=pd.DataFrame(index=genes),
        layers={"mean": mean, "fraction_expressed": frac},
    )
    return a


def _mspec(stats=(), targets=("CM",), name="cm_spec"):
    from hvantk.algorithms.annotation.spec import MatrixSpec, SpecificitySpec

    return MatrixSpec(
        group_axis="celltype",
        atlas="asp",
        tissue_tag="cardiac",
        stats=tuple(stats),
        drop_groups=(),
        specificity=SpecificitySpec(
            method="ewce_fraction", targets=tuple(targets), combine="max", name=name
        ),
    )


def test_reduce_emits_specificity_and_strips_group_whitespace():
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene

    df = reduce_matrix_to_gene(_summary_adata(), _mspec())
    assert list(df.columns) == ["symbol", "asp_cm_spec"]
    row = df.set_index("symbol")
    # TNNT2: mean 10 in CM, 0 in Other -> spec 1.0. ACTB: 5/5 -> 0.5.
    assert row.loc["TNNT2", "asp_cm_spec"] == pytest.approx(1.0)
    assert row.loc["ACTB", "asp_cm_spec"] == pytest.approx(0.5)


def test_reduce_collapses_duplicate_symbols_with_max():
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene

    a = _summary_adata()
    a.var_names = ["TNNT2", "ACTB", "TNNT2"]  # DUP renamed to a duplicate TNNT2
    df = reduce_matrix_to_gene(a, _mspec()).set_index("symbol")
    assert df.index.is_unique  # one row per symbol
    # the two TNNT2 rows (spec 1.0 and 0.5) collapse to max 1.0
    assert df.loc["TNNT2", "asp_cm_spec"] == pytest.approx(1.0)


def test_reduce_sum_combine_pools_subtype_fractions_for_a_cell_class():
    # A cell CLASS (e.g. cardiomyocytes) split across subtypes: combine='sum' pools the
    # per-subtype EWCE fractions so a pan-class gene reads high; 'max' would take only the
    # single strongest subtype and under-read it.
    import anndata as ad

    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene
    from hvantk.algorithms.annotation.spec import MatrixSpec, SpecificitySpec

    genes = ["PAN", "SPECIFIC"]
    groups = ["CM_a", "CM_b", "Other"]
    # PAN: expressed equally in both CM subtypes, nowhere else -> per-subtype spec 0.5, sum 1.0.
    # SPECIFIC: only in CM_a -> per-subtype spec 1.0/0.0, sum 1.0, but max is also 1.0.
    mean = np.array([[5.0, 5.0], [5.0, 0.0], [0.0, 0.0]])  # groups x genes
    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": groups}, index=groups),
        var=pd.DataFrame(index=genes),
        layers={"mean": mean},
    )
    mspec = MatrixSpec(
        group_axis="celltype",
        atlas="asp",
        specificity=SpecificitySpec(
            method="ewce_fraction",
            targets=("CM_a", "CM_b"),
            combine="sum",
            name="cm_spec",
        ),
    )
    row = reduce_matrix_to_gene(a, mspec).set_index("symbol")
    # PAN pooled across both CM subtypes reads fully CM (0.5 + 0.5); 'max' would give only 0.5.
    assert row.loc["PAN", "asp_cm_spec"] == pytest.approx(1.0)
    assert row.loc["SPECIFIC", "asp_cm_spec"] == pytest.approx(1.0)


def test_reduce_emits_per_group_stats_when_declared():
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene

    df = reduce_matrix_to_gene(_summary_adata(), _mspec(stats=("mean",))).set_index(
        "symbol"
    )
    # per-group mean columns, group label sanitized (' CM ' -> 'cm'); plus the specificity col
    assert "asp_cm_mean" in df.columns and "asp_other_mean" in df.columns
    assert df.loc["TNNT2", "asp_cm_mean"] == pytest.approx(10.0)
