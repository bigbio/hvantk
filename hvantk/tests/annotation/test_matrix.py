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
    # vector-by-default: naming targets adds the roll-up, it does not replace the vector.
    assert sorted(df.columns) == ["asp_cm", "asp_cm_spec", "asp_other", "symbol"]
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


def test_reduce_raises_when_a_declared_specificity_target_is_absent():
    # A partially-present target set must fail loudly, not silently pool over the survivors
    # (which would under-report every pan-class gene's cell-class specificity).
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene

    # _summary_adata has groups 'CM' and 'Other'; 'GHOST' is not present.
    with pytest.raises(ValueError, match="GHOST"):
        reduce_matrix_to_gene(_summary_adata(), _mspec(targets=("CM", "GHOST")))


def test_reduce_raises_on_sanitized_group_label_collision():
    import anndata as ad

    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene
    from hvantk.algorithms.annotation.spec import MatrixSpec

    # 'T-cell' and 'T cell' both sanitize to 't_cell' -> the per-group stat columns would
    # otherwise silently overwrite each other.
    groups = ["T-cell", "T cell"]
    a = ad.AnnData(
        X=None,
        obs=pd.DataFrame({"celltype": groups}, index=groups),
        var=pd.DataFrame(index=["A", "B"]),
        layers={"mean": np.array([[1.0, 2.0], [3.0, 4.0]])},
    )
    mspec = MatrixSpec(
        group_axis="celltype", atlas="asp", stats=("mean",), specificity=None
    )
    with pytest.raises(ValueError, match="sanitize"):
        reduce_matrix_to_gene(a, mspec)


def test_reduce_emits_per_group_stats_when_declared():
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene

    df = reduce_matrix_to_gene(_summary_adata(), _mspec(stats=("mean",))).set_index(
        "symbol"
    )
    # per-group mean columns, group label sanitized (' CM ' -> 'cm'); plus the specificity col
    assert "asp_cm_mean" in df.columns and "asp_other_mean" in df.columns
    assert df.loc["TNNT2", "asp_cm_mean"] == pytest.approx(10.0)


# --- vector emission -------------------------------------------------------------------
# Reducing an atlas to one summed scalar discards the cross-cell-type contrast entirely: the
# non-target groups are computed, used as the denominator, then thrown away. Measured cost on
# real data (analysis/rerank-homogenised): collapsing 33 cortical cell types to one number
# took an epilepsy axis from +0.0795 to +0.0187 and from significant to not, while the
# multiplicity penalty for keeping the vector was +0.0009. So the vector is the default and
# the roll-up is additive.


def test_specificity_emits_one_column_per_group_by_default():
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene
    from hvantk.algorithms.annotation.spec import MatrixSpec, SpecificitySpec

    mspec = MatrixSpec(
        group_axis="celltype",
        atlas="asp",
        specificity=SpecificitySpec(method="ewce_fraction", targets=()),
    )
    df = reduce_matrix_to_gene(_summary_adata(), mspec).set_index("symbol")

    # one column per surviving group, sanitized, atlas-prefixed -- and NO roll-up, because
    # no targets were named.
    assert sorted(df.columns) == ["asp_cm", "asp_other"]
    # TNNT2 is 10 in CM, 0 in Other -> the vector carries the contrast the scalar hid.
    assert df.loc["TNNT2", "asp_cm"] == pytest.approx(1.0)
    assert df.loc["TNNT2", "asp_other"] == pytest.approx(0.0)
    assert df.loc["ACTB", "asp_cm"] == pytest.approx(0.5)
    assert df.loc["ACTB", "asp_other"] == pytest.approx(0.5)


def test_emit_rollup_suppresses_the_vector():
    """The old single-column behaviour stays reachable for callers that want only the class.

    Vector-by-default changes what you get without asking; it must not remove the ability to
    ask for the scalar. `emit="rollup"` is that opt-out.
    """
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene
    from hvantk.algorithms.annotation.spec import MatrixSpec, SpecificitySpec

    mspec = MatrixSpec(
        group_axis="celltype",
        atlas="asp",
        specificity=SpecificitySpec(
            method="ewce_fraction", targets=("CM",), combine="max",
            name="cm_spec", emit="rollup",
        ),
    )
    df = reduce_matrix_to_gene(_summary_adata(), mspec)
    assert list(df.columns) == ["symbol", "asp_cm_spec"]


def test_emit_rollup_without_targets_is_rejected():
    """An empty target set summed to a silent all-zero column before this change."""
    from hvantk.algorithms.annotation.spec import SpecificitySpec

    with pytest.raises(ValueError, match="needs targets"):
        SpecificitySpec(method="ewce_fraction", targets=(), emit="rollup")


def test_unknown_emit_is_rejected():
    from hvantk.algorithms.annotation.spec import SpecificitySpec

    with pytest.raises(ValueError, match="emit must be"):
        SpecificitySpec(method="ewce_fraction", targets=("CM",), emit="scalar")
