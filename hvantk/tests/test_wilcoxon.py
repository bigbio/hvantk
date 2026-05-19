"""Unit tests for Wilcoxon rank-sum marker gene detection.

All tests use synthetic numpy data — no Hail dependency.
Runs in the default ``pytest -q`` configuration.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

from hvantk.algorithms.statistics.wilcoxon import (
    WilcoxonParams,
    _compute_rank_matrix,
    _compute_tie_correction,
    _wilcoxon_one_vs_rest,
    rank_genes_groups,
    results_to_gene_set_collection,
)
from hvantk.core.utils.gene_sets import GeneSetCollection


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def two_group_data():
    """Two groups (A=30, B=20) with 6 genes.

    G0, G1 are markers for A (high in A, low in B).
    G2, G3 are markers for B (high in B, low in A).
    G4 is housekeeping.  G5 is noise.
    """
    rng = np.random.RandomState(42)
    n_A, n_B, n_genes = 30, 20, 6
    n_cells = n_A + n_B

    X = np.zeros((n_cells, n_genes), dtype=np.float64)
    # A-markers
    X[:n_A, 0] = rng.uniform(5, 10, n_A)
    X[:n_A, 1] = rng.uniform(4, 8, n_A)
    X[n_A:, 0] = rng.uniform(0, 1, n_B)
    X[n_A:, 1] = rng.uniform(0, 1, n_B)
    # B-markers
    X[:n_A, 2] = rng.uniform(0, 1, n_A)
    X[:n_A, 3] = rng.uniform(0, 1, n_A)
    X[n_A:, 2] = rng.uniform(5, 10, n_B)
    X[n_A:, 3] = rng.uniform(4, 8, n_B)
    # Housekeeping
    X[:, 4] = rng.uniform(3, 5, n_cells)
    # Noise
    X[:, 5] = rng.uniform(0, 0.1, n_cells)

    labels = np.array(["A"] * n_A + ["B"] * n_B)
    gene_ids = np.array([f"ENSG{i:04d}" for i in range(n_genes)])
    gene_names = np.array(["G0", "G1", "G2", "G3", "G4", "G5"])

    return X, labels, gene_ids, gene_names


@pytest.fixture()
def three_group_data():
    """Three groups (X=20, Y=15, Z=15) with 4 genes."""
    rng = np.random.RandomState(123)
    nX, nY, nZ, n_genes = 20, 15, 15, 4
    n_cells = nX + nY + nZ

    X = np.zeros((n_cells, n_genes), dtype=np.float64)
    # Gene 0: marker for X
    X[:nX, 0] = rng.uniform(5, 10, nX)
    X[nX:, 0] = rng.uniform(0, 1, nY + nZ)
    # Gene 1: marker for Y
    X[nX : nX + nY, 1] = rng.uniform(5, 10, nY)
    X[:nX, 1] = rng.uniform(0, 1, nX)
    X[nX + nY :, 1] = rng.uniform(0, 1, nZ)
    # Gene 2: marker for Z
    X[nX + nY :, 2] = rng.uniform(5, 10, nZ)
    X[: nX + nY, 2] = rng.uniform(0, 1, nX + nY)
    # Gene 3: housekeeping
    X[:, 3] = rng.uniform(3, 5, n_cells)

    labels = np.array(["X"] * nX + ["Y"] * nY + ["Z"] * nZ)
    gene_ids = np.array([f"ENSG{i:04d}" for i in range(n_genes)])
    gene_names = np.array(["GX", "GY", "GZ", "GH"])

    return X, labels, gene_ids, gene_names


# ---------------------------------------------------------------------------
# TestComputeRankMatrix
# ---------------------------------------------------------------------------


class TestComputeRankMatrix:
    def test_no_ties(self):
        X = np.array([[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]])
        ranks = _compute_rank_matrix(X)
        np.testing.assert_array_equal(ranks[:, 0], [1, 2, 3, 4])
        np.testing.assert_array_equal(ranks[:, 1], [4, 3, 2, 1])

    def test_ties(self):
        X = np.array([[1.0], [1.0], [3.0], [3.0]])
        ranks = _compute_rank_matrix(X)
        # average tie-breaking: (1+2)/2=1.5, (3+4)/2=3.5
        np.testing.assert_array_equal(ranks[:, 0], [1.5, 1.5, 3.5, 3.5])

    def test_zero_heavy(self):
        """Many zeros should produce correct average ranks."""
        X = np.array([[0.0], [0.0], [0.0], [1.0], [2.0]])
        ranks = _compute_rank_matrix(X)
        # zeros get average rank (1+2+3)/3 = 2.0
        np.testing.assert_allclose(ranks[:3, 0], [2.0, 2.0, 2.0])
        assert ranks[3, 0] == 4.0
        assert ranks[4, 0] == 5.0

    def test_shape_preserved(self):
        X = np.random.rand(10, 5)
        ranks = _compute_rank_matrix(X)
        assert ranks.shape == X.shape


# ---------------------------------------------------------------------------
# TestComputeTieCorrection
# ---------------------------------------------------------------------------


class TestComputeTieCorrection:
    def test_no_ties_gives_one(self):
        X = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])
        tc = _compute_tie_correction(X, 3)
        np.testing.assert_allclose(tc, [1.0, 1.0])

    def test_all_same_gives_zero(self):
        X = np.array([[5.0], [5.0], [5.0], [5.0]])
        tc = _compute_tie_correction(X, 4)
        np.testing.assert_allclose(tc, [0.0])

    def test_partial_ties(self):
        # 3 cells: two share a value, one different
        X = np.array([[1.0], [1.0], [2.0]])
        tc = _compute_tie_correction(X, 3)
        # t=2 for the tie group, t=1 for the singleton
        # correction = 1 - ((8-2) + (1-1)) / (27-3) = 1 - 6/24 = 0.75
        np.testing.assert_allclose(tc, [0.75])


# ---------------------------------------------------------------------------
# TestWilcoxonOneVsRest
# ---------------------------------------------------------------------------


class TestWilcoxonOneVsRest:
    def test_perfect_separation(self):
        """Groups with no overlap should give very small p-values."""
        X = np.array(
            [
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],  # group 0
                [10.0],
                [10.0],
                [10.0],
                [10.0],
                [10.0],
            ]  # group 1
        )
        ranks = _compute_rank_matrix(X)
        tc = _compute_tie_correction(X, 10)
        mask = np.array([False] * 5 + [True] * 5)

        U, z, p = _wilcoxon_one_vs_rest(ranks, mask, 10, tc)
        assert p[0] < 0.01

    def test_identical_groups(self):
        """Identical distributions → p ≈ 1."""
        rng = np.random.RandomState(99)
        X = rng.uniform(0, 1, (40, 3))
        ranks = _compute_rank_matrix(X)
        tc = _compute_tie_correction(X, 40)
        mask = np.array([True] * 20 + [False] * 20)

        _, _, p = _wilcoxon_one_vs_rest(ranks, mask, 40, tc)
        # p-values should be large (no significant difference)
        assert np.all(p > 0.05)

    def test_versus_scipy(self, two_group_data):
        """Compare vectorised implementation against scipy.stats.mannwhitneyu."""
        X, labels, _, _ = two_group_data
        n_cells = X.shape[0]
        mask = labels == "A"

        ranks = _compute_rank_matrix(X)
        tc = _compute_tie_correction(X, n_cells)
        U_ours, _, p_ours = _wilcoxon_one_vs_rest(ranks, mask, n_cells, tc)

        for j in range(X.shape[1]):
            u_scipy, p_scipy = mannwhitneyu(
                X[mask, j], X[~mask, j], alternative="two-sided"
            )
            # scipy returns U for sample x (first arg) directly
            np.testing.assert_allclose(U_ours[j], u_scipy, rtol=1e-6)
            # p-values should be in the same ballpark (normal approx vs exact)
            if p_scipy < 0.01:
                assert p_ours[j] < 0.05  # both small
            if p_scipy > 0.5:
                assert p_ours[j] > 0.01  # both large

    def test_empty_group(self):
        """Empty focal group should return NaN/1.0."""
        X = np.random.rand(10, 2)
        ranks = _compute_rank_matrix(X)
        tc = _compute_tie_correction(X, 10)
        mask = np.zeros(10, dtype=bool)

        U, z, p = _wilcoxon_one_vs_rest(ranks, mask, 10, tc)
        assert np.all(np.isnan(U))
        assert np.all(p == 1.0)


# ---------------------------------------------------------------------------
# TestRankGenesGroups
# ---------------------------------------------------------------------------


class TestRankGenesGroups:
    def test_two_groups(self, two_group_data):
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0, min_fraction_expressed=0.0, alpha=0.05
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)

        assert isinstance(results, pd.DataFrame)
        assert "group" in results.columns
        assert "gene_id" in results.columns
        assert "gene_name" in results.columns
        assert "pvalue_adj" in results.columns
        assert "fold_change" in results.columns
        assert len(results) > 0

    def test_three_groups(self, three_group_data):
        X, labels, gene_ids, gene_names = three_group_data
        params = WilcoxonParams(min_fold_change=1.0, min_fraction_expressed=0.0)
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        groups_found = set(results["group"].unique())
        assert groups_found == {"X", "Y", "Z"}

    def test_markers_make_sense(self, two_group_data):
        """G0/G1 should be top markers for A; G2/G3 for B."""
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0, min_fraction_expressed=0.0, alpha=0.05
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)

        a_sig = results[(results["group"] == "A") & (results["pvalue_adj"] <= 0.05)]
        assert "G0" in a_sig["gene_name"].values
        assert "G1" in a_sig["gene_name"].values

        b_sig = results[(results["group"] == "B") & (results["pvalue_adj"] <= 0.05)]
        assert "G2" in b_sig["gene_name"].values
        assert "G3" in b_sig["gene_name"].values

    def test_output_columns_no_gene_names(self, two_group_data):
        X, labels, gene_ids, _ = two_group_data
        params = WilcoxonParams(min_fold_change=1.0, min_fraction_expressed=0.0)
        results = rank_genes_groups(X, labels, gene_ids, gene_names=None, params=params)
        assert "gene_name" not in results.columns
        assert "gene_id" in results.columns

    def test_all_genes_tested(self, two_group_data):
        """All input genes should be tested for all groups (no per-group pre-filter)."""
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams()
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        n_groups = len(np.unique(labels))
        n_genes = len(gene_ids)
        assert len(results) == n_groups * n_genes

    def test_skip_small_group(self):
        """Groups with < 2 cells should be skipped → empty if all too small."""
        X = np.random.rand(3, 2)
        labels = np.array(["A", "B", "C"])  # 1 cell each
        gene_ids = np.array(["G0", "G1"])
        results = rank_genes_groups(X, labels, gene_ids)
        assert len(results) == 0
        assert "group" in results.columns

    def test_n_total_genes_correction(self, two_group_data):
        """Passing n_total_genes should make adjusted p-values more
        conservative (Seurat-style: p.adjust uses total gene count)."""
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        # Without n_total_genes (correction uses n_candidates only)
        results_default = rank_genes_groups(
            X,
            labels,
            gene_ids,
            gene_names,
            params,
        )
        # With n_total_genes >> n_genes (simulating pre-filtered subset)
        results_total = rank_genes_groups(
            X,
            labels,
            gene_ids,
            gene_names,
            params,
            n_total_genes=30000,
        )

        # Merge on (group, gene_id) to compare adjusted p-values
        merged = results_default.merge(
            results_total,
            on=["group", "gene_id"],
            suffixes=("_default", "_total"),
        )
        assert len(merged) > 0
        # Raw p-values should be identical
        np.testing.assert_allclose(merged["pvalue_default"], merged["pvalue_total"])
        # Adjusted p-values with n_total_genes should be >= default
        assert (
            merged["pvalue_adj_total"] >= merged["pvalue_adj_default"] - 1e-12
        ).all()


# ---------------------------------------------------------------------------
# TestResultsToGeneSetCollection
# ---------------------------------------------------------------------------


class TestResultsToGeneSetCollection:
    def test_basic_conversion(self, two_group_data):
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            alpha=0.05,
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        bg = set(gene_names)

        collection = results_to_gene_set_collection(
            results,
            bg,
            top_n=10,
            alpha=0.05,
            gene_col="gene_name",
        )
        assert isinstance(collection, GeneSetCollection)
        assert len(collection.background_genes) == len(bg)

    def test_alpha_filtering(self, two_group_data):
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        bg = set(gene_names)

        strict = results_to_gene_set_collection(
            results,
            bg,
            top_n=100,
            alpha=1e-10,
        )
        lenient = results_to_gene_set_collection(
            results,
            bg,
            top_n=100,
            alpha=0.5,
        )
        strict_total = sum(gs.n_genes for gs in strict)
        lenient_total = sum(gs.n_genes for gs in lenient)
        assert strict_total <= lenient_total

    def test_top_n(self, two_group_data):
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        bg = set(gene_names)

        collection = results_to_gene_set_collection(
            results,
            bg,
            top_n=1,
            alpha=1.0,
        )
        for gs in collection:
            assert gs.n_genes <= 1

    def test_fallback_to_gene_id(self, two_group_data):
        X, labels, gene_ids, _ = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names=None, params=params)
        bg = set(gene_ids)

        # gene_col="gene_name" not in results → should fall back to gene_id
        collection = results_to_gene_set_collection(
            results,
            bg,
            top_n=10,
            alpha=1.0,
            gene_col="gene_name",
        )
        # Verify genes are gene_ids
        for gs in collection:
            assert all(g.startswith("ENSG") for g in gs.genes)

    def test_metadata_fields(self, two_group_data):
        X, labels, gene_ids, gene_names = two_group_data
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        results = rank_genes_groups(X, labels, gene_ids, gene_names, params)
        bg = set(gene_names)

        collection = results_to_gene_set_collection(
            results,
            bg,
            top_n=10,
            alpha=1.0,
        )
        for gs in collection:
            assert "fold_changes" in gs.metadata
            assert "adjusted_pvalues" in gs.metadata
            assert gs.metadata["method"] == "wilcoxon"
