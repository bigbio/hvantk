"""Tests for AnnData-based expression helpers and visualizations.

The Phase A ``build_*_ad`` builders (Expression Atlas, CPTAC, UCSC Cell
Browser) were retired with issue #114; the per-skill tests now exercise the
Phase B ``build_<x>_<dataset>`` builders via ``run_builder_for_spec`` (see
``hvantk/skills/<plugin>/tests/test_builder.py``). What stays here:

- Coverage for ``create_anndata_from_expression_atlas`` (a shared helper
  used by the Expression Atlas Phase B builder) — sanity checks that the
  TSV→AnnData transposition is correct.
- Coverage for ``visualize_expression_distribution`` — generic AnnData
  visualization utility, unrelated to the build path.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


class TestExpressionAtlasHelpers:
    """Tests for the shared ``create_anndata_from_expression_atlas`` helper."""

    def _write_expression_tsv(self, path, gene_ids, gene_names, samples, values):
        """Write an Expression Atlas-style expression TSV."""
        with open(path, "w") as fh:
            fh.write("Gene ID\tGene Name\t" + "\t".join(samples) + "\n")
            for gid, gname, row in zip(gene_ids, gene_names, values):
                fh.write(gid + "\t" + gname + "\t" + "\t".join(str(v) for v in row) + "\n")

    def test_builds_without_sdrf_metadata(self, tmp_path):
        """Build AnnData from expression TSV only, no SDRF."""
        gene_ids = ["ENSG00000000001", "ENSG00000000002"]
        gene_names = ["GENE_A", "GENE_B"]
        samples = ["sample_1", "sample_2"]
        values = [[1.0, 2.0], [3.0, 4.0]]

        expr_path = str(tmp_path / "expr.tsv")
        self._write_expression_tsv(expr_path, gene_ids, gene_names, samples, values)

        from hvantk.skills.expression_atlas.shared.expression_atlas import (
            create_anndata_from_expression_atlas,
        )

        adata = create_anndata_from_expression_atlas(
            expression_matrix_path=expr_path,
        )

        assert isinstance(adata, ad.AnnData)
        assert adata.shape == (2, 2)
        assert set(adata.obs.index) == set(samples)
        assert set(adata.var.index) == set(gene_ids)

    def test_expression_values_correct(self, tmp_path):
        """Verify expression values are correctly transposed."""
        gene_ids = ["G1", "G2"]
        gene_names = ["Gene1", "Gene2"]
        samples = ["s1", "s2"]
        values = [[1.0, 2.0], [3.0, 4.0]]

        expr_path = str(tmp_path / "expr.tsv")
        self._write_expression_tsv(expr_path, gene_ids, gene_names, samples, values)

        from hvantk.skills.expression_atlas.shared.expression_atlas import (
            create_anndata_from_expression_atlas,
        )

        adata = create_anndata_from_expression_atlas(expression_matrix_path=expr_path)

        # X should be (samples x genes), so X[0] = [1.0, 3.0] (s1 across G1, G2)
        np.testing.assert_array_almost_equal(adata.X[0], [1.0, 3.0])
        np.testing.assert_array_almost_equal(adata.X[1], [2.0, 4.0])


class TestVisualizeExpressionAd:
    def test_returns_matplotlib_figure(self):
        from hvantk.algorithms.visualization.expression.anndata import (
            visualize_expression_distribution,
        )
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend
        import matplotlib.pyplot as plt

        adata = ad.AnnData(
            X=np.random.rand(50, 20).astype(np.float32),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(50)]),
            var=pd.DataFrame(index=[f"g_{i}" for i in range(20)]),
        )
        fig = visualize_expression_distribution(adata)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_returns_matplotlib_figure_sparse(self):
        from hvantk.algorithms.visualization.expression.anndata import (
            visualize_expression_distribution,
        )
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        adata = ad.AnnData(
            X=sparse.csr_matrix(np.array([[0.0, 1.0], [0.0, 2.0]], dtype=np.float32)),
            obs=pd.DataFrame(index=["c1", "c2"]),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        fig = visualize_expression_distribution(adata)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
