"""Tests for AnnData-based expression builders (UCSC, Expression Atlas, CPTAC)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.mark.parametrize("backed", [False, True])
class TestBuildUcscAd:
    """Tests for build_ucsc_ad in matrix_builders."""

    def _write_expression_matrix(self, path, genes, cells, values):
        """Write a genes-x-cells TSV expression matrix."""
        with open(path, "w") as fh:
            fh.write("gene\t" + "\t".join(cells) + "\n")
            for gene, row in zip(genes, values):
                fh.write(gene + "\t" + "\t".join(str(v) for v in row) + "\n")

    def _write_metadata(self, path, cell_ids, cell_types):
        """Write a metadata TSV with cell_id index and cell_type column."""
        with open(path, "w") as fh:
            fh.write("cell_id\tcell_type\n")
            for cid, ct in zip(cell_ids, cell_types):
                fh.write(f"{cid}\t{ct}\n")

    def test_builds_anndata_from_tsv(self, tmp_path, backed):
        """Build AnnData from minimal expression matrix + metadata."""
        cells = ["cell_A", "cell_B", "cell_C", "cell_D"]
        genes = ["TP53", "BRCA1", "EGFR"]
        values = [
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
        ]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        out_path = str(tmp_path / "output.h5ad")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia", "neuron", "glia"])

        from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            output_path=out_path,
            backed=backed,
            overwrite=True,
        )

        assert isinstance(adata, ad.AnnData)
        # 4 cells (obs) x 3 genes (var)
        assert adata.shape == (4, 3)
        assert "cell_type" in adata.obs.columns
        assert set(adata.obs.index) == set(cells)
        assert set(adata.var.index) == set(genes)
        # Output file written
        assert (tmp_path / "output.h5ad").exists()
        # Provenance metadata stored
        assert "hvantk_metadata" in adata.uns
        assert adata.uns["hvantk_metadata"]["source_name"] == "UCSC"
        # Column summary computed (skipped in backed mode by design)
        if not backed:
            assert "column_summary" in adata.uns

    def test_minimal_two_by_two(self, tmp_path, backed):
        """Build AnnData from a minimal 2x2 matrix and verify shape."""
        cells = ["cell_A", "cell_B"]
        genes = ["TP53", "BRCA1"]
        values = [[1.0, 2.0], [3.0, 4.0]]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        out_path = str(tmp_path / "output.h5ad")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia"])

        from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            output_path=out_path,
            backed=backed,
            overwrite=True,
        )

        assert isinstance(adata, ad.AnnData)
        assert adata.shape == (2, 2)

    def test_split_gene_field(self, tmp_path, backed):
        """Pipe-separated gene names are split to first element."""
        cells = ["cell_A", "cell_B"]
        genes = ["TP53|TP53L1", "BRCA1|BRCA1P1"]
        values = [
            [1.0, 2.0],
            [3.0, 4.0],
        ]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        out_path = str(tmp_path / "output.h5ad")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia"])

        from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            split_gene_field=True,
            output_path=out_path,
            backed=backed,
            overwrite=True,
        )

        assert "TP53" in adata.var.index.tolist()
        assert "BRCA1" in adata.var.index.tolist()
        # Pipe-separated versions should NOT be present
        assert "TP53|TP53L1" not in adata.var.index.tolist()

    def test_ucsc_respects_delimiter_for_metadata(self, tmp_path, backed):
        expr_path = str(tmp_path / "expr.csv")
        meta_path = str(tmp_path / "meta.csv")
        out_path = str(tmp_path / "output.h5ad")
        with open(expr_path, "w") as fh:
            fh.write("gene,cell_A,cell_B\n")
            fh.write("TP53,1.0,2.0\n")
            fh.write("BRCA1,3.0,4.0\n")
        with open(meta_path, "w") as fh:
            fh.write("cell_id,cell_type\n")
            fh.write("cell_A,neuron\n")
            fh.write("cell_B,glia\n")

        from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            delimiter=",",
            output_path=out_path,
            backed=backed,
            overwrite=True,
        )
        assert adata.shape == (2, 2)
        assert "cell_type" in adata.obs.columns


class TestBuildExpressionAtlasAd:
    """Tests for build_expression_atlas_ad in hvantk.skills.expression_atlas.builder."""

    def _write_expression_tsv(self, path, gene_ids, gene_names, samples, values):
        """Write an Expression Atlas-style expression TSV."""
        with open(path, "w") as fh:
            fh.write("Gene ID\tGene Name\t" + "\t".join(samples) + "\n")
            for gid, gname, row in zip(gene_ids, gene_names, values):
                fh.write(gid + "\t" + gname + "\t" + "\t".join(str(v) for v in row) + "\n")

    def _write_sdrf(self, path, rows):
        """Write an SDRF file (tab-delimited, no header).

        Each row is a tuple of (accession, unused, sample_id, column_type, column_name, column_value).
        """
        with open(path, "w") as fh:
            for row in rows:
                fh.write("\t".join(str(c) for c in row) + "\n")

    def test_builds_anndata_from_tsv_and_sdrf(self, tmp_path):
        """Build AnnData from Expression Atlas expression TSV + SDRF metadata."""
        gene_ids = ["ENSG00000000001", "ENSG00000000002"]
        gene_names = ["GENE_A", "GENE_B"]
        samples = ["sample_1", "sample_2"]
        values = [
            [1.5, 2.5],
            [3.5, 4.5],
        ]

        expr_path = str(tmp_path / "expr.tsv")
        sdrf_path = str(tmp_path / "sdrf.tsv")
        out_path = str(tmp_path / "output.h5ad")

        self._write_expression_tsv(expr_path, gene_ids, gene_names, samples, values)
        self._write_sdrf(sdrf_path, [
            ("E-MTAB-0001", "", "sample_1", "characteristic", "organism", "Homo sapiens"),
            ("E-MTAB-0001", "", "sample_1", "characteristic", "tissue", "brain"),
            ("E-MTAB-0001", "", "sample_2", "characteristic", "organism", "Homo sapiens"),
            ("E-MTAB-0001", "", "sample_2", "characteristic", "tissue", "liver"),
        ])

        from hvantk.skills.expression_atlas.builder import build_expression_atlas_ad

        adata = build_expression_atlas_ad(
            expression_matrix_path=expr_path,
            sdrf_file=sdrf_path,
            output_path=out_path,
        )

        assert isinstance(adata, ad.AnnData)
        # 2 samples (obs) x 2 genes (var)
        assert adata.shape == (2, 2)
        # Gene Name stored in var
        assert "Gene Name" in adata.var.columns
        assert adata.var["Gene Name"].tolist() == gene_names
        # Output file written
        assert (tmp_path / "output.h5ad").exists()
        # Provenance metadata
        assert "hvantk_metadata" in adata.uns
        assert adata.uns["hvantk_metadata"]["source_name"] == "ExpressionAtlas"
        # Column summary computed
        assert "column_summary" in adata.uns
        # SDRF metadata joined into obs
        assert "organism" in adata.obs.columns
        assert "tissue" in adata.obs.columns

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


class TestBuildCptacAd:
    """Tests for build_cptac_ad in the cptac expression plugin."""

    def test_builds_anndata_from_long_format(self, tmp_path):
        """Build AnnData from long-format CPTAC expression + metadata."""
        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        out_path = str(tmp_path / "output.h5ad")

        # Long-format expression: GeneID, Gene Name, SampleID, Expression
        with open(expr_path, "w") as fh:
            fh.write("GeneID\tGene Name\tSampleID\tExpression\n")
            fh.write("G1\tTP53\tS1\t1.5\n")
            fh.write("G1\tTP53\tS2\t2.5\n")
            fh.write("G2\tBRCA1\tS1\t3.5\n")
            fh.write("G2\tBRCA1\tS2\t4.5\n")

        # Metadata with tumor_type
        with open(meta_path, "w") as fh:
            fh.write("SampleID\ttumor_type\n")
            fh.write("S1\tLUAD\n")
            fh.write("S2\tBRCA\n")

        from hvantk.skills.cptac.expression.builder import build_cptac_ad

        adata = build_cptac_ad(
            expression_path=expr_path,
            metadata_path=meta_path,
            output_path=out_path,
        )

        assert isinstance(adata, ad.AnnData)
        # 2 samples x 2 genes
        assert adata.shape == (2, 2)
        assert "tumor_type" in adata.obs.columns
        assert (tmp_path / "output.h5ad").exists()
        assert "hvantk_metadata" in adata.uns
        assert adata.uns["hvantk_metadata"]["source_name"] == "CPTAC"
        assert "column_summary" in adata.uns

    def test_builds_from_csv_long_format(self, tmp_path):
        expr_path = str(tmp_path / "expr.csv")
        meta_path = str(tmp_path / "meta.csv")

        with open(expr_path, "w") as fh:
            fh.write("GeneID,Gene Name,SampleID,Expression\n")
            fh.write("G1,TP53,S1,1.5\n")
            fh.write("G1,TP53,S2,2.5\n")
            fh.write("G2,BRCA1,S1,3.5\n")
            fh.write("G2,BRCA1,S2,4.5\n")

        with open(meta_path, "w") as fh:
            fh.write("SampleID,tumor_type\n")
            fh.write("S1,LUAD\n")
            fh.write("S2,BRCA\n")

        from hvantk.skills.cptac.expression.builder import build_cptac_ad

        adata = build_cptac_ad(expression_path=expr_path, metadata_path=meta_path)
        assert adata.shape == (2, 2)


class TestBuildCptacPhosphoAd:
    """Tests for build_cptac_phospho_ad in the cptac phospho plugin."""

    def test_builds_anndata_from_sites_matrix(self, tmp_path):
        """Build AnnData from wide-format CPTAC phospho matrix + metadata."""
        expr_path = str(tmp_path / "phospho.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        out_path = str(tmp_path / "output.h5ad")

        # Wide format: SiteID, S1, S2
        with open(expr_path, "w") as fh:
            fh.write("SiteID\tS1\tS2\n")
            fh.write("TP53_S315\t100.0\t200.0\n")
            fh.write("EGFR_Y1068\t300.0\t400.0\n")

        # Metadata with tumor_type
        with open(meta_path, "w") as fh:
            fh.write("SampleID\ttumor_type\n")
            fh.write("S1\tLUAD\n")
            fh.write("S2\tBRCA\n")

        from hvantk.skills.cptac.phospho.builder import build_cptac_phospho_ad

        adata = build_cptac_phospho_ad(
            expression_path=expr_path,
            metadata_path=meta_path,
            output_path=out_path,
        )

        assert isinstance(adata, ad.AnnData)
        # 2 samples x 2 sites
        assert adata.shape == (2, 2)
        assert "gene_symbol" in adata.var.columns
        assert "tumor_type" in adata.obs.columns
        assert (tmp_path / "output.h5ad").exists()
        assert "hvantk_metadata" in adata.uns
        assert adata.uns["hvantk_metadata"]["source_name"] == "CPTAC"

    def test_builds_anndata_from_csv_sites_matrix(self, tmp_path):
        expr_path = str(tmp_path / "phospho.csv")
        meta_path = str(tmp_path / "meta.csv")

        with open(expr_path, "w") as fh:
            fh.write("SiteID,S1,S2\n")
            fh.write("TP53_S315,100.0,200.0\n")
            fh.write("EGFR_Y1068,300.0,400.0\n")

        with open(meta_path, "w") as fh:
            fh.write("SampleID,tumor_type\n")
            fh.write("S1,LUAD\n")
            fh.write("S2,BRCA\n")

        from hvantk.skills.cptac.phospho.builder import build_cptac_phospho_ad

        adata = build_cptac_phospho_ad(expression_path=expr_path, metadata_path=meta_path)
        assert adata.shape == (2, 2)


class TestMkmatrixCli:
    """CLI integration tests for mkmatrix commands."""

    def _write_expression_matrix(self, path, genes, cells, values):
        """Write a genes-x-cells TSV expression matrix."""
        with open(path, "w") as fh:
            fh.write("gene\t" + "\t".join(cells) + "\n")
            for gene, row in zip(genes, values):
                fh.write(gene + "\t" + "\t".join(str(v) for v in row) + "\n")

    def _write_metadata(self, path, cell_ids, cell_types):
        """Write a metadata TSV with cell_id index and cell_type column."""
        with open(path, "w") as fh:
            fh.write("cell_id\tcell_type\n")
            for cid, ct in zip(cell_ids, cell_types):
                fh.write(f"{cid}\t{ct}\n")

    def test_ucsc_produces_h5ad(self, tmp_path):
        from click.testing import CliRunner
        from hvantk.tools.build.make_matrix_cli import mkmatrix_group

        cells = ["cell_A", "cell_B", "cell_C", "cell_D"]
        genes = ["TP53", "BRCA1", "EGFR"]
        values = [
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
        ]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        output = str(tmp_path / "output.h5ad")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia", "neuron", "glia"])

        runner = CliRunner()
        result = runner.invoke(
            mkmatrix_group,
            ["ucsc", "-e", expr_path, "-m", meta_path, "-o", output],
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "output.h5ad").exists()

        # Verify it's a valid h5ad
        loaded = ad.read_h5ad(output)
        assert loaded.shape == (4, 3)

    def test_cptac_produces_h5ad(self, tmp_path):
        from click.testing import CliRunner
        from hvantk.tools.build.make_matrix_cli import mkmatrix_group

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")
        output = str(tmp_path / "output.h5ad")

        with open(expr_path, "w") as fh:
            fh.write("GeneID\tGene Name\tSampleID\tExpression\n")
            fh.write("G1\tTP53\tS1\t1.5\n")
            fh.write("G1\tTP53\tS2\t2.5\n")
            fh.write("G2\tBRCA1\tS1\t3.5\n")
            fh.write("G2\tBRCA1\tS2\t4.5\n")

        with open(meta_path, "w") as fh:
            fh.write("SampleID\ttumor_type\n")
            fh.write("S1\tLUAD\n")
            fh.write("S2\tBRCA\n")

        runner = CliRunner()
        result = runner.invoke(
            mkmatrix_group,
            ["cptac", "-e", expr_path, "-m", meta_path, "-o", output],
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "output.h5ad").exists()

        loaded = ad.read_h5ad(output)
        assert loaded.shape == (2, 2)


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
