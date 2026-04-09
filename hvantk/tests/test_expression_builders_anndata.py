"""Tests for AnnData-based expression builders (UCSC Cell Browser)."""

import anndata as ad
import numpy as np
import pytest


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

    def test_builds_anndata_from_tsv(self, tmp_path):
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

        from hvantk.tables.matrix_builders import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            output_path=out_path,
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
        # Column summary computed
        assert "column_summary" in adata.uns

    def test_no_output_path(self, tmp_path):
        """Build AnnData without writing to disk."""
        cells = ["cell_A", "cell_B"]
        genes = ["TP53", "BRCA1"]
        values = [[1.0, 2.0], [3.0, 4.0]]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia"])

        from hvantk.tables.matrix_builders import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
        )

        assert isinstance(adata, ad.AnnData)
        assert adata.shape == (2, 2)

    def test_split_gene_field(self, tmp_path):
        """Pipe-separated gene names are split to first element."""
        cells = ["cell_A", "cell_B"]
        genes = ["TP53|TP53L1", "BRCA1|BRCA1P1"]
        values = [
            [1.0, 2.0],
            [3.0, 4.0],
        ]

        expr_path = str(tmp_path / "expr.tsv")
        meta_path = str(tmp_path / "meta.tsv")

        self._write_expression_matrix(expr_path, genes, cells, values)
        self._write_metadata(meta_path, cells, ["neuron", "glia"])

        from hvantk.tables.matrix_builders import build_ucsc_ad

        adata = build_ucsc_ad(
            expression_matrix_path=expr_path,
            metadata_path=meta_path,
            split_gene_field=True,
        )

        assert "TP53" in adata.var.index.tolist()
        assert "BRCA1" in adata.var.index.tolist()
        # Pipe-separated versions should NOT be present
        assert "TP53|TP53L1" not in adata.var.index.tolist()


class TestBuildExpressionAtlasAd:
    """Tests for build_expression_atlas_ad in matrix_builders."""

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

        from hvantk.tables.matrix_builders import build_expression_atlas_ad

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

        from hvantk.tables.expression_atlas import create_anndata_from_expression_atlas

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

        from hvantk.tables.expression_atlas import create_anndata_from_expression_atlas

        adata = create_anndata_from_expression_atlas(expression_matrix_path=expr_path)

        # X should be (samples x genes), so X[0] = [1.0, 3.0] (s1 across G1, G2)
        np.testing.assert_array_almost_equal(adata.X[0], [1.0, 3.0])
        np.testing.assert_array_almost_equal(adata.X[1], [2.0, 4.0])
