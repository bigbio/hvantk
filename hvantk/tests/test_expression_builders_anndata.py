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
