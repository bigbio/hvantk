"""Tests for Hail <-> AnnData converters."""

import numpy as np
import pytest


@pytest.fixture()
def simple_mt(hail_session):
    """Build a minimal MatrixTable: 5 rows x 3 cols."""
    hl = hail_session
    mt = hl.utils.range_matrix_table(5, 3)
    mt = mt.annotate_rows(gene_id="GENE_" + hl.str(mt.row_idx))
    mt = mt.annotate_cols(sample_id="SAMPLE_" + hl.str(mt.col_idx))
    mt = mt.annotate_entries(x=hl.float64(mt.row_idx * 10 + mt.col_idx))
    mt = mt.key_rows_by("gene_id")
    mt = mt.key_cols_by("sample_id")
    return mt


@pytest.mark.hail
class TestHailMtToAnndata:
    """Tests for hail_mt_to_anndata."""

    def test_shape(self, simple_mt):
        from hvantk.core.utils.converters import hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        # AnnData convention: obs=samples (cols), var=genes (rows)
        assert adata.shape == (3, 5)

    def test_obs_index(self, simple_mt):
        from hvantk.core.utils.converters import hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        expected = {"SAMPLE_0", "SAMPLE_1", "SAMPLE_2"}
        assert set(adata.obs.index) == expected

    def test_var_index(self, simple_mt):
        from hvantk.core.utils.converters import hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        expected = {"GENE_0", "GENE_1", "GENE_2", "GENE_3", "GENE_4"}
        assert set(adata.var.index) == expected

    def test_values(self, simple_mt):
        from hvantk.core.utils.converters import hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        # Entry value = row_idx * 10 + col_idx
        # In AnnData: obs=samples (col_idx), var=genes (row_idx)
        # So adata.X[sample_i, gene_j] = gene_j_row_idx * 10 + sample_i_col_idx
        # SAMPLE_0 is col_idx=0, GENE_2 is row_idx=2 => value = 20
        sample_0_idx = list(adata.obs.index).index("SAMPLE_0")
        gene_2_idx = list(adata.var.index).index("GENE_2")
        assert adata.X[sample_0_idx, gene_2_idx] == pytest.approx(20.0)

    def test_dtype_float32(self, simple_mt):
        from hvantk.core.utils.converters import hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        assert adata.X.dtype == np.float32


@pytest.mark.hail
class TestAnndataToHailMt:
    """Tests for anndata_to_hail_mt (roundtrip)."""

    def test_roundtrip_shape(self, simple_mt, hail_session):
        from hvantk.core.utils.converters import anndata_to_hail_mt, hail_mt_to_anndata

        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        mt2 = anndata_to_hail_mt(
            adata, row_key="gene_id", col_key="sample_id", entry_field="x"
        )
        assert mt2.count() == (5, 3)

    def test_roundtrip_values(self, simple_mt, hail_session):
        from hvantk.core.utils.converters import anndata_to_hail_mt, hail_mt_to_anndata

        hl = hail_session
        adata = hail_mt_to_anndata(simple_mt, entry_field="x")
        mt2 = anndata_to_hail_mt(
            adata, row_key="gene_id", col_key="sample_id", entry_field="x"
        )
        # Check a specific value: GENE_3, SAMPLE_1 => 3*10+1 = 31
        row = mt2.filter_rows(mt2.gene_id == "GENE_3")
        row = row.filter_cols(row.sample_id == "SAMPLE_1")
        val = row.x.collect()[0]
        assert val == pytest.approx(31.0, abs=0.1)
