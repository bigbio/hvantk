"""
Hail integration tests for expression table builders (UCSC + Expression Atlas).

Merged from: test_ucsc_matrix_table.py, test_expression_atlas_table.py.
"""

import shutil
import tempfile
from pathlib import Path

import pytest

from hvantk.core.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN
from hvantk.tables.ucsc import (
    convert_ucsc_metadata_to_hail_table,
    create_mt_from_ucsc_expression_matrix,
)
from hvantk.tables.expression_atlas import (
    convert_sdrf_to_hail_table,
    create_mt_from_expression_atlas_matrix,
)

pytestmark = [pytest.mark.hail, pytest.mark.slow]

TESTDATA_DIR = Path(__file__).parent / "testdata" / "raw"
UCSC_DIR = TESTDATA_DIR / "ucsc"
ATLAS_DIR = TESTDATA_DIR / "expression_atlas"


@pytest.fixture
def temp_dir():
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


# --- UCSC ---


def test_ucsc_metadata_to_hail_table(temp_dir):
    """Test conversion of UCSC metadata to Hail Table."""
    ht = convert_ucsc_metadata_to_hail_table(
        str(UCSC_DIR / "meta.test.tsv"), sep="\t", index_col=0, index_name=UCSC_CELL_ID_COLUMN
    )
    assert ht.count() == 9999


def test_ucsc_expression_matrix_to_mt(temp_dir):
    """Test UCSC expression matrix → MatrixTable with metadata."""
    output_path = Path(temp_dir) / "ucsc_expression_matrix.mt"
    metadata_ht = convert_ucsc_metadata_to_hail_table(
        str(UCSC_DIR / "meta.test.tsv"), sep="\t", index_col=0, index_name=UCSC_CELL_ID_COLUMN
    )
    mt = create_mt_from_ucsc_expression_matrix(
        expression_matrix_path=str(UCSC_DIR / "exprMatrix.test.tsv.bgz"),
        output_path=str(output_path),
        delimiter="\t",
        row_fields=None,
        row_key=UCSC_GENE_COLUMN,
        split_gene_field=True,
        min_partitions=5,
        force_bgz=True,
        overwrite=True,
        metadata_ht=metadata_ht,
    )
    assert mt.count_rows() == 199
    assert mt.count_cols() == 9999
    assert (output_path / "_SUCCESS").exists()


# --- Expression Atlas ---


def test_expression_atlas_sdrf_to_hail_table(temp_dir):
    """Test conversion of SDRF metadata to Hail Table."""
    output_path = Path(temp_dir) / "expression_atlas_table.ht"
    ht = convert_sdrf_to_hail_table(
        str(ATLAS_DIR / "E-MTAB-6798.condensed-sdrf.tsv"), output_file=str(output_path)
    )
    assert ht.count() == 317


def test_expression_atlas_matrix_to_mt(temp_dir):
    """Test Expression Atlas expression matrix → MatrixTable."""
    metadata_ht = convert_sdrf_to_hail_table(
        str(ATLAS_DIR / "E-MTAB-6798.condensed-sdrf.tsv"),
    )
    mt = create_mt_from_expression_atlas_matrix(
        expression_matrix_path=str(ATLAS_DIR / "E-MTAB-6798-transcripts-tpms.tsv.bgz"),
        metadata_ht=metadata_ht,
    )
    assert mt.count_rows() == 116643
    assert mt.count_cols() == 317
