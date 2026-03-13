"""
Consolidated Hail integration tests for HGC core operations.

Merged from: test_file_utils, test_gvcf_combiner, test_converters,
             test_mts_combiner, test_sort_matrix_table.
Tests are ordered to reflect the pipeline: combine → convert → combine MTs → sort.
"""

import random
import shutil
from pathlib import Path

import pytest
import hail as hl

from hvantk.hgc.file_utils import validate_vcfs_paths, sort_mts_cols
from hvantk.hgc.combiners import (
    combine_gvcfs,
    combine_matrix_table_rows,
    combine_matrix_table_cols,
)
from hvantk.hgc.converters import (
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
    GNOMAD_AVAILABLE,
)
from hvantk.hgc.constants import GVCF_EXTENSION
from hvantk.data.file_utils import compress_files, decompress_files

TESTS_DIR = Path(__file__).parent.parent / "testdata" / "hgc_data"


def test_validate_vcfs_paths():
    """Test GVCF file discovery from directory."""
    directory = TESTS_DIR / "gvcfs"
    vcfs = validate_vcfs_paths(str(directory), pattern=f"*{GVCF_EXTENSION}")
    assert len(vcfs) == 3


@pytest.mark.hail
@pytest.mark.order2
def test_combine_gvcfs():
    """Test combining GVCFs into a VDS."""
    gvcf_dir = TESTS_DIR / "gvcfs"
    vds_output_path = TESTS_DIR / "vds/cohort.vds"
    tmp_path = TESTS_DIR.parent / "local/tmp"
    plan_path = tmp_path / "combiner_plan.json"
    combine_gvcfs(
        gvcf_dir=str(gvcf_dir),
        vds_output_path=str(vds_output_path),
        tmp_path=str(tmp_path),
        save_path=str(plan_path),
        vdses=[],
        kwargs={},
    )
    assert vds_output_path.joinpath("reference_data/_SUCCESS").exists()
    assert vds_output_path.joinpath("variant_data/_SUCCESS").exists()

    compress_files(
        source_dir=str(vds_output_path),
        output_zip=str(vds_output_path.with_suffix(".vds.zip")),
        remove_originals=True,
    )
    if vds_output_path.exists():
        shutil.rmtree(vds_output_path)
    if tmp_path.exists():
        shutil.rmtree(tmp_path)


@pytest.mark.order3
@pytest.mark.skipif(not GNOMAD_AVAILABLE, reason="gnomad package not installed")
def test_convert_vds_to_mt():
    """Test VDS → MatrixTable conversion."""
    decompress_files(
        zip_path=str(TESTS_DIR / "vds/cohort.vds.zip"),
        extract_to=str(TESTS_DIR / "vds/cohort.vds"),
        remove_originals=False,
    )
    vds_path = TESTS_DIR / "vds/cohort.vds"
    mt_output_path = TESTS_DIR / "mts/cohort.mt"
    convert_vds_to_mt(
        vds_path=str(vds_path),
        output_path=str(mt_output_path),
        adjust_genotypes=True,
        skip_validation=False,
        skip_split_multi=False,
        skip_keying_by_cols=False,
        overwrite=True,
    )
    assert mt_output_path.joinpath("_SUCCESS").exists()

    compress_files(
        source_dir=str(mt_output_path),
        output_zip=str(mt_output_path.with_suffix(".mt.zip")),
        remove_originals=True,
    )
    if vds_path.exists():
        shutil.rmtree(vds_path)
    if mt_output_path.exists():
        shutil.rmtree(mt_output_path)


@pytest.mark.hail
@pytest.mark.order4
def test_convert_mt_to_cvcf():
    """Test MatrixTable → multi-sample VCF conversion."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(TESTS_DIR / "mts/cohort.mt"),
        remove_originals=False,
    )
    mt_path = TESTS_DIR / "mts/cohort.mt"
    vcf_output_path = TESTS_DIR / "vcf/cohort.vcf.bgz"
    convert_mt_to_multi_sample_vcf(
        mt_path=str(mt_path),
        vcf_path=str(vcf_output_path),
        filter_adj_genotypes=True,
        min_ac=1,
        split_multi=True,
    )
    assert vcf_output_path.exists()
    if mt_path.exists():
        shutil.rmtree(mt_path)


@pytest.mark.hail
@pytest.mark.order5
def test_combine_matrix_table_rows():
    """Test combining rows of two MatrixTables."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(TESTS_DIR / "mts/test.mt"),
        remove_originals=False,
    )
    mt1_path = TESTS_DIR / "mts/test.mt"
    mt_output_path = TESTS_DIR / "mts/combined.mt"
    combine_matrix_table_rows(
        mt_paths=[str(mt1_path), str(mt1_path)],
        output_path=str(mt_output_path),
        n_partitions=20,
        overwrite=True,
    )
    assert mt_output_path.joinpath("_SUCCESS").exists()
    assert (
        hl.read_matrix_table(str(mt_output_path)).count_rows()
        == 2 * hl.read_matrix_table(str(mt1_path)).count_rows()
    )
    shutil.rmtree(mt1_path)
    shutil.rmtree(mt_output_path)


@pytest.mark.hail
@pytest.mark.order6
def test_combine_matrix_table_cols():
    """Test combining columns of two MatrixTables."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(TESTS_DIR / "mts/test.mt"),
        remove_originals=False,
    )
    mt1_path = TESTS_DIR / "mts/test.mt"
    mt_output_path = TESTS_DIR / "mts/combined.mt"
    combine_matrix_table_cols(
        mt_paths=[str(mt1_path), str(mt1_path)],
        output_path=str(mt_output_path),
        n_partitions=20,
        overwrite=True,
    )
    assert mt_output_path.joinpath("_SUCCESS").exists()
    assert (
        hl.read_matrix_table(str(mt_output_path)).count_cols()
        == 2 * hl.read_matrix_table(str(mt1_path)).count_cols()
    )
    shutil.rmtree(mt1_path)
    shutil.rmtree(mt_output_path)


@pytest.mark.hail
@pytest.mark.order7
def test_sort_mts_cols():
    """Test sorting column order of MatrixTables for union_rows."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(TESTS_DIR / "mts/cohort.mt"),
        remove_originals=False,
    )
    mt1_path = TESTS_DIR / "mts/cohort.mt"
    mt1 = hl.read_matrix_table(str(mt1_path))

    idx = list(range(mt1.count_cols()))
    random.shuffle(idx)
    mt2 = mt1.choose_cols(idx)

    sorted_mts = sort_mts_cols([mt1, mt2])
    mt = hl.MatrixTable.union_rows(*sorted_mts)

    assert mt.count_rows() == mt1.count_rows() + mt2.count_rows()
    assert mt.count_cols() == mt1.count_cols() == mt2.count_cols()

    if mt1_path.exists():
        shutil.rmtree(mt1_path)
