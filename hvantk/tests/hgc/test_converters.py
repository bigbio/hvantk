import pytest
import shutil
from pathlib import Path

from hvantk.hgc.converters import (
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
    GNOMAD_AVAILABLE,
)
from hvantk.data.file_utils import compress_files, decompress_files

TESTS_DIR = Path(__file__).parent.parent / "testdata" / "hgc_data"


# @pytest.mark.hail
@pytest.mark.order3
@pytest.mark.skipif(not GNOMAD_AVAILABLE, reason="gnomad package not installed")
def test_convert_vds_to_mt():
    """
    Test the convert_vds_to_mt function by converting a VDS to a MatrixTable.

    This test uses the session-scoped hail_session fixture for Hail initialization.
    """
    # Decompress the VDS zip file
    decompress_files(
        zip_path=str(TESTS_DIR / "vds/cohort.vds.zip"),
        extract_to=str(TESTS_DIR / "vds/cohort.vds"),
        remove_originals=False,
    )

    # Convert a VDS to a MatrixTable
    vds_path = TESTS_DIR / "vds/cohort.vds"
    mt_output_path = TESTS_DIR / "mts/cohort.mt"
    convert_vds_to_mt(
        vds_path=str(vds_path),
        output_path=str(mt_output_path),
        adjust_genotypes=True,
        skip_validation=False,  # Run validation in tests
        skip_split_multi=False,
        skip_keying_by_cols=False,
        overwrite=True,
    )

    # Check if the mt_output_path / '_SUCCESS' file exists
    assert mt_output_path.joinpath("_SUCCESS").exists()

    # Compress the MatrixTable into a zip file and remove the original MatrixTable
    compress_files(
        source_dir=str(mt_output_path),
        output_zip=str(mt_output_path.with_suffix(".mt.zip")),
        remove_originals=True,
    )

    # Cleanup temporary files or directories
    # Both VDS and MatrixTable(s) contain numerous files, so we remove the entire directory after compressing
    if vds_path.exists():
        shutil.rmtree(vds_path)
    if mt_output_path.exists():
        shutil.rmtree(mt_output_path)


@pytest.mark.hail
@pytest.mark.order4
def test_convert_mt_to_cvcf():
    """
    Test the convert_mt_to_multi_sample_vcf function by converting a MatrixTable to a multi-sample (cohort) VCF.

    This test uses the session-scoped hail_session fixture for Hail initialization.
    """
    # Decompress the MatrixTable zip file
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(TESTS_DIR / "mts/cohort.mt"),
        remove_originals=False,
    )

    # Convert a MatrixTable to a multi-sample VCF
    mt_path = TESTS_DIR / "mts/cohort.mt"
    vcf_output_path = TESTS_DIR / "vcf/cohort.vcf.bgz"
    convert_mt_to_multi_sample_vcf(
        mt_path=str(mt_path),
        vcf_path=str(vcf_output_path),
        filter_adj_genotypes=True,
        min_ac=1,
        split_multi=True,
    )

    # Check if the vcf_output_path file exists
    assert vcf_output_path.exists()

    # Cleanup temporary files or directories
    if mt_path.exists():
        shutil.rmtree(mt_path)
