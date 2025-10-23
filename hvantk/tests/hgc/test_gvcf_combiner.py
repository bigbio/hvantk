import shutil
from pathlib import Path

import pytest

from hvantk.hgc.combiners import combine_gvcfs
from hvantk.hgc.file_utils import compress_files

TESTS_DIR = Path(__file__).parent.parent / "testdata" / "hgc_data"


@pytest.mark.hail
@pytest.mark.order2
def test_combine_gvcfs():
    """
    Test the combine_gvcfs function by combining GVCFs into a VDS.

    This test uses the session-scoped hail_session fixture for Hail initialization.
    """
    # Combine GVCFs into a VDS
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

    # Check if the vds_output_path / 'reference_data/_SUCCESS' file exists
    assert vds_output_path.joinpath("reference_data/_SUCCESS").exists()
    # Check if the vds_output_path / 'variant_data/_SUCCESS' file exists
    assert vds_output_path.joinpath("variant_data/_SUCCESS").exists()

    # Compress the VDS into a zip file and remove the original VDS
    compress_files(
        source_dir=str(vds_output_path),
        output_zip=str(vds_output_path.with_suffix(".vds.zip")),
        remove_originals=True,
    )

    # cleanup temporary files or directories
    if vds_output_path.exists():
        shutil.rmtree(vds_output_path)

    # remove temporary dir
    if tmp_path.exists():
        shutil.rmtree(tmp_path)
