"""
Consolidated Hail integration tests for HGC core operations.

Merged from: test_file_utils, test_gvcf_combiner, test_converters,
             test_mts_combiner, test_sort_matrix_table.
Tests are ordered to reflect the pipeline: combine → convert → combine MTs → sort.
"""

import random
from pathlib import Path

import pytest
import hail as hl

from hvantk.algorithms.hgc.file_utils import validate_vcfs_paths, sort_mts_cols
from hvantk.algorithms.hgc.combiners import (
    combine_gvcfs,
    combine_matrix_table_rows,
    combine_matrix_table_cols,
)
from hvantk.algorithms.hgc.converters import (
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
)
from hvantk.algorithms.hgc.constants import GVCF_EXTENSION
from hvantk.core.utils.file_utils import compress_files, decompress_files

TESTS_DIR = Path(__file__).parent.parent / "testdata" / "hgc_data"


def test_validate_vcfs_paths():
    """Test GVCF file discovery from directory."""
    directory = TESTS_DIR / "gvcfs"
    vcfs = validate_vcfs_paths(str(directory), pattern=f"*{GVCF_EXTENSION}")
    assert len(vcfs) == 3


@pytest.mark.hail
@pytest.mark.order2
def test_combine_gvcfs(tmp_path):
    """Test combining GVCFs into a VDS."""
    gvcf_dir = TESTS_DIR / "gvcfs"
    vds_output_path = tmp_path / "cohort.vds"
    combiner_tmp = tmp_path / "tmp"
    plan_path = combiner_tmp / "combiner_plan.json"
    combine_gvcfs(
        gvcf_dir=str(gvcf_dir),
        vds_output_path=str(vds_output_path),
        tmp_path=str(combiner_tmp),
        save_path=str(plan_path),
        vdses=[],
        kwargs={},
    )
    assert vds_output_path.joinpath("reference_data/_SUCCESS").exists()
    assert vds_output_path.joinpath("variant_data/_SUCCESS").exists()

    compress_files(
        source_dir=str(vds_output_path),
        output_zip=str(tmp_path / "cohort.vds.zip"),
        remove_originals=True,
    )


# No longer skipped on a missing gnomad install: annotate_adj is ported into
# hvantk.algorithms.hgc.adj, so adjusted genotypes always work.
@pytest.mark.hail
@pytest.mark.order3
def test_convert_vds_to_mt(tmp_path):
    """Test VDS → MatrixTable conversion."""
    decompress_files(
        zip_path=str(TESTS_DIR / "vds/cohort.vds.zip"),
        extract_to=str(tmp_path / "cohort.vds"),
        remove_originals=False,
    )
    vds_path = tmp_path / "cohort.vds"
    mt_output_path = tmp_path / "cohort.mt"
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


@pytest.mark.hail
@pytest.mark.order3
def test_convert_vds_to_mt_honours_n_partitions(tmp_path):
    """#207/#208: --n-partitions must actually change the written layout.

    #208's acceptance criteria require that the flag either changes the output's
    partitioning or does not exist; before this it was accepted, printed in the run plan,
    and read by no stage. The default is asserted too, so a regression that silently
    coalesced every run would fail rather than pass quietly.

    The fixture VDS carries 2,586 partitions, so 4 is a genuine reduction rather than a
    no-op -- `naive_coalesce` does nothing when the current count is already lower, which
    would make a badly-chosen target vacuously "pass". That precondition is asserted from
    the VDS's own partition count rather than by running a second, uncoalesced conversion:
    it is free metadata, where the extra conversion cost ~4 minutes of a job that runs on
    every push. The uncoalesced path is already covered by test_convert_vds_to_mt above.
    """
    import hail as hl

    decompress_files(
        zip_path=str(TESTS_DIR / "vds/cohort.vds.zip"),
        extract_to=str(tmp_path / "cohort.vds"),
        remove_originals=False,
    )
    vds_path = str(tmp_path / "cohort.vds")

    source_parts = hl.vds.read_vds(vds_path).variant_data.n_partitions()
    assert source_parts > 4, (
        f"fixture VDS has only {source_parts} partitions; coalescing to 4 would be a "
        f"no-op and the assertion below would prove nothing"
    )

    out = tmp_path / "coalesced.mt"
    convert_vds_to_mt(
        vds_path=vds_path,
        output_path=str(out),
        adjust_genotypes=False,
        skip_validation=True,
        skip_split_multi=False,
        skip_keying_by_cols=False,
        overwrite=True,
        n_partitions=4,
    )

    written_parts = hl.read_matrix_table(str(out)).n_partitions()
    assert written_parts == 4, f"expected 4 partitions, got {written_parts}"


@pytest.mark.hail
@pytest.mark.order3
def test_convert_vds_to_mt_rejects_a_nonsense_partition_count(tmp_path):
    """A bad count must fail before Hail does, with a message naming the parameter."""
    decompress_files(
        zip_path=str(TESTS_DIR / "vds/cohort.vds.zip"),
        extract_to=str(tmp_path / "cohort.vds"),
        remove_originals=False,
    )
    with pytest.raises(ValueError, match="n_partitions must be >= 1"):
        convert_vds_to_mt(
            vds_path=str(tmp_path / "cohort.vds"),
            output_path=str(tmp_path / "never.mt"),
            n_partitions=0,
        )


@pytest.mark.hail
@pytest.mark.order4
def test_convert_mt_to_cvcf(tmp_path):
    """Test MatrixTable → multi-sample VCF conversion."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(tmp_path / "cohort.mt"),
        remove_originals=False,
    )
    mt_path = tmp_path / "cohort.mt"
    vcf_output_path = tmp_path / "cohort.vcf.bgz"
    convert_mt_to_multi_sample_vcf(
        mt_path=str(mt_path),
        vcf_path=str(vcf_output_path),
        filter_adj_genotypes=True,
        min_ac=1,
        split_multi=True,
    )
    assert vcf_output_path.exists()


@pytest.mark.hail
@pytest.mark.order5
def test_combine_matrix_table_rows(tmp_path):
    """Test combining rows of two MatrixTables."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(tmp_path / "test.mt"),
        remove_originals=False,
    )
    mt1_path = tmp_path / "test.mt"
    mt_output_path = tmp_path / "combined.mt"
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


@pytest.mark.hail
@pytest.mark.order6
def test_combine_matrix_table_cols(tmp_path):
    """Test combining columns of two MatrixTables."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(tmp_path / "test.mt"),
        remove_originals=False,
    )
    mt1_path = tmp_path / "test.mt"
    mt_output_path = tmp_path / "combined.mt"
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


@pytest.mark.hail
@pytest.mark.order7
def test_sort_mts_cols(tmp_path):
    """Test sorting column order of MatrixTables for union_rows."""
    decompress_files(
        zip_path=str(TESTS_DIR / "mts/cohort.mt.zip"),
        extract_to=str(tmp_path / "cohort.mt"),
        remove_originals=False,
    )
    mt1_path = tmp_path / "cohort.mt"
    mt1 = hl.read_matrix_table(str(mt1_path))

    n = mt1.count_cols()
    if n <= 1:
        idx = list(range(n))
    else:
        identity = list(range(n))
        idx = identity.copy()
        while idx == identity:
            random.shuffle(idx)
    mt2 = mt1.choose_cols(idx)

    sorted_mts = sort_mts_cols([mt1, mt2])
    mt = hl.MatrixTable.union_rows(*sorted_mts)

    assert mt.count_rows() == mt1.count_rows() + mt2.count_rows()
    assert mt.count_cols() == mt1.count_cols() == mt2.count_cols()
