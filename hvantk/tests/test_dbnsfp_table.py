import shutil
from pathlib import Path

import pytest

from hvantk.core.builders.table import create_dbnsfp_tb

# Mark as Hail-dependent and slow
pytestmark = [pytest.mark.hail, pytest.mark.slow]

TEST_DIR = Path(__file__).parent / "testdata"
TMP_DIR = Path(__file__).parent / "tmp_dbnsfp"


@pytest.fixture(autouse=True)
def setup_teardown():
    TMP_DIR.mkdir(exist_ok=True, parents=True)
    yield
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)


def test_create_dbnsfp_tb_builds_and_keys_variants():
    input_path = TEST_DIR / "raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz"
    output_path = TMP_DIR / "dbnsfp_example.ht"

    tb = create_dbnsfp_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        overwrite=True,
        export_tsv=False,
        min_partitions=10,
        reference_genome="GRCh38",
        parse_transcript_scores=True,
    )

    # sanity: row count equals number of input records
    assert tb.count() == 4999

    # output exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()
