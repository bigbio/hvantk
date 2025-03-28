import os
from hvantk.utils.file_utils import download_file
from hvantk.utils.constants import (UCSC_CELL_BROWSER_BASE_URL,
                                    EXPRESSION_MATRIX_FILE_NAME,
                                    METADATA_FILE_NAME)

dataset = "adultPancreas"
url_expression_matrix = f"{UCSC_CELL_BROWSER_BASE_URL}/{dataset}/{EXPRESSION_MATRIX_FILE_NAME}"
url_metadata = f"{UCSC_CELL_BROWSER_BASE_URL}/{dataset}/{METADATA_FILE_NAME}"

def test_download_expression_matrix(tmp_path):
    """
    Test download expression matrix file from UCSC Cell Browser

    :return: None
    """
    out_dir = str(tmp_path / "data")
    download_file(url_expression_matrix, out_dir, EXPRESSION_MATRIX_FILE_NAME)
    assert os.path.exists(os.path.join(out_dir, EXPRESSION_MATRIX_FILE_NAME))


def test_download_metadata():
    """
    Test download metadata file from UCSC Cell Browser

    :return: None
    """
    download_file(url_metadata, "data", METADATA_FILE_NAME)
    assert os.path.exists(f"data/{METADATA_FILE_NAME}")

