import zipfile
import pytest
from hvantk.data.file_utils import compress_files, decompress_files


@pytest.fixture
def setup_test_environment(tmp_path):
    test_dir = tmp_path / "test_dir"
    test_zip = tmp_path / "test.zip"
    extracted_dir = tmp_path / "extracted_dir"
    test_dir.mkdir()
    (test_dir / "file1.txt").write_text("This is a test file.")
    (test_dir / "file2.txt").write_text("This is another test file.")
    return test_dir, test_zip, extracted_dir


def test_compress_files(setup_test_environment):
    test_dir, test_zip, _ = setup_test_environment
    compress_files(str(test_dir), str(test_zip))
    assert test_zip.exists()
    with zipfile.ZipFile(test_zip, "r") as zipf:
        assert len(zipf.namelist()) == 2
        assert "file1.txt" in zipf.namelist()
        assert "file2.txt" in zipf.namelist()


def test_compress_files_with_removal(setup_test_environment):
    test_dir, test_zip, _ = setup_test_environment
    compress_files(str(test_dir), str(test_zip), remove_originals=True)
    assert test_zip.exists()
    assert not test_dir.exists()


def test_decompress_files(setup_test_environment):
    test_dir, test_zip, extracted_dir = setup_test_environment
    compress_files(str(test_dir), str(test_zip))
    decompress_files(str(test_zip), str(extracted_dir))
    assert extracted_dir.exists()
    assert (extracted_dir / "file1.txt").exists()
    assert (extracted_dir / "file2.txt").exists()


def test_decompress_files_with_removal(setup_test_environment):
    test_dir, test_zip, extracted_dir = setup_test_environment
    compress_files(str(test_dir), str(test_zip))
    decompress_files(str(test_zip), str(extracted_dir), remove_originals=True)
    assert extracted_dir.exists()
    assert not test_zip.exists()
