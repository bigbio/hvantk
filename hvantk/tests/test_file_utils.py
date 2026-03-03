import gzip
import shutil
import zipfile
from unittest import mock

import pytest

from hvantk.data.file_utils import (
    compress_files,
    convert_gz_to_bgz,
    decompress_files,
    detect_compression,
    is_bgzf,
    is_gzipped,
    resolve_compression,
    _get_conversion_backend,
    _make_bgzf_block,
    _convert_with_python,
)


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


# ---------------------------------------------------------------------------
# Fixtures for compression tests
# ---------------------------------------------------------------------------

_SAMPLE_CONTENT = b"Hello, this is test data for BGZF conversion.\n" * 100


@pytest.fixture
def plain_gz_file(tmp_path):
    """Create a standard gzip file."""
    gz_path = tmp_path / "sample.tsv.gz"
    with gzip.open(str(gz_path), "wb") as f:
        f.write(_SAMPLE_CONTENT)
    return str(gz_path)


@pytest.fixture
def bgzf_file(tmp_path):
    """Create a BGZF file using pure Python block writer."""
    bgz_path = tmp_path / "sample.tsv.bgz"
    # Write a proper BGZF file with our helper
    with open(str(bgz_path), "wb") as fout:
        fout.write(_make_bgzf_block(_SAMPLE_CONTENT))
        fout.write(_make_bgzf_block(b""))  # EOF block
    return str(bgz_path)


@pytest.fixture
def plain_text_file(tmp_path):
    """Create an uncompressed text file."""
    txt_path = tmp_path / "sample.tsv"
    txt_path.write_bytes(_SAMPLE_CONTENT)
    return str(txt_path)


@pytest.fixture
def empty_file(tmp_path):
    """Create an empty file."""
    empty_path = tmp_path / "empty.gz"
    empty_path.write_bytes(b"")
    return str(empty_path)


# ---------------------------------------------------------------------------
# Detection tests
# ---------------------------------------------------------------------------


class TestIsGzipped:
    def test_plain_gz(self, plain_gz_file):
        assert is_gzipped(plain_gz_file) is True

    def test_bgzf(self, bgzf_file):
        assert is_gzipped(bgzf_file) is True  # BGZF is also gzip

    def test_plain_text(self, plain_text_file):
        assert is_gzipped(plain_text_file) is False

    def test_empty_file(self, empty_file):
        assert is_gzipped(empty_file) is False

    def test_nonexistent_file(self):
        with pytest.raises(FileNotFoundError):
            is_gzipped("/nonexistent/file.gz")


class TestIsBgzf:
    def test_plain_gz(self, plain_gz_file):
        assert is_bgzf(plain_gz_file) is False

    def test_bgzf(self, bgzf_file):
        assert is_bgzf(bgzf_file) is True

    def test_plain_text(self, plain_text_file):
        assert is_bgzf(plain_text_file) is False

    def test_empty_file(self, empty_file):
        assert is_bgzf(empty_file) is False

    def test_nonexistent_file(self):
        with pytest.raises(FileNotFoundError):
            is_bgzf("/nonexistent/file.gz")


class TestDetectCompression:
    def test_plain_gz(self, plain_gz_file):
        assert detect_compression(plain_gz_file) == "gzip"

    def test_bgzf(self, bgzf_file):
        assert detect_compression(bgzf_file) == "bgzf"

    def test_plain_text(self, plain_text_file):
        assert detect_compression(plain_text_file) == "none"

    def test_empty_file(self, empty_file):
        assert detect_compression(empty_file) == "none"


# ---------------------------------------------------------------------------
# Backend detection tests
# ---------------------------------------------------------------------------


class TestGetConversionBackend:
    def test_prefers_bgzip_when_available(self):
        with mock.patch("shutil.which", return_value="/usr/bin/bgzip"):
            assert _get_conversion_backend() == "bgzip"

    def test_falls_back_to_pysam(self):
        with mock.patch("shutil.which", return_value=None):
            # pysam should be importable in our test env
            backend = _get_conversion_backend()
            assert backend in ("pysam", "python")

    def test_falls_back_to_python(self):
        with mock.patch("shutil.which", return_value=None), mock.patch.dict(
            "sys.modules", {"pysam": None}
        ):
            assert _get_conversion_backend() == "python"


# ---------------------------------------------------------------------------
# Conversion tests
# ---------------------------------------------------------------------------


class TestConvertGzToBgz:
    def test_conversion_produces_valid_bgzf(self, plain_gz_file, tmp_path):
        output = str(tmp_path / "output.bgz")
        result = convert_gz_to_bgz(plain_gz_file, output_path=output)
        assert result == output
        assert is_bgzf(output)

    def test_content_roundtrip(self, plain_gz_file, tmp_path):
        output = str(tmp_path / "output.bgz")
        convert_gz_to_bgz(plain_gz_file, output_path=output)
        # Decompress and verify content matches
        with gzip.open(output, "rb") as f:
            decompressed = f.read()
        assert decompressed == _SAMPLE_CONTENT

    def test_default_output_path(self, plain_gz_file):
        result = convert_gz_to_bgz(plain_gz_file)
        assert result.endswith(".bgz")
        assert is_bgzf(result)

    def test_nonexistent_input(self):
        with pytest.raises(FileNotFoundError):
            convert_gz_to_bgz("/nonexistent/file.gz")

    def test_pure_python_backend(self, plain_gz_file, tmp_path):
        """Directly test the pure Python fallback."""
        output = str(tmp_path / "python_output.bgz")
        _convert_with_python(plain_gz_file, output)
        assert is_bgzf(output)
        with gzip.open(output, "rb") as f:
            assert f.read() == _SAMPLE_CONTENT

    def test_pysam_backend(self, plain_gz_file, tmp_path):
        """Test pysam backend if available."""
        try:
            import pysam  # noqa: F401
        except ImportError:
            pytest.skip("pysam not installed")

        from hvantk.data.file_utils import _convert_with_pysam

        output = str(tmp_path / "pysam_output.bgz")
        _convert_with_pysam(plain_gz_file, output)
        assert is_bgzf(output)
        with gzip.open(output, "rb") as f:
            assert f.read() == _SAMPLE_CONTENT

    def test_bgzip_backend(self, plain_gz_file, tmp_path):
        """Test system bgzip backend if available."""
        if not shutil.which("bgzip"):
            pytest.skip("bgzip not on PATH")

        from hvantk.data.file_utils import _convert_with_bgzip

        output = str(tmp_path / "bgzip_output.bgz")
        _convert_with_bgzip(plain_gz_file, output, threads=2)
        assert is_bgzf(output)
        with gzip.open(output, "rb") as f:
            assert f.read() == _SAMPLE_CONTENT


# ---------------------------------------------------------------------------
# resolve_compression tests
# ---------------------------------------------------------------------------


class TestResolveCompression:
    def test_bgzf_passthrough(self, bgzf_file):
        path, force = resolve_compression(bgzf_file)
        assert path == bgzf_file
        assert force is True

    def test_plain_gz_warns_and_disables_force(self, plain_gz_file):
        path, force = resolve_compression(plain_gz_file, auto_convert=False)
        assert path == plain_gz_file
        assert force is False

    def test_plain_gz_auto_convert(self, plain_gz_file):
        path, force = resolve_compression(plain_gz_file, auto_convert=True)
        assert path != plain_gz_file
        assert path.endswith(".bgz")
        assert force is True
        assert is_bgzf(path)

    def test_uncompressed_file(self, plain_text_file):
        path, force = resolve_compression(plain_text_file)
        assert path == plain_text_file
        assert force is False
