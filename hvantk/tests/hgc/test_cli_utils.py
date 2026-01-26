"""
Tests for HGC CLI utilities.
"""

import os
from unittest.mock import patch

from hvantk.commands.hgc.utils import (
    setup_logging_for_hgc,
    expand_file_patterns,
    validate_input_files,
    validate_output_path,
    estimate_resource_requirements,
    DEFAULT_TEMP_DIR,
)


def test_setup_logging_for_hgc():
    """Test logging setup."""
    # Should not raise any errors
    setup_logging_for_hgc("INFO")
    setup_logging_for_hgc("DEBUG")
    setup_logging_for_hgc("WARNING")
    setup_logging_for_hgc("ERROR")


def test_expand_file_patterns(tmp_path):
    """Test file pattern expansion."""
    # Create test files
    (tmp_path / "file1.txt").touch()
    (tmp_path / "file2.txt").touch()
    (tmp_path / "file3.vcf").touch()

    # Test wildcard expansion
    pattern = str(tmp_path / "*.txt")
    expanded = expand_file_patterns([pattern])
    assert len(expanded) == 2
    assert all(f.endswith(".txt") for f in expanded)

    # Test multiple patterns
    patterns = [str(tmp_path / "*.txt"), str(tmp_path / "*.vcf")]
    expanded = expand_file_patterns(patterns)
    assert len(expanded) == 3

    # Test no matches
    pattern = str(tmp_path / "*.nonexistent")
    expanded = expand_file_patterns([pattern])
    assert len(expanded) == 0


@patch("hvantk.commands.hgc.utils.check_path_exists_and_readable")
def test_validate_input_files_gvcf(mock_check):
    """Test input file validation for GVCF files."""
    mock_check.return_value = True

    # Valid files
    is_valid, errors = validate_input_files(["/path/to/file.vcf"], "gvcf")
    assert is_valid
    assert len(errors) == 0
    mock_check.assert_called_once_with("/path/to/file.vcf")

    # Invalid files
    mock_check.side_effect = Exception("File not found")
    is_valid, errors = validate_input_files(["/path/to/missing.vcf"], "gvcf")
    assert not is_valid
    assert len(errors) == 1
    assert "File not found" in errors[0]


@patch("hvantk.commands.hgc.utils.validate_vds_paths")
def test_validate_input_files_vds(mock_validate):
    """Test input file validation for VDS files."""
    mock_validate.return_value = True

    # Valid VDS
    is_valid, errors = validate_input_files(["/path/to/file.vds"], "vds")
    assert is_valid
    assert len(errors) == 0
    mock_validate.assert_called_once_with(["/path/to/file.vds"])

    # Invalid VDS
    mock_validate.side_effect = Exception("Invalid VDS")
    is_valid, errors = validate_input_files(["/path/to/bad.vds"], "vds")
    assert not is_valid
    assert len(errors) == 1
    assert "Invalid VDS" in errors[0]


def test_validate_output_path(tmp_path):
    """Test output path validation."""
    # Valid path with existing directory
    output_path = str(tmp_path / "output.txt")
    assert validate_output_path(output_path)

    # Path with non-existent parent directory (no create)
    output_path = str(tmp_path / "nonexistent" / "output.txt")
    assert not validate_output_path(output_path, create_dirs=False)

    # Path with non-existent parent directory (create)
    output_path = str(tmp_path / "newdir" / "output.txt")
    assert validate_output_path(output_path, create_dirs=True)
    assert os.path.exists(str(tmp_path / "newdir"))


def test_estimate_resource_requirements(tmp_path):
    """Test resource estimation."""
    # Create test files
    file1 = tmp_path / "file1.txt"
    file1.write_text("x" * 1024 * 1024)  # 1 MB
    file2 = tmp_path / "file2.txt"
    file2.write_text("y" * 1024 * 1024 * 5)  # 5 MB

    # Estimate for files
    results = estimate_resource_requirements([str(file1), str(file2)])
    assert "memory" in results
    assert "partitions" in results
    assert "estimated_runtime_minutes" in results
    assert "total_size_gb" in results
    assert results["total_size_gb"] > 0

    # Estimate for directory
    results = estimate_resource_requirements([str(tmp_path)])
    assert results["total_size_gb"] > 0


def test_estimate_resource_requirements_empty():
    """Test resource estimation with empty list."""
    results = estimate_resource_requirements([])
    assert results["total_size_gb"] == 0
    assert results["memory"] == "4g"  # Minimum
    assert results["partitions"] == 100  # Minimum


def test_estimate_resource_requirements_nonexistent():
    """Test resource estimation with nonexistent paths."""
    # Should not raise, just skip nonexistent paths
    results = estimate_resource_requirements(["/nonexistent/path"])
    assert results["total_size_gb"] == 0


def test_default_temp_dir():
    """Test DEFAULT_TEMP_DIR constant."""
    assert DEFAULT_TEMP_DIR is not None
    assert isinstance(DEFAULT_TEMP_DIR, str)
