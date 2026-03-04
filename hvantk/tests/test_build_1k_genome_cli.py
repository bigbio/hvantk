"""
Unit tests for the build-1k-genome CLI and file-discovery helpers.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from hvantk.commands.build_1k_genome_cli import build_1k_genome_cmd
from hvantk.tables.genome_builders import _extract_chrom_token, discover_vcf_files


def _make_vcf_dir(tmp_path: Path, chroms: list[str], include_index: bool = True) -> Path:
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    for chrom in chroms:
        vcf = vcf_dir / f"1kg_{chrom}.vcf.gz"
        vcf.touch()
        if include_index:
            (vcf_dir / f"1kg_{chrom}.vcf.gz.tbi").touch()
    return vcf_dir


# ---------------------------------------------------------------------------
# discover_vcf_files
# ---------------------------------------------------------------------------


def test_discover_raises_if_no_vcfs(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="No \\*.vcf.gz files found"):
        discover_vcf_files(str(empty))


def test_discover_raises_if_missing_tbi(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1"], include_index=False)
    with pytest.raises(FileNotFoundError, match="missing a .tbi tabix index"):
        discover_vcf_files(str(vcf_dir))


def test_discover_sorted_by_chromosome(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chrX", "chr22", "chr1", "chr10", "chr2"])
    files = discover_vcf_files(str(vcf_dir))
    tokens = [_extract_chrom_token(f) for f in files]
    assert tokens == ["chr1", "chr2", "chr10", "chr22", "chrX"]


def test_discover_chromosome_filter(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1", "chr2", "chrX"])
    files = discover_vcf_files(str(vcf_dir), chromosomes=["chr1", "chrX"])
    tokens = {_extract_chrom_token(f) for f in files}
    assert tokens == {"chr1", "chrX"}


def test_discover_filter_raises_no_match(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1", "chr2"])
    with pytest.raises(ValueError, match="No VCF files matched"):
        discover_vcf_files(str(vcf_dir), chromosomes=["chrY"])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_basic_invocation(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1"])
    mock_mt = MagicMock()
    mock_mt.count_rows.return_value = 1000
    mock_mt.count_cols.return_value = 2504

    runner = CliRunner()
    with patch(
        "hvantk.commands.build_1k_genome_cli._build_1k_genome_mt",
        return_value=mock_mt,
    ) as mock_build:
        result = runner.invoke(
            build_1k_genome_cmd,
            ["--input-vcfs", str(vcf_dir), "--output-mt", str(tmp_path / "out.mt")],
        )

    assert result.exit_code == 0, result.output
    mock_build.assert_called_once_with(
        input_vcfs=str(vcf_dir),
        output_mt=str(tmp_path / "out.mt"),
        sample_annotations=None,
        sample_annotations_delimiter=None,
        reference_genome="GRCh38",
        chromosomes=None,
        overwrite=False,
    )
    assert "1000 Genomes MatrixTable created" in result.output


def test_cli_all_options(tmp_path):
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1", "chr2"])
    annot_file = tmp_path / "annotations.tsv"
    annot_file.write_text("sample_id\tpopulation\nNA12878\tCEU\n")

    mock_mt = MagicMock()
    mock_mt.count_rows.return_value = 50
    mock_mt.count_cols.return_value = 10

    runner = CliRunner()
    with patch(
        "hvantk.commands.build_1k_genome_cli._build_1k_genome_mt",
        return_value=mock_mt,
    ) as mock_build:
        result = runner.invoke(
            build_1k_genome_cmd,
            [
                "--input-vcfs", str(vcf_dir),
                "--output-mt", str(tmp_path / "out.mt"),
                "--sample-annotations", str(annot_file),
                "--reference-genome", "GRCh37",
                "--chromosomes", "chr1, chr2",
                "--overwrite",
            ],
        )

    assert result.exit_code == 0, result.output
    kw = mock_build.call_args.kwargs
    assert kw["reference_genome"] == "GRCh37"
    assert kw["chromosomes"] == ["chr1", "chr2"]
    assert kw["overwrite"] is True
    assert kw["sample_annotations"] == str(annot_file)
    assert kw["sample_annotations_delimiter"] is None


def test_cli_sample_annotations_delimiter(tmp_path):
    """The --sample-annotations-delimiter option is passed through."""
    vcf_dir = _make_vcf_dir(tmp_path, ["chr1"])
    annot_file = tmp_path / "annotations.ped"
    annot_file.write_text("sample_id population\nNA12878 CEU\n")

    mock_mt = MagicMock()
    mock_mt.count_rows.return_value = 100
    mock_mt.count_cols.return_value = 5

    runner = CliRunner()
    with patch(
        "hvantk.commands.build_1k_genome_cli._build_1k_genome_mt",
        return_value=mock_mt,
    ) as mock_build:
        result = runner.invoke(
            build_1k_genome_cmd,
            [
                "--input-vcfs", str(vcf_dir),
                "--output-mt", str(tmp_path / "out.mt"),
                "--sample-annotations", str(annot_file),
                "--sample-annotations-delimiter", " ",
            ],
        )

    assert result.exit_code == 0, result.output
    kw = mock_build.call_args.kwargs
    assert kw["sample_annotations"] == str(annot_file)
    assert kw["sample_annotations_delimiter"] == " "
