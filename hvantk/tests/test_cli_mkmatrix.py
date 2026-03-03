from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from hvantk.commands.make_matrix_cli import mkmatrix_group


def test_mkmatrix_ucsc_short_options(tmp_path: Path):
    expr = tmp_path / "expr.tsv.bgz"
    meta = tmp_path / "meta.tsv"
    expr.write_text("gene\tcell1\nG1\t1\n")
    meta.write_text("Cell\ttype\ncell1\ta\n")
    out_mt = tmp_path / "out.mt"

    runner = CliRunner()
    with patch("hvantk.commands.make_matrix_cli._build_ucsc_mt") as mock_build:
        result = runner.invoke(
            mkmatrix_group,
            [
                "ucsc",
                "-e",
                str(expr),
                "-m",
                str(meta),
                "-o",
                str(out_mt),
                "-g",
                "symbol",
                "--ix",
                "1",
                "-d",
                ",",
                "-p",
                "12",
                "-w",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_build.assert_called_once()
    kwargs = mock_build.call_args.kwargs
    assert kwargs["expression_matrix_path"] == str(expr)
    assert kwargs["metadata_path"] == str(meta)
    assert kwargs["output_mt"] == str(out_mt)
    assert kwargs["gene_column"] == "symbol"
    assert kwargs["metadata_index_col"] == 1
    assert kwargs["delimiter"] == ","
    assert kwargs["min_partitions"] == 12
    assert kwargs["overwrite"] is True
    assert "MatrixTable created at" in result.output


def test_mkmatrix_expression_atlas_short_options(tmp_path: Path):
    expr = tmp_path / "atlas.tsv"
    sdrf = tmp_path / "atlas.sdrf.tsv"
    expr.write_text("Gene ID\tsample_1\nENSG1\t1\n")
    sdrf.write_text("source_name\tcomment[data file]\nsample_1\tsample_1\n")
    out_mt = tmp_path / "atlas.mt"

    runner = CliRunner()
    with patch("hvantk.commands.make_matrix_cli._build_expression_atlas_mt") as mock_build:
        result = runner.invoke(
            mkmatrix_group,
            [
                "expression-atlas",
                "-e",
                str(expr),
                "-s",
                str(sdrf),
                "-o",
                str(out_mt),
                "-g",
                "GeneSymbol",
                "--sid",
                "sample",
                "-d",
                ",",
                "-p",
                "25",
                "-w",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_build.assert_called_once()
    kwargs = mock_build.call_args.kwargs
    assert kwargs["expression_matrix_path"] == str(expr)
    assert kwargs["sdrf_file"] == str(sdrf)
    assert kwargs["output_mt"] == str(out_mt)
    assert kwargs["gene_column"] == "GeneSymbol"
    assert kwargs["sample_id_column"] == "sample"
    assert kwargs["delimiter"] == ","
    assert kwargs["min_partitions"] == 25
    assert kwargs["overwrite"] is True


def test_mkmatrix_cptac_short_options(tmp_path: Path):
    expr = tmp_path / "cptac_expr.tsv"
    meta = tmp_path / "cptac_meta.tsv"
    expr.write_text("GeneID\tGene Name\tSampleID\tExpression\n1\tA\tS1\t10\n")
    meta.write_text("SampleID\tCondition\nS1\tCase\n")
    out_mt = tmp_path / "cptac.mt"

    runner = CliRunner()
    with patch("hvantk.commands.make_matrix_cli._build_cptac_mt") as mock_build:
        result = runner.invoke(
            mkmatrix_group,
            [
                "cptac",
                "-e",
                str(expr),
                "-m",
                str(meta),
                "-o",
                str(out_mt),
                "-g",
                "GeneIDCol",
                "-n",
                "GeneNameCol",
                "--sid",
                "SampleIDCol",
                "-x",
                "ExprCol",
                "-c",
                "CancerType,Stage",
                "-u",
                "Age,Score",
                "-w",
            ],
        )

    assert result.exit_code == 0, result.output
    mock_build.assert_called_once()
    kwargs = mock_build.call_args.kwargs
    assert kwargs["expression_path"] == str(expr)
    assert kwargs["metadata_path"] == str(meta)
    assert kwargs["output_mt"] == str(out_mt)
    assert kwargs["gene_id_col"] == "GeneIDCol"
    assert kwargs["gene_name_col"] == "GeneNameCol"
    assert kwargs["sample_id_col"] == "SampleIDCol"
    assert kwargs["expression_col"] == "ExprCol"
    assert kwargs["categorical_cols"] == ["CancerType", "Stage"]
    assert kwargs["numeric_cols"] == ["Age", "Score"]
    assert kwargs["overwrite"] is True
