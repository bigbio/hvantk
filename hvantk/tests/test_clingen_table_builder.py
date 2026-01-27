"""
Unit tests for ClinGen Gene-Disease Validity CLI commands.

These tests mock the builder functions to test CLI argument parsing
without requiring Hail.
"""

from unittest.mock import patch

from click.testing import CliRunner

from hvantk.commands.make_table_cli import mktable_group


def test_mktable_clingen_gene_disease_default_options():
    """Test CLI with default options."""
    runner = CliRunner()
    with patch(
        "hvantk.commands.make_table_cli._create_clingen_gene_disease_tb"
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "clingen-gene-disease",
                "--raw-input",
                "/path/to/clingen.csv",
                "--output-ht",
                "/out/clingen.ht",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/clingen.csv",
            output_path="/out/clingen.ht",
            key_by="gene_disease",
            min_classification=None,
            fields=None,
            overwrite=False,
            export_tsv=False,
        )
        assert "ClinGen Gene-Disease table created" in result.output


def test_mktable_clingen_gene_disease_key_by_gene():
    """Test CLI with gene-level keying."""
    runner = CliRunner()
    with patch(
        "hvantk.commands.make_table_cli._create_clingen_gene_disease_tb"
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "clingen-gene-disease",
                "--raw-input",
                "/path/to/clingen.csv",
                "--output-ht",
                "/out/clingen.ht",
                "--key-by",
                "gene",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once()
        kwargs = mock_create.call_args.kwargs
        assert kwargs["key_by"] == "gene"


def test_mktable_clingen_gene_disease_min_classification():
    """Test CLI with minimum classification filter."""
    runner = CliRunner()
    with patch(
        "hvantk.commands.make_table_cli._create_clingen_gene_disease_tb"
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "clingen-gene-disease",
                "--raw-input",
                "/path/to/clingen.csv",
                "--output-ht",
                "/out/clingen.ht",
                "--min-classification",
                "Moderate",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once()
        kwargs = mock_create.call_args.kwargs
        assert kwargs["min_classification"] == "Moderate"


def test_mktable_clingen_gene_disease_all_options():
    """Test CLI with all options specified."""
    runner = CliRunner()
    with patch(
        "hvantk.commands.make_table_cli._create_clingen_gene_disease_tb"
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "clingen-gene-disease",
                "--raw-input",
                "/path/to/clingen.csv",
                "--output-ht",
                "/out/clingen.ht",
                "--key-by",
                "gene",
                "--min-classification",
                "Strong",
                "--fields",
                "hgnc_id, gene_symbol, classification",
                "--overwrite",
                "--export-tsv",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/clingen.csv",
            output_path="/out/clingen.ht",
            key_by="gene",
            min_classification="Strong",
            fields=["hgnc_id", "gene_symbol", "classification"],
            overwrite=True,
            export_tsv=True,
        )


def test_mktable_clingen_gene_disease_invalid_classification():
    """Test CLI rejects invalid classification level."""
    runner = CliRunner()
    result = runner.invoke(
        mktable_group,
        [
            "clingen-gene-disease",
            "--raw-input",
            "/path/to/clingen.csv",
            "--output-ht",
            "/out/clingen.ht",
            "--min-classification",
            "Invalid",
        ],
    )
    assert result.exit_code != 0
    assert "Invalid value" in result.output


def test_mktable_clingen_gene_disease_invalid_key_by():
    """Test CLI rejects invalid key_by value."""
    runner = CliRunner()
    result = runner.invoke(
        mktable_group,
        [
            "clingen-gene-disease",
            "--raw-input",
            "/path/to/clingen.csv",
            "--output-ht",
            "/out/clingen.ht",
            "--key-by",
            "invalid",
        ],
    )
    assert result.exit_code != 0
    assert "Invalid value" in result.output
