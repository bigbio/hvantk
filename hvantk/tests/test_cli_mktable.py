import warnings
from unittest.mock import patch, MagicMock

import pytest
from click.testing import CliRunner

from hvantk.tools.build.make_table_cli import mktable_group


def test_mktable_clinvar_cli():
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_clinvar_tb") as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "clinvar",
                "--raw-input",
                "/path/to/clinvar.vcf.bgz",
                "--output-ht",
                "/out/clinvar.ht",
                "--ref-genome",
                "GRCh38",
                "--overwrite",
                "--export-tsv",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/clinvar.vcf.bgz",
            output_path="/out/clinvar.ht",
            overwrite=True,
            export_tsv=True,
            reference_genome="GRCh38",
        )
        assert "ClinVar table created" in result.output


def test_mktable_interactome_cli():
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_interactome_tb") as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "interactome",
                "--raw-input",
                "/path/to/interactome.bed.bgz",
                "--output-ht",
                "/out/interactome.ht",
                "--ref-genome",
                "GRCh38",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/interactome.bed.bgz",
            output_path="/out/interactome.ht",
            overwrite=False,
            export_tsv=False,
            reference_genome="GRCh38",
        )
        assert "Interactome table created" in result.output


def test_mktable_gevir_cli_with_fields():
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_gevir_tb") as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "gevir",
                "--raw-input",
                "/path/to/gevir.tsv.bgz",
                "--output-ht",
                "/out/gevir.ht",
                "--fields",
                "oe_syn_upper, oe_mis_upper",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/gevir.tsv.bgz",
            output_path="/out/gevir.ht",
            fields=["oe_syn_upper", "oe_mis_upper"],
            overwrite=False,
            export_tsv=False,
        )
        assert "GEVIR table created" in result.output


def test_mktable_gnomad_metrics_cli():
    runner = CliRunner()
    with patch(
        "hvantk.tools.build.make_table_cli._create_gnomad_constraint_gene_metrics_tb"
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "gnomad-metrics",
                "--raw-input",
                "/path/to/gnomad.tsv.bgz",
                "--output-ht",
                "/out/gnomad.ht",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/gnomad.tsv.bgz",
            output_path="/out/gnomad.ht",
            fields=None,
            overwrite=False,
            export_tsv=False,
        )
        assert "gnomAD metrics table created" in result.output


def test_mktable_ensembl_gene_cli_no_canonical():
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_ensembl_gene_tb") as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "ensembl-gene",
                "--raw-input",
                "/path/to/biomart.tsv.bgz",
                "--output-ht",
                "/out/ensembl.ht",
                "--no-canonical",
                "--export-tsv",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            input_path="/path/to/biomart.tsv.bgz",
            output_path="/out/ensembl.ht",
            fields=None,
            canonical=False,
            overwrite=False,
            export_tsv=True,
        )
        assert "Ensembl gene table created" in result.output


def test_mktable_dbnsfp_cli_invokes_builder():
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_dbnsfp_tb") as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "dbnsfp",
                "--raw-input",
                "/data/dbNSFP4.bgz",
                "--output-ht",
                "/out/dbnsfp.ht",
                "--ref-genome",
                "GRCh38",
                "--group-prefixes",
                "gnomAD,ExAC",
                "--overwrite",
            ],
        )
        assert result.exit_code == 0
        mock_create.assert_called_once()
        kwargs = mock_create.call_args.kwargs
        assert kwargs["input_path"] == "/data/dbNSFP4.bgz"
        assert kwargs["output_path"] == "/out/dbnsfp.ht"
        assert kwargs["reference_genome"].upper() == "GRCH38"
        assert kwargs["overwrite"] is True
        assert kwargs["group_prefixes"] == ["gnomAD", "ExAC"]
        # parse_transcript_scores defaults to True
        assert kwargs["parse_transcript_scores"] is True
        assert "dbNSFP table created" in result.output


def test_mktable_clingen_gene_disease_default_options():
    """Test ClinGen CLI with default options."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.build.make_table_cli._create_clingen_gene_disease_tb"
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


def test_mktable_clingen_gene_disease_all_options():
    """Test ClinGen CLI with all options specified."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.build.make_table_cli._create_clingen_gene_disease_tb"
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


@pytest.mark.parametrize(
    "option,value",
    [
        ("--min-classification", "Invalid"),
        ("--key-by", "invalid"),
    ],
)
def test_mktable_clingen_gene_disease_rejects_invalid_values(option, value):
    """Test ClinGen CLI rejects invalid option values."""
    runner = CliRunner()
    result = runner.invoke(
        mktable_group,
        [
            "clingen-gene-disease",
            "--raw-input",
            "/path/to/clingen.csv",
            "--output-ht",
            "/out/clingen.ht",
            option,
            value,
        ],
    )
    assert result.exit_code != 0
    assert "Invalid value" in result.output


# ---------------------------------------------------------------------------
# Deprecation warnings
# ---------------------------------------------------------------------------

def test_mktable_clinvar_emits_deprecation_warning():
    """mktable clinvar should emit a DeprecationWarning directing to reprocess."""
    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_cli._create_clinvar_tb") as mock_create:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", DeprecationWarning)
            result = runner.invoke(
                mktable_group,
                [
                    "clinvar",
                    "--raw-input", "/tmp/fake.vcf.bgz",
                    "--output-ht", "/tmp/fake_out.ht",
                ],
            )
    deprecation_warnings = [
        w for w in caught
        if issubclass(w.category, DeprecationWarning)
        and "hvantk reprocess" in str(w.message)
    ]
    assert len(deprecation_warnings) >= 1, (
        "Expected a DeprecationWarning mentioning 'hvantk reprocess' "
        f"but got: {[str(w.message) for w in caught]}"
    )


def test_mktable_group_help_mentions_deprecated():
    """The mktable group --help output should mention deprecation."""
    runner = CliRunner()
    result = runner.invoke(mktable_group, ["--help"])
    assert result.exit_code == 0
    assert "Deprecated" in result.output or "deprecated" in result.output


def test_mkmatrix_ucsc_emits_deprecation_warning():
    """mkmatrix ucsc should emit a DeprecationWarning directing to reprocess."""
    from hvantk.tools.build.make_matrix_cli import mkmatrix_group
    import os

    runner = CliRunner()
    # Create minimal stub files so click's exists=True check passes
    with runner.isolated_filesystem():
        open("expr.tsv", "w").close()
        open("meta.tsv", "w").close()
        with patch("hvantk.tools.build.make_matrix_cli._build_ucsc_ad") as mock_build:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", DeprecationWarning)
                result = runner.invoke(
                    mkmatrix_group,
                    [
                        "ucsc",
                        "-e", "expr.tsv",
                        "-m", "meta.tsv",
                        "-o", "out.h5ad",
                    ],
                )
    deprecation_warnings = [
        w for w in caught
        if issubclass(w.category, DeprecationWarning)
        and "hvantk reprocess" in str(w.message)
    ]
    assert len(deprecation_warnings) >= 1, (
        "Expected a DeprecationWarning mentioning 'hvantk reprocess' "
        f"but got: {[str(w.message) for w in caught]}"
    )
