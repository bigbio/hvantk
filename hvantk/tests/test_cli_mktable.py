from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from hvantk.commands.make_table_cli import mktable_group


def test_mktable_clinvar_cli():
    runner = CliRunner()
    with patch("hvantk.commands.make_table_cli._create_clinvar_tb") as mock_create:
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
    with patch("hvantk.commands.make_table_cli._create_interactome_tb") as mock_create:
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
    with patch("hvantk.commands.make_table_cli._create_gevir_tb") as mock_create:
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
        "hvantk.commands.make_table_cli._create_gnomad_constraint_gene_metrics_tb"
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
    with patch("hvantk.commands.make_table_cli._create_ensembl_gene_tb") as mock_create:
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
    with patch("hvantk.commands.make_table_cli._create_dbnsfp_tb") as mock_create:
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
