from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from hvantk.commands.make_table_cli import mktable_group


def test_mktable_eqtl_cli_default():
    runner = CliRunner()
    mock_ht = MagicMock()
    with patch("hvantk.commands.make_table_cli._create_eqtl_tb", return_value=mock_ht):
        result = runner.invoke(
            mktable_group,
            [
                "eqtl",
                "--raw-input",
                "/data/gtex_v11/signif_pairs/",
                "--output-ht",
                "/out/eqtl.ht",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "eQTL table created" in result.output


def test_mktable_eqtl_cli_with_options():
    runner = CliRunner()
    mock_ht = MagicMock()
    with patch(
        "hvantk.commands.make_table_cli._create_eqtl_tb", return_value=mock_ht
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "eqtl",
                "--raw-input",
                "/data/gtex_v8/",
                "--output-ht",
                "/out/eqtl_v8.ht",
                "--source",
                "gtex_v8",
                "--tissue",
                "Liver",
                "--p-threshold",
                "0",
                "--ref-genome",
                "GRCh38",
                "--overwrite",
                "--export-tsv",
            ],
        )
        assert result.exit_code == 0, result.output
        mock_create.assert_called_once_with(
            input_path="/data/gtex_v8/",
            output_path="/out/eqtl_v8.ht",
            reference_genome="GRCh38",
            source="gtex_v8",
            tissue="Liver",
            p_threshold=0.0,
            overwrite=True,
            export_tsv=True,
        )


def test_mktable_pqtl_cli_default():
    runner = CliRunner()
    mock_ht = MagicMock()
    with patch("hvantk.commands.make_table_cli._create_pqtl_tb", return_value=mock_ht):
        result = runner.invoke(
            mktable_group,
            [
                "pqtl",
                "--raw-input",
                "/data/fang_pqtl/",
                "--output-ht",
                "/out/pqtl.ht",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "pQTL table created" in result.output


def test_mktable_pqtl_cli_with_gene_map():
    runner = CliRunner()
    mock_ht = MagicMock()
    with patch(
        "hvantk.commands.make_table_cli._create_pqtl_tb", return_value=mock_ht
    ) as mock_create:
        result = runner.invoke(
            mktable_group,
            [
                "pqtl",
                "--raw-input",
                "/data/fang_pqtl/Liver.txt.gz",
                "--output-ht",
                "/out/pqtl_liver.ht",
                "--tissue",
                "Liver",
                "--gene-map-ht",
                "/data/ensembl_gene.ht",
                "--overwrite",
            ],
        )
        assert result.exit_code == 0, result.output
        mock_create.assert_called_once_with(
            input_path="/data/fang_pqtl/Liver.txt.gz",
            output_path="/out/pqtl_liver.ht",
            reference_genome="GRCh38",
            source="gtex_fang",
            tissue="Liver",
            gene_map_ht="/data/ensembl_gene.ht",
            p_threshold=None,
            overwrite=True,
            export_tsv=False,
        )
