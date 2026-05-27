import sys
from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from hvantk.tools.qtl.qtlcascade_cli import qtlcascade_group

_MOCK_HAIL_CONTEXT = MagicMock()


def test_cascade_cmd():
    runner = CliRunner()
    mock_ht = MagicMock()
    with patch.dict(
        sys.modules,
        {"hail": MagicMock(), "hvantk.core.hail_context": _MOCK_HAIL_CONTEXT},
    ), patch(
        "hvantk.algorithms.qtlcascade.cascade.build_cascade", return_value=mock_ht
    ) as mock_build:
        result = runner.invoke(
            qtlcascade_group,
            [
                "cascade",
                "--eqtl-ht",
                "/data/eqtl.ht",
                "--pqtl-ht",
                "/data/pqtl.ht",
                "-o",
                "/out/cascade.ht",
                "--tissue",
                "Liver",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Cascade table written" in result.output
        mock_build.assert_called_once()


def test_coloc_cmd(tmp_path):
    runner = CliRunner()
    genes_file = tmp_path / "genes.txt"
    genes_file.write_text("ENSG00000000003\nENSG00000000005\n")

    import pandas as pd

    mock_df = pd.DataFrame(
        {
            "gene_id": ["ENSG00000000003"],
            "tissue": ["Liver"],
            "H0": [0.01],
            "H1": [0.02],
            "H2": [0.02],
            "H3": [0.05],
            "H4": [0.90],
            "n_variants": [50],
        }
    )

    with patch.dict(
        sys.modules,
        {"hail": MagicMock(), "hvantk.core.hail_context": _MOCK_HAIL_CONTEXT},
    ), patch(
        "hvantk.algorithms.qtlcascade.coloc.run_coloc_per_gene", return_value=mock_df
    ):
        result = runner.invoke(
            qtlcascade_group,
            [
                "coloc",
                "--eqtl-allpairs",
                "/data/eqtl_allpairs.ht",
                "--pqtl-allpairs",
                "/data/pqtl_allpairs.ht",
                "--cascade-genes",
                str(genes_file),
                "--tissue",
                "Liver",
                "-o",
                str(tmp_path / "coloc.tsv"),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Coloc results written" in result.output


def test_run_cmd_dry_run():
    runner = CliRunner()
    result = runner.invoke(
        qtlcascade_group,
        [
            "run",
            "--eqtl-ht",
            "/data/eqtl.ht",
            "--pqtl-ht",
            "/data/pqtl.ht",
            "-o",
            "/tmp/cascade_out",
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "EXECUTION PLAN" in result.output


def test_run_cmd_single_tissue():
    runner = CliRunner()
    mock_result = MagicMock()
    mock_result.n_cascade_genes = 42
    mock_result.class_counts = {"eqtl_mediated": 30, "discordant": 12}

    with patch("hvantk.algorithms.qtlcascade.pipeline.CascadePipeline") as MockPipeline:
        MockPipeline.return_value.run.return_value = mock_result
        result = runner.invoke(
            qtlcascade_group,
            [
                "run",
                "--eqtl-ht",
                "/data/eqtl.ht",
                "--pqtl-ht",
                "/data/pqtl.ht",
                "-o",
                "/tmp/cascade_out",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "42" in result.output


def test_report_cmd(tmp_path):
    runner = CliRunner()
    with patch("hvantk.algorithms.qtlcascade.report.generate_report") as mock_report:
        result = runner.invoke(
            qtlcascade_group,
            [
                "report",
                "-o",
                str(tmp_path / "report.html"),
                "--title",
                "Test Report",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Report written" in result.output
        mock_report.assert_called_once()
