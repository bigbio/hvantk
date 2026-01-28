from pathlib import Path

from click.testing import CliRunner
import pandas as pd

from hvantk.commands.enrichex_cli import enrichex_group


def _write_enrichment(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes"],
            "n_overlap": [12, 8],
            "odds_ratio": [2.5, 1.7],
            "p_value": [1e-6, 5e-4],
            "p_adjusted": [5e-6, 1e-3],
            "significant": [True, True],
        }
    )
    path = tmp_path / "enrichment.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def _write_burden(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes"],
            "odds_ratio": [1.6, 1.2],
            "ci_lower": [1.2, 0.9],
            "ci_upper": [2.1, 1.4],
            "p_value": [0.001, 0.08],
            "p_adjusted": [0.002, 0.12],
            "significant": [True, False],
        }
    )
    path = tmp_path / "burden.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_cli_plot_dotplot(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        enrichex_group,
        [
            "plot",
            "dotplot",
            "-i",
            str(_write_enrichment(tmp_path)),
            "-o",
            str(tmp_path / "dot.png"),
            "--top-n",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "dot.png").exists()


def test_cli_report_command(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        enrichex_group,
        [
            "report",
            "--overlap-results",
            str(_write_enrichment(tmp_path)),
            "--burden-results",
            str(_write_burden(tmp_path)),
            "-o",
            str(tmp_path / "report.html"),
            "--title",
            "CLI Report",
        ],
    )
    assert result.exit_code == 0, result.output
    contents = (tmp_path / "report.html").read_text()
    assert "CLI Report" in contents
