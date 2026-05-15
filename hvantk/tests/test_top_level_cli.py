"""Smoke test: top-level CLI exposes plugins and drift subcommands."""

from click.testing import CliRunner

from hvantk.hvantk import cli


def test_plugins_subcommand_is_attached():
    runner = CliRunner()
    result = runner.invoke(cli, ["plugins", "--help"])
    assert result.exit_code == 0
    assert "list" in result.output


def test_drift_subcommand_is_attached():
    runner = CliRunner()
    result = runner.invoke(cli, ["drift", "--help"])
    assert result.exit_code == 0
    assert "--all" in result.output
