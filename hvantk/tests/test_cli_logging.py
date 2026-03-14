"""Tests for centralized logging configuration in the hvantk CLI."""

import logging

import pytest
from click.testing import CliRunner

from hvantk.hvantk import cli, setup_logging


class TestSetupLogging:
    """Tests for the setup_logging() function."""

    def setup_method(self):
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)
            h.close()
        root.setLevel(logging.WARNING)

    @pytest.mark.parametrize(
        "verbosity,expected_level",
        [
            (0, logging.WARNING),
            (1, logging.INFO),
            (2, logging.DEBUG),
            (5, logging.DEBUG),
        ],
        ids=["default-WARNING", "v-INFO", "vv-DEBUG", "vvvvv-DEBUG"],
    )
    def test_verbosity_levels(self, verbosity, expected_level):
        """Verbosity maps to correct log level."""
        setup_logging(verbosity=verbosity)
        assert logging.getLogger().level == expected_level

    def test_log_file_receives_messages(self, tmp_path):
        """Log messages are written to the specified log file."""
        log_path = tmp_path / "test.log"
        setup_logging(verbosity=1, log_file=str(log_path))
        test_logger = logging.getLogger("test_cli_logging.file_test")
        test_logger.info("test message for file")
        for h in logging.getLogger().handlers:
            h.flush()
        content = log_path.read_text()
        assert "test message for file" in content


class TestCLIVerboseOption:
    """Tests for the -v/--verbose and --log-file CLI options."""

    def test_help_shows_verbose_and_log_file_options(self):
        """The --help output includes -v/--verbose and --log-file."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "--verbose" in result.output
        assert "-v" in result.output
        assert "--log-file" in result.output

    def test_verbose_flags_accepted(self):
        """The -v and -vv flags are accepted without error."""
        runner = CliRunner()
        for flags in ["-v", "-vv"]:
            result = runner.invoke(cli, [flags, "--help"])
            assert result.exit_code == 0
