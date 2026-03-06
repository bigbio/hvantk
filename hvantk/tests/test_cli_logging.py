"""
Tests for centralized logging configuration in the hvantk CLI.
"""

import logging

from click.testing import CliRunner

from hvantk.hvantk import cli, setup_logging


class TestSetupLogging:
    """Tests for the setup_logging() function."""

    def setup_method(self):
        """Reset root logger handlers before each test."""
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)
            h.close()
        root.setLevel(logging.WARNING)

    def test_default_verbosity_sets_warning(self):
        """Default verbosity (0) should set WARNING level."""
        setup_logging(verbosity=0)
        assert logging.getLogger().level == logging.WARNING

    def test_verbosity_one_sets_info(self):
        """Verbosity 1 (-v) should set INFO level."""
        setup_logging(verbosity=1)
        assert logging.getLogger().level == logging.INFO

    def test_verbosity_two_sets_debug(self):
        """Verbosity 2 (-vv) should set DEBUG level."""
        setup_logging(verbosity=2)
        assert logging.getLogger().level == logging.DEBUG

    def test_verbosity_greater_than_two_sets_debug(self):
        """Verbosity > 2 should still set DEBUG level."""
        setup_logging(verbosity=5)
        assert logging.getLogger().level == logging.DEBUG

    def test_no_log_file_creates_single_handler(self):
        """Without log_file, only a StreamHandler should be added."""
        setup_logging(verbosity=0, log_file=None)
        root = logging.getLogger()
        stream_handlers = [
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]
        assert len(stream_handlers) >= 1

    def test_log_file_creates_file_handler(self, tmp_path):
        """With log_file, a FileHandler should be added."""
        log_path = str(tmp_path / "test.log")
        setup_logging(verbosity=1, log_file=log_path)
        root = logging.getLogger()
        file_handlers = [
            h for h in root.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) >= 1

    def test_log_file_receives_messages(self, tmp_path):
        """Log messages should be written to the specified log file."""
        log_path = tmp_path / "test.log"
        setup_logging(verbosity=1, log_file=str(log_path))
        test_logger = logging.getLogger("test_cli_logging.file_test")
        test_logger.info("test message for file")
        # Flush handlers
        for h in logging.getLogger().handlers:
            h.flush()
        content = log_path.read_text()
        assert "test message for file" in content


class TestCLIVerboseOption:
    """Tests for the -v/--verbose CLI option."""

    def setup_method(self):
        """Reset root logger handlers before each test."""
        root = logging.getLogger()
        for h in root.handlers[:]:
            root.removeHandler(h)
            h.close()
        root.setLevel(logging.WARNING)

    def test_help_shows_verbose_option(self):
        """The --help output should include -v/--verbose."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "--verbose" in result.output
        assert "-v" in result.output

    def test_help_shows_log_file_option(self):
        """The --help output should include --log-file."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "--log-file" in result.output

    def test_default_no_flags_shows_help(self):
        """Running without flags or subcommand shows usage (Click group behavior)."""
        runner = CliRunner()
        result = runner.invoke(cli, [])
        # Click groups without invoke_without_command show usage info
        assert "Usage:" in result.output

    def test_verbose_flag_accepted_with_help(self):
        """The -v flag before --help should be accepted without error."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-v", "--help"])
        assert result.exit_code == 0
        assert "--verbose" in result.output

    def test_double_verbose_flag_accepted_with_help(self):
        """The -vv flag before --help should be accepted without error."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-vv", "--help"])
        assert result.exit_code == 0

    def test_log_file_option_accepted_with_help(self, tmp_path):
        """The --log-file option before --help should be accepted without error."""
        runner = CliRunner()
        log_path = str(tmp_path / "cli_test.log")
        result = runner.invoke(cli, ["--log-file", log_path, "--help"])
        assert result.exit_code == 0
