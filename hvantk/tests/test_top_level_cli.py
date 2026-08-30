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


def test_lazy_command_registry_matches_real_commands():
    """``_LAZY_COMMANDS`` must stay in sync with the commands it stands in for.

    The registry duplicates each subcommand's name and short help so that
    listing commands (and ``hvantk --help``) does not have to import them --
    that indirection is what keeps startup at ~0.1 s instead of ~10 s. The
    duplication is only safe if something notices when it drifts, which is
    what this test is. It imports all 18 subcommand modules, so it is the one
    place that pays the full import cost on purpose.
    """
    import importlib

    from hvantk.hvantk import _LAZY_COMMANDS

    mismatches = []
    for name, (module, attr, short_help) in sorted(_LAZY_COMMANDS.items()):
        command = getattr(importlib.import_module(module), attr)
        if command.name != name:
            mismatches.append(
                f"{module}.{attr}: registry key {name!r} != command name "
                f"{command.name!r}"
            )
        actual = command.get_short_help_str(limit=200)
        if actual != short_help:
            mismatches.append(
                f"{module}.{attr}: registry help {short_help!r} != actual {actual!r}"
            )

    assert not mismatches, "stale _LAZY_COMMANDS entries:\n  " + "\n  ".join(mismatches)
