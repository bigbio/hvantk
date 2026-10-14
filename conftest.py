"""Repo-root pytest configuration.

Holds only options that must be registered no matter which subtree is under test.

``--regenerate-snapshots`` previously lived in ``hvantk/tests/conftest.py``, which
pytest loads only when the selected paths sit under ``hvantk/tests/`` -- conftest
discovery walks a test file's *ancestor* directories, and ``hvantk/skills/`` is a
sibling of ``hvantk/tests/``, not a descendant. Running a plugin's suite on its own,
exactly as each ``plugin.yaml`` documents::

    pytest hvantk/skills/hgnc/tests -m hail

therefore never registered the option, and every snapshot test died at setup with
``ValueError: no option named '--regenerate-snapshots'`` -- the fixture asks for an
option nothing had added. The rootdir is this directory (pytest.ini lives here), so a
conftest here is an ancestor of both trees and is always loaded.

Keep this file minimal: fixtures belong in ``hvantk/tests/conftest.py`` (re-exported by
each plugin's ``tests/conftest.py``), not here.

``--require-round-trips`` and the hook below live here for the same reason: a dataset's
round trip may run scoped to its own plugin directory or as part of the whole suite from
here, and the enforcement needs to see it either way (issue #432). The companion
collection-time check lives in ``hvantk/tests/test_plugin_contract_artifacts.py``; this
hook covers what that one structurally cannot -- a ``pytest.skip()`` (or
``pytest.importorskip``) raised from inside a round trip's own body, which only shows up
once the test actually runs.
"""

from __future__ import annotations

import pytest

_ROUND_TRIP_VIOLATIONS: pytest.StashKey[list] = pytest.StashKey()


def pytest_addoption(parser):
    parser.addoption(
        "--regenerate-snapshots",
        action="store_true",
        default=False,
        help="Regenerate skill snapshot files in place instead of asserting against them.",
    )
    parser.addoption(
        "--require-round-trips",
        action="store_true",
        default=False,
        help=(
            "Fail the run if any @pytest.mark.round_trip(...) test ends skipped or "
            "xfailed -- including a pytest.skip() call inside the test body, which no "
            "collection-time check can see. Ignored under --regenerate-snapshots, since "
            "every round trip deliberately self-skips there after writing snapshots."
        ),
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Record a round_trip-marked item whose setup/call phase ends skipped or xfailed."""
    outcome = yield
    config = item.config
    if not config.getoption("--require-round-trips") or config.getoption(
        "--regenerate-snapshots"
    ):
        return
    report = outcome.get_result()
    if report.when not in ("setup", "call"):
        return
    round_trip_ids = [mark.args[0] for mark in item.iter_markers("round_trip")]
    if not round_trip_ids:
        return
    wasxfail = getattr(report, "wasxfail", None)
    if not (report.skipped or wasxfail):
        return
    reason = wasxfail or _skip_reason(report)
    violations = config.stash.setdefault(_ROUND_TRIP_VIOLATIONS, [])
    violations.append(
        f"{item.nodeid} ({', '.join(round_trip_ids)}) [{report.when}]: {reason}"
    )


def _skip_reason(report) -> str:
    longrepr = report.longrepr
    if isinstance(longrepr, tuple) and len(longrepr) == 3:
        return str(longrepr[2])
    return str(longrepr)


def pytest_sessionfinish(session):
    """Fail the session if any round_trip test was recorded as skipped/xfailed."""
    if session.config.stash.get(_ROUND_TRIP_VIOLATIONS, None):
        session.exitstatus = 1


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    violations = config.stash.get(_ROUND_TRIP_VIOLATIONS, None)
    if not violations:
        return
    terminalreporter.section("round_trip tests that did not really run")
    for violation in violations:
        terminalreporter.line(violation, red=True)
    terminalreporter.line(
        f"{len(violations)} round_trip-marked test(s) ended skipped or xfailed; "
        "--require-round-trips fails the run on that (pass --regenerate-snapshots to "
        "exempt a snapshot-regeneration run).",
        red=True,
    )
