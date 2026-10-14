"""Standalone pytest plugin: report each collected item's round-trip marker state.

Not part of the normal conftest chain -- it is loaded only by the per-manifest
collection check in ``test_plugin_contract_artifacts.py``, which runs each
dataset's own ``tests.command`` as a subprocess with ``--collect-only -p
hvantk.tests._round_trip_introspection --round-trip-report=<path>``.

Writing JSON from a subprocess, rather than parsing ``--collect-only`` text output
or re-deriving marker state from source, means the check sees exactly what pytest
itself decided: items already filtered by that command's own ``-k``/``-m``, and
skip/xfail state resolved through pytest's own marker evaluation (so a decorator,
a ``pytestmark`` list, an aliased marker, or a class-based test all resolve the
same way, and a ``skipif(False, ...)`` is correctly not flagged).
"""

from __future__ import annotations

import json
from pathlib import Path

from _pytest.skipping import evaluate_skip_marks, evaluate_xfail_marks


def pytest_addoption(parser):
    parser.addoption(
        "--round-trip-report",
        action="store",
        default=None,
        help="(internal) write collected items' round_trip marker state as JSON here.",
    )


def pytest_collection_finish(session):
    """Dump one row per item pytest decided to keep, after all deselection."""
    report_path = session.config.getoption("--round-trip-report")
    if not report_path:
        return
    rows = []
    for item in session.items:
        round_trip_ids = [mark.args[0] for mark in item.iter_markers("round_trip")]
        if not round_trip_ids:
            continue
        skip = evaluate_skip_marks(item)
        xfail = evaluate_xfail_marks(item)
        rows.append(
            {
                "nodeid": item.nodeid,
                "round_trip": round_trip_ids,
                "skip_or_xfail": skip is not None or xfail is not None,
                "hail": next(item.iter_markers("hail"), None) is not None,
            }
        )
    Path(report_path).write_text(json.dumps(rows))
