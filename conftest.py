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
"""

from __future__ import annotations


def pytest_addoption(parser):
    parser.addoption(
        "--regenerate-snapshots",
        action="store_true",
        default=False,
        help="Regenerate skill snapshot files in place instead of asserting against them.",
    )
