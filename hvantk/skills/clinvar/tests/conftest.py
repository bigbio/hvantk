"""Local conftest for the ClinVar plugin tests.

Re-exports the session-level fixtures (``hail_session``, ``regenerate_snapshots``)
and the auto-marker hook from ``hvantk/tests/conftest.py``. Pytest's conftest
discovery only walks ancestor directories, and ``hvantk/skills/`` is a sibling
of ``hvantk/tests/`` -- without this shim, hail-marked tests in this folder
would not receive the shared session fixture.

``pytest_addoption`` is deliberately NOT re-exported here: it lives in the
top-level test conftest and pytest only needs it registered once across the
session. Re-exporting it would raise ``ValueError: option names already added``
when both conftests load (e.g., when running the full suite).
"""

from hvantk.tests.conftest import (  # noqa: F401
    hail_session,
    pytest_collection_modifyitems,
    regenerate_snapshots,
)
