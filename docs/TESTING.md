# Testing Guide

This project’s test suite is organized by behavior and runtime characteristics using pytest markers, so day-to-day development runs are fast and focused, while full validation remains available.

## Test categories (pytest markers)

- fast (default): unit-style tests without external systems (no marker needed)
- slow: heavier CPU/IO or long-running
- network: real network/FTP/HTTP activity
- hail: requires Hail/Spark and sizable data
- integration: multi-component flows, end-to-end
- llm: LLM-related tests (mocked or real providers)

Markers are declared in pytest.ini and the default addopts skips slow, network, hail, integration, and llm by default.

## Common commands

- Fast local run (default):
  pytest

- Explicit fast subset:
  pytest -m "not slow and not network and not hail and not integration and not llm"

- Run Hail-only tests:
  pytest -m hail

- Run network-only tests:
  pytest -m network

- Full suite:
  pytest -m "slow or network or hail or integration or llm"

- Parallelize (optional, if pytest-xdist installed):
  pytest -n auto -m hail

## What changed in this refactor

- Introduced markers and default filtering via pytest.ini; constrained discovery to hvantk/tests for faster collection.
- Marked heavy tests accordingly (hail, network, slow, integration, llm).
- Refactored network tests to use mocks where practical (e.g., requests.get, ftplib.FTP) and added a fast downloader test.
- Moved heavy Hail imports inside tests/fixtures to reduce import-time overhead.
- Set a headless matplotlib backend in tests/conftest.py.

## Suggested directory structure (optional)

Markers are sufficient, but if you prefer folders:

- hvantk/tests/fast/
- hvantk/tests/hail/
- hvantk/tests/network/
- hvantk/tests/integration/
- hvantk/tests/llm/

You can auto-apply markers by folder using pytest_collection_modifyitems in tests/conftest.py if you adopt this layout.

## CI recommendations

- PR checks: run the default fast subset.
- Nightly/cron: run the full suite with -n auto.
- For flaky/network tests, add retries or timeouts as needed.

## Notes on Hail tests

- Prefer tiny fixtures or small testdata excerpts to keep runtime low.
- Consider a session-scoped Hail init fixture if you add many Hail tests; reuse it across tests to avoid repeated initialization overhead.

