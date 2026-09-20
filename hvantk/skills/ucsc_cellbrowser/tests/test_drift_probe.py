"""UCSC Cell Browser drift probe: catalog + upstream index fingerprint.

Runs OFFLINE - the upstream request is stubbed via requests_mock; the local
catalog is read from the package resources.

Two things this file deliberately does, both of which it previously did not
(issue #341, item 2):

1. **The expected index URL is pinned as a literal here**, not imported from the
   module under test. Importing ``UCSC_DATASET_INDEX_URL`` and feeding it to the
   mock made the URL assertion circular: a probe pointing at
   ``https://cells.ucsc.edu/datasets.json`` (plural - a real mistake made while
   re-authoring this probe) registered its own wrong URL with the mocker and
   passed green. The literal below is the one upstream actually serves; if the
   module constant moves, that is a deliberate change and this test should fail
   until the literal is updated too.

2. **Both HEAD and GET are stubbed.** The committed probe fingerprints the index
   with a single ``HEAD`` (ETag/Content-Length/Last-Modified), but checksumming
   the fetched body with a ``GET`` is an equally reasonable design - and under a
   HEAD-only mock it died with ``requests_mock.exceptions.NoMockAddress`` rather
   than on its merits. The probe *contract* (``SKILL.md`` s 9) is "fingerprint the
   index"; the verb is an implementation detail, so the test allows either.
"""

from __future__ import annotations

import json

import requests
import requests_mock
import pytest

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.ucsc_cellbrowser.drift_probe import (
    fetch_fingerprint,
    UCSC_DATASET_INDEX_URL,
)

# Pinned independently of the module under test - see the note above. Note the
# SINGULAR "dataset.json"; "datasets.json" does not exist upstream.
EXPECTED_INDEX_URL = "https://cells.ucsc.edu/dataset.json"

_INDEX_HEADERS = {
    "Last-Modified": "Sat, 16 May 2026 01:18:32 GMT",
    "ETag": '"21878-651e51b95b435"',
    "Content-Length": "137336",
}

# Minimal stand-in for the upstream index body, so a probe that GETs and hashes
# the payload has something deterministic to hash.
_INDEX_BODY = json.dumps(
    {"datasets": [{"name": "cortex-dev", "shortLabel": "Cortex development"}]}
)


def test_probe_targets_the_documented_index_url():
    """The probe must point at the real upstream index, not a plausible-looking one."""
    assert UCSC_DATASET_INDEX_URL == EXPECTED_INDEX_URL


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        # Register against the pinned literal, NOT the imported constant, so a
        # probe aimed elsewhere raises NoMockAddress instead of passing.
        m.head(EXPECTED_INDEX_URL, headers=_INDEX_HEADERS)
        m.get(EXPECTED_INDEX_URL, headers=_INDEX_HEADERS, text=_INDEX_BODY)
        fp = fetch_fingerprint()
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Sat, 16 May 2026 01:18:32 GMT"
    assert "cells_ucsc_datasets.json" in fp["checksums"]
    assert "dataset.json" in fp["checksums"]
    # The local-catalog hash is deterministic per ship; just assert it is non-empty.
    assert fp["checksums"]["cells_ucsc_datasets.json"]
    assert fp["checksums"]["dataset.json"]
    assert "fetched_at" in fp


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.head(EXPECTED_INDEX_URL, exc=requests.ConnectionError("boom"))
        m.get(EXPECTED_INDEX_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
