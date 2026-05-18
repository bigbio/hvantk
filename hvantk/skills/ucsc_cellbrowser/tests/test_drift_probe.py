"""UCSC Cell Browser drift probe: catalog + upstream HEAD fingerprint.

Runs OFFLINE - the upstream HEAD is stubbed via requests_mock; the local
catalog is read from the package resources.
"""

from __future__ import annotations

import requests
import requests_mock
import pytest

from hvantk.core.plugin_api import DriftProbeError
from hvantk.skills.ucsc_cellbrowser.drift_probe import (
    fetch_fingerprint,
    UCSC_DATASET_INDEX_URL,
)


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(
            UCSC_DATASET_INDEX_URL,
            headers={
                "Last-Modified": "Sat, 16 May 2026 01:18:32 GMT",
                "ETag": '"21878-651e51b95b435"',
                "Content-Length": "137336",
            },
        )
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
        m.head(UCSC_DATASET_INDEX_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
