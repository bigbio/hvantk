"""Expression Atlas drift probe: experiment-index JSON fingerprint.

Runs OFFLINE - the gxa experiments endpoint is stubbed via requests_mock.
"""

from __future__ import annotations

import hashlib
import json

import requests
import requests_mock
import pytest

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.expression_atlas.drift_probe import (
    fetch_fingerprint,
    EXPRESSION_ATLAS_INDEX_URL,
)


def test_fetch_fingerprint_shape():
    payload = {
        "experiments": [
            {"experimentAccession": "E-MTAB-4045", "lastUpdate": "04-02-2016"},
            {"experimentAccession": "E-GEOD-61857", "lastUpdate": "15-05-2019"},
        ]
    }
    with requests_mock.Mocker() as m:
        m.get(EXPRESSION_ATLAS_INDEX_URL, json=payload)
        fp = fetch_fingerprint()

    expected_projection = sorted(
        (
            ("E-MTAB-4045", "04-02-2016"),
            ("E-GEOD-61857", "15-05-2019"),
        )
    )
    expected_checksum = hashlib.sha256(
        json.dumps(expected_projection, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    assert fp["probe_version"] == 1
    assert fp["source_version"] == "2 experiments"
    assert fp["headers"]["experiments.json"] == ["experimentAccession", "lastUpdate"]
    assert fp["checksums"]["experiments.json"] == expected_checksum
    assert "fetched_at" in fp


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.get(EXPRESSION_ATLAS_INDEX_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()


def test_fetch_fingerprint_raises_on_invalid_json():
    with requests_mock.Mocker() as m:
        m.get(EXPRESSION_ATLAS_INDEX_URL, text="not json")
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
