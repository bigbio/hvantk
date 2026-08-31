"""alphagenome drift probe should fingerprint the published SDK release set.

Runs OFFLINE via requests_mock. The prediction service itself stays credentialed
(issue #177 was right that there is no static artifact); what this probe reaches
is the openly published SDK release stream.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.alphagenome.drift_probe import (
    ALPHAGENOME_PYPI_URL,
    fetch_fingerprint,
)


def _payload(current="0.8.0", releases=("0.7.0", "0.8.0")):
    return {
        "info": {"name": "alphagenome", "version": current},
        "releases": {v: [] for v in releases},
    }


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.get(ALPHAGENOME_PYPI_URL, json=_payload())
        fp = fetch_fingerprint()

    assert fp["source_version"] == "0.8.0"
    assert fp["extras"]["releases_found"] == ["0.7.0", "0.8.0"]


def test_new_sdk_release_moves_the_checksum():
    with requests_mock.Mocker() as m:
        m.get(ALPHAGENOME_PYPI_URL, json=_payload())
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(
            ALPHAGENOME_PYPI_URL,
            json=_payload(current="0.9.0", releases=("0.7.0", "0.8.0", "0.9.0")),
        )
        after = fetch_fingerprint()

    assert before["checksums"] != after["checksums"]
    assert after["source_version"] == "0.9.0"


def test_missing_version_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(ALPHAGENOME_PYPI_URL, json={"info": {}, "releases": {}})
        with pytest.raises(DriftProbeError, match="no version or no releases"):
            fetch_fingerprint()


def test_non_json_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(ALPHAGENOME_PYPI_URL, text="<html>not json</html>")
        with pytest.raises(DriftProbeError, match="non-JSON"):
            fetch_fingerprint()
