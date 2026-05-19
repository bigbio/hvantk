"""CPTAC protein-expression drift probe test.

Runs OFFLINE: cptac package version and the GitHub releases endpoint
are both stubbed. The shared helper at
``hvantk.skills.cptac.shared.drift`` does the real work.
"""

from __future__ import annotations

from unittest import mock

import pytest
import requests
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.cptac.shared import drift as shared_drift
from hvantk.skills.cptac.expression import drift_probe


def _mock_release_response(m, *, tag_name="v1.6.0", published_at="2026-04-01T00:00:00Z"):
    m.get(
        shared_drift.GITHUB_LATEST_RELEASE_URL,
        json={"tag_name": tag_name, "published_at": published_at},
    )


def test_fingerprint_shape():
    with mock.patch.object(shared_drift, "cptac_version", return_value="1.5.13"):
        with requests_mock.Mocker() as m:
            _mock_release_response(m)
            fp = drift_probe.fetch_fingerprint()

    assert fp["probe_version"] == drift_probe.PROBE_VERSION
    assert fp["source_version"] == "v1.6.0"
    assert fp["headers"] == {
        "protein_expression": [
            "installed_cptac_version",
            "latest_release_tag",
            "latest_release_published_at",
        ]
    }
    assert fp["checksums"]["protein_expression"]
    assert fp["extras"]["installed_cptac_version"] == "1.5.13"
    assert fp["extras"]["latest_release_tag"] == "v1.6.0"
    assert "fetched_at" in fp


def test_fingerprint_raises_when_package_missing():
    with mock.patch.object(shared_drift, "cptac_version", return_value=None):
        with pytest.raises(DriftProbeError):
            drift_probe.fetch_fingerprint()


def test_fingerprint_raises_on_http_failure():
    with mock.patch.object(shared_drift, "cptac_version", return_value="1.5.13"):
        with requests_mock.Mocker() as m:
            m.get(
                shared_drift.GITHUB_LATEST_RELEASE_URL,
                exc=requests.ConnectionError("boom"),
            )
            with pytest.raises(DriftProbeError):
                drift_probe.fetch_fingerprint()
