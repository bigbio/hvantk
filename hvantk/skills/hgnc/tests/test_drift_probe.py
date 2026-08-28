"""HGNC drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live HGNC endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.hgnc.drift_probe import fetch_fingerprint
from hvantk.skills.hgnc.shared.constants import HGNC_DOWNLOAD_URL

_FIRST_LINE = "hgnc_id\tsymbol\tname\tstatus\tlocus_type\n"


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(
            HGNC_DOWNLOAD_URL,
            headers={
                "Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT",
                "Content-Length": "52428800",
            },
        )
        m.get(HGNC_DOWNLOAD_URL, text=_FIRST_LINE)
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 2
    assert fp["headers"]["hgnc_complete_set.txt"] == [
        "hgnc_id", "symbol", "name", "status", "locus_type"
    ]
    assert "checksums" in fp
    assert "fetched_at" in fp


def test_last_modified_is_demoted_out_of_the_compared_surface():
    """Last-Modified moved on all 8 historical hgnc drift PRs while the content
    checksum never did. It is recorded for humans, but must not be a drift trigger,
    so it belongs in `informational` and `source_version` must be None."""
    with requests_mock.Mocker() as m:
        m.head(
            HGNC_DOWNLOAD_URL,
            headers={
                "Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT",
                "Content-Length": "52428800",
            },
        )
        m.get(HGNC_DOWNLOAD_URL, text=_FIRST_LINE)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert fp["informational"]["last_modified"] == "Wed, 01 Jan 2026 00:00:00 GMT"


def test_content_length_is_recorded_as_the_content_signal():
    with requests_mock.Mocker() as m:
        m.head(
            HGNC_DOWNLOAD_URL,
            headers={"Last-Modified": "x", "Content-Length": "52428800"},
        )
        m.get(HGNC_DOWNLOAD_URL, text=_FIRST_LINE)
        fp = fetch_fingerprint()

    assert fp["extras"]["content_length"] == "52428800"


def test_missing_content_length_fails_closed():
    """Recording a fingerprint with no content signal would be baked in permanently:
    the scheduled bot regenerates drifted baselines automatically, so one transient
    omission would retire content detection for good. Match ClinGen and raise."""
    with requests_mock.Mocker() as m:
        m.head(HGNC_DOWNLOAD_URL, headers={"Last-Modified": "x"})
        m.get(HGNC_DOWNLOAD_URL, text=_FIRST_LINE)
        with pytest.raises(DriftProbeError, match="Content-Length"):
            fetch_fingerprint()
