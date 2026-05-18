"""GTEx eQTL drift probe: portal HEAD + version-meta scrape.

Runs OFFLINE - the gtexportal endpoint is stubbed via requests_mock.
"""

from __future__ import annotations

import hashlib

import requests
import requests_mock
import pytest

from hvantk.core.plugin_api import DriftProbeError
from hvantk.skills.gtex_eqtl.drift_probe import (
    fetch_fingerprint,
    GTEX_PORTAL_DOWNLOAD_URL,
)


def test_fetch_fingerprint_shape():
    body = (
        '<html><head><meta name="gtex-version" content="cl361-2-gd31acd02">'
        '<meta name="gtex-build-date" content="3/13/2026, 2:21:21 PM">'
        "</head><body>...</body></html>"
    )
    with requests_mock.Mocker() as m:
        m.get(
            GTEX_PORTAL_DOWNLOAD_URL,
            text=body,
            headers={
                "Last-Modified": "Fri, 13 Mar 2026 18:22:51 GMT",
                "ETag": '"69b455fb-2424"',
            },
        )
        fp = fetch_fingerprint()

    expected_blob = (
        b'cl361-2-gd31acd02|3/13/2026, 2:21:21 PM|"69b455fb-2424"'
    )
    expected_checksum = hashlib.sha256(expected_blob).hexdigest()

    assert fp["probe_version"] == 1
    assert fp["source_version"] == "cl361-2-gd31acd02"
    assert fp["headers"]["gtex-portal-qtl-downloads"] == [
        "gtex-version",
        "gtex-build-date",
        "ETag",
        "Last-Modified",
    ]
    assert fp["checksums"]["gtex-portal-qtl-downloads"] == expected_checksum
    assert fp["extras"]["gtex_version"] == "cl361-2-gd31acd02"
    assert fp["extras"]["gtex_build_date"] == "3/13/2026, 2:21:21 PM"
    assert fp["extras"]["etag"] == '"69b455fb-2424"'
    assert "fetched_at" in fp


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.get(GTEX_PORTAL_DOWNLOAD_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()


def test_fetch_fingerprint_raises_when_no_version_meta():
    with requests_mock.Mocker() as m:
        m.get(GTEX_PORTAL_DOWNLOAD_URL, text="<html>no meta here</html>")
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
