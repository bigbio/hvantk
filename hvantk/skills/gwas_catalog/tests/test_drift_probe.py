"""GWAS Catalog drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD so CI
never hits the live EBI FTP endpoint.
"""

from __future__ import annotations

import requests
import requests_mock
import pytest

from hvantk.core.plugin_api import DriftProbeError
from hvantk.skills.gwas_catalog.drift_probe import (
    fetch_fingerprint,
    GWAS_CATALOG_FULL_URL,
)


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(
            GWAS_CATALOG_FULL_URL,
            headers={
                "Last-Modified": "Tue, 28 Apr 2026 11:02:27 GMT",
                "ETag": '"382ae5a-650832ab65939"',
                "Content-Length": "58895962",
            },
        )
        fp = fetch_fingerprint()
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Tue, 28 Apr 2026 11:02:27 GMT"
    assert fp["headers"]["gwas-catalog-associations-full.zip"] == [
        "ETag",
        "Content-Length",
        "Last-Modified",
    ]
    assert "checksums" in fp
    assert fp["checksums"]["gwas-catalog-associations-full.zip"]
    assert "fetched_at" in fp


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.head(GWAS_CATALOG_FULL_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
