"""HGNC drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live HGNC endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.skills.hgnc.drift_probe import fetch_fingerprint
from hvantk.core.constants import HGNC_DOWNLOAD_URL


def test_fetch_fingerprint_shape():
    fake_first_line = "hgnc_id\tsymbol\tname\tstatus\tlocus_type\n"
    with requests_mock.Mocker() as m:
        m.head(HGNC_DOWNLOAD_URL, headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"})
        m.get(HGNC_DOWNLOAD_URL, text=fake_first_line)
        fp = fetch_fingerprint()
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"]["hgnc_complete_set.txt"] == [
        "hgnc_id", "symbol", "name", "status", "locus_type"
    ]
    assert "checksums" in fp
    assert "fetched_at" in fp
