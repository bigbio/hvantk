"""ClinVar drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD so CI
never hits the live NCBI FTP endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.skills.clinvar.drift_probe import fetch_fingerprint, _URL


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(
            _URL,
            headers={
                "Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT",
                "Content-Length": "12345",
            },
        )
        fp = fetch_fingerprint()
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"]["clinvar.vcf.gz"]["content_length"] == "12345"
    assert "checksums" in fp
    assert "fetched_at" in fp
