"""ClinGen drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live ClinGen endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.core.constants import CLINGEN_BASE_URL, CLINGEN_FILE_PREFIX
from hvantk.skills.clingen.drift_probe import fetch_fingerprint


def test_fetch_fingerprint_shape():
    fake_body = (
        "CLINGEN GENE VALIDITY CURATIONS\n"
        "FILE CREATED: 2026-01-15\n"
        "GENOME BUILD: GRCh38\n"
        "FORMAT: CSV\n"
        "For more information visit: https://search.clinicalgenome.org/\n"
        "----\n"
        '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)",'
        '"MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
        '"BRCA1","HGNC:1100","breast cancer","MONDO:0005144","AD","SOP8",'
        '"Definitive","https://example/","2025-06-15","Hereditary Cancer GCEP"\n'
    )
    with requests_mock.Mocker() as m:
        m.head(
            CLINGEN_BASE_URL,
            headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"},
        )
        m.get(CLINGEN_BASE_URL, text=fake_body)
        fp = fetch_fingerprint()

    expected_filename = f"{CLINGEN_FILE_PREFIX}.csv"
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"][expected_filename] == [
        "GENE SYMBOL",
        "GENE ID (HGNC)",
        "DISEASE LABEL",
        "DISEASE ID (MONDO)",
        "MOI",
        "SOP",
        "CLASSIFICATION",
        "ONLINE REPORT",
        "CLASSIFICATION DATE",
        "GCEP",
    ]
    assert expected_filename in fp["checksums"]
    assert "fetched_at" in fp
