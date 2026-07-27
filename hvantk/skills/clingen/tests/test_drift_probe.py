"""ClinGen drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live ClinGen endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.core.plugin.drift_runner import _compare_fingerprints
from hvantk.skills.clingen.shared.constants import CLINGEN_BASE_URL, CLINGEN_FILE_PREFIX
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
            headers={
                "Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT",
                "Content-Length": "1113685",
            },
        )
        m.get(CLINGEN_BASE_URL, text=fake_body)
        fp = fetch_fingerprint()

    expected_filename = f"{CLINGEN_FILE_PREFIX}.csv"
    assert fp["probe_version"] == 2
    # Last-Modified is deliberately NOT recorded: the endpoint renders the CSV
    # per request, so the server returns the request time and every run would
    # otherwise report drift. See the probe module docstring.
    assert fp["source_version"] is None
    assert fp["extras"]["content_length"] == "1113685"
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


def test_fingerprint_is_stable_across_repeated_probes():
    """Two probes of an unchanged source must compare equal.

    Regression test for the request-time ``Last-Modified``: ClinGen renders the
    export on demand, so consecutive HEADs return different timestamps for
    byte-identical content. Recording that value made ``hvantk drift`` report
    clingen:gene-disease as drifted on every run, which is what the scheduled
    drift bot turned into a daily no-op pull request.
    """
    header_row = (
        '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)",'
        '"MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
    )
    body = f"CLINGEN GENE VALIDITY CURATIONS\nFILE CREATED: 2026-01-15\n{header_row}"

    def probe_with(last_modified: str) -> dict:
        with requests_mock.Mocker() as m:
            m.head(
                CLINGEN_BASE_URL,
                headers={
                    "Last-Modified": last_modified,
                    "Content-Length": "1113685",
                },
            )
            m.get(CLINGEN_BASE_URL, text=body)
            return fetch_fingerprint()

    first = probe_with("Mon, 27 Jul 2026 17:12:51 GMT")
    second = probe_with("Mon, 27 Jul 2026 17:12:59 GMT")

    assert _compare_fingerprints(first, second) is None
