"""GenCC drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live GenCC endpoint. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.skills.gencc.shared.constants import GENCC_BASE_URL, GENCC_FILE_PREFIX
from hvantk.skills.gencc.drift_probe import fetch_fingerprint


def test_fetch_fingerprint_shape():
    fake_body = (
        "sgc_id\tgene_curie\tgene_symbol\tdisease_curie\tdisease_title\t"
        "disease_original_curie\tdisease_original_title\tclassification_title\t"
        "moi_title\tsubmitter_title\tsubmitted_as_date\t"
        "submitted_as_public_report_url\tsubmitted_as_pmids\n"
        "SGC-001\tHGNC:1100\tBRCA1\tMONDO:0005144\tbreast-ovarian cancer 1\t"
        "MONDO:0005144\tbreast-ovarian cancer 1\tDefinitive\tAutosomal dominant\t"
        "ClinGen\t2025-06-15\thttps://example/\t12345678\n"
    )
    with requests_mock.Mocker() as m:
        m.head(
            GENCC_BASE_URL,
            headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"},
        )
        m.get(GENCC_BASE_URL, text=fake_body)
        fp = fetch_fingerprint()

    expected_filename = f"{GENCC_FILE_PREFIX}.tsv"
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"][expected_filename] == [
        "sgc_id",
        "gene_curie",
        "gene_symbol",
        "disease_curie",
        "disease_title",
        "disease_original_curie",
        "disease_original_title",
        "classification_title",
        "moi_title",
        "submitter_title",
        "submitted_as_date",
        "submitted_as_public_report_url",
        "submitted_as_pmids",
    ]
    assert expected_filename in fp["checksums"]
    assert "fetched_at" in fp
