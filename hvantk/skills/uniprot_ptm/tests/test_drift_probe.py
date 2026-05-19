"""UniProt PTM drift probe should return a fingerprint with the documented shape.

This test runs OFFLINE - it uses requests_mock to stub the HTTP HEAD/GET so
CI never hits the live UniProt REST API. A live integration test would need
network access; we do not run it in this suite.
"""

from __future__ import annotations

import requests_mock

from hvantk.core.ptm_constants import (
    UNIPROT_API_FIELDS,
    UNIPROT_API_URL,
    UNIPROT_HUMAN_PTM_QUERY,
)
from hvantk.skills.uniprot_ptm.drift_probe import fetch_fingerprint
from hvantk.skills.uniprot_ptm.shared.datasets import _build_search_url


def test_fetch_fingerprint_shape():
    probe_url = _build_search_url(
        UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1
    )
    fake_body = {
        "results": [
            {
                "primaryAccession": "P38398",
                "genes": [{"geneName": {"value": "BRCA1"}}],
                "sequence": {"value": "MDLSALRVEEV", "length": 11},
                "uniProtKBCrossReferences": [
                    {"database": "Ensembl", "id": "ENST00000357654"}
                ],
                "features": [
                    {
                        "type": "Modified residue",
                        "description": "Phosphoserine",
                        "location": {"start": {"value": 2}},
                    }
                ],
            }
        ]
    }
    with requests_mock.Mocker() as m:
        m.head(
            probe_url,
            headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"},
        )
        m.get(probe_url, json=fake_body)
        fp = fetch_fingerprint()

    expected_filename = "uniprot-ptm-human.tsv"
    assert fp["probe_version"] == 1
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"][expected_filename] == [
        "accession",
        "gene_symbol",
        "position",
        "description",
        "amino_acid",
        "ensembl_xrefs",
        "sequence_length",
    ]
    assert fp["headers"]["uniprot_entry_keys"] == [
        "features",
        "genes",
        "primaryAccession",
        "sequence",
        "uniProtKBCrossReferences",
    ]
    assert expected_filename in fp["checksums"]
    assert "fetched_at" in fp
