"""UniProt PTM drift probe: a real release signal, fail-closed (#271).

Offline via requests_mock. Probe version 1 recorded `source_version: null` (the search
endpoint sends no Last-Modified) and hashed only the keys of the first result -- the
shape of the API, not the data -- so a UniProt release that revised every PTM annotation
compared clean. UniProt exposes its release on every response (`X-UniProt-Release`,
`X-UniProt-Release-Date`) plus the result count of the exact query (`X-Total-Results`).
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.uniprot_ptm.drift_probe import PROBE_VERSION, fetch_fingerprint
from hvantk.skills.uniprot_ptm.shared.constants import (
    UNIPROT_API_FIELDS,
    UNIPROT_API_URL,
    UNIPROT_HUMAN_PTM_QUERY,
)
from hvantk.skills.uniprot_ptm.shared.datasets import _TSV_COLUMNS, _build_search_url

PROBE_URL = _build_search_url(UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1)

_LIVE_HEADERS = {
    "X-UniProt-Release": "2026_03",
    "X-UniProt-Release-Date": "02-September-2026",
    "X-Total-Results": "9493",
    "Content-Type": "application/json",
}

_BODY = {
    "results": [
        {
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "extraAttributes": {"uniParcId": "UPI0000160563"},
            "primaryAccession": "P38398",
            "genes": [{"geneName": {"value": "BRCA1"}}],
            "sequence": {"value": "MDLSALRVEEV", "length": 11},
            "uniProtKBCrossReferences": [{"database": "Ensembl", "id": "ENST00000357654"}],
            "features": [{"type": "Modified residue", "description": "Phosphoserine",
                          "location": {"start": {"value": 2}}}],
        }
    ]
}


def _probe(headers=_LIVE_HEADERS, body=_BODY, status=200):
    with requests_mock.Mocker() as m:
        m.get(PROBE_URL, json=body, headers=headers, status_code=status)
        fp = fetch_fingerprint()
        return fp, m


def test_fingerprint_carries_the_release_the_count_and_the_schema():
    fp, m = _probe()
    assert fp["probe_version"] == PROBE_VERSION == 2
    assert fp["source_version"] == "2026_03"
    assert fp["extras"] == {"total_results": 9493}
    assert fp["informational"] == {"release_date": "02-September-2026"}
    assert fp["headers"]["uniprot-ptm-human.tsv"] == list(_TSV_COLUMNS)
    assert fp["headers"]["uniprot_entry_keys"] == [
        "entryType", "extraAttributes", "features", "genes",
        "primaryAccession", "sequence", "uniProtKBCrossReferences",
    ]
    assert "uniprot-ptm-human.tsv" in fp["checksums"]
    assert "fetched_at" in fp
    assert [r.method for r in m.request_history] == ["GET"], "one GET; the HEAD is gone"


def test_a_new_release_moves_the_compared_surface():
    from hvantk.core.plugin.drift_runner import _compare_fingerprints

    before, _ = _probe()
    after, _ = _probe(headers={**_LIVE_HEADERS, "X-UniProt-Release": "2026_04",
                               "X-UniProt-Release-Date": "15-October-2026",
                               "X-Total-Results": "9512"})
    diff = _compare_fingerprints(before, after)
    assert diff is not None
    assert set(diff["changed"]) == {"source_version", "extras"}


def test_release_date_alone_never_signals_drift():
    from hvantk.core.plugin.drift_runner import _compare_fingerprints

    a, _ = _probe()
    b, _ = _probe(headers={**_LIVE_HEADERS, "X-UniProt-Release-Date": "03-September-2026"})
    assert _compare_fingerprints(a, b) is None


@pytest.mark.parametrize("missing", ["X-UniProt-Release", "X-Total-Results"])
def test_missing_release_or_count_header_fails_closed(missing):
    headers = {k: v for k, v in _LIVE_HEADERS.items() if k != missing}
    with pytest.raises(DriftProbeError, match=missing):
        _probe(headers=headers)


def test_non_integer_count_fails_closed():
    with pytest.raises(DriftProbeError, match="X-Total-Results"):
        _probe(headers={**_LIVE_HEADERS, "X-Total-Results": "many"})


def test_zero_results_fails_closed():
    with pytest.raises(DriftProbeError, match="zero results"):
        _probe(body={"results": []})


def test_http_error_is_a_probe_error():
    # 503 is in RETRY_STATUSES (request_with_retry backs off across 4 attempts, ~14s of
    # real sleep) -- 404 is not retried, so this stays a fast, offline check of the same
    # fail-closed path (HTTP failure -> DriftProbeError) without patching time.sleep.
    with pytest.raises(DriftProbeError, match="HTTP"):
        _probe(status=404)
