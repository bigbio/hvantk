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
from hvantk.skills.uniprot_ptm.shared.datasets import _build_search_url

PROBE_URL = _build_search_url(
    UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1
)

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
    # Pinned literally, not against `_TSV_COLUMNS`: comparing against the same
    # constant the probe reads back cannot catch a change to that constant.
    assert fp["headers"]["uniprot-ptm-human.tsv"] == [
        "accession",
        "gene_symbol",
        "position",
        "description",
        "amino_acid",
        "ensembl_xrefs",
        "sequence_length",
    ]
    assert fp["headers"]["uniprot_entry_keys"] == [
        "entryType",
        "extraAttributes",
        "features",
        "genes",
        "primaryAccession",
        "sequence",
        "uniProtKBCrossReferences",
    ]
    assert "uniprot-ptm-human.tsv" in fp["checksums"]
    assert "fetched_at" in fp
    assert [r.method for r in m.request_history] == ["GET"], "one GET; the HEAD is gone"


def test_a_new_release_moves_the_compared_surface():
    from hvantk.core.plugin.drift_runner import _compare_fingerprints

    before, _ = _probe()
    after, _ = _probe(
        headers={
            **_LIVE_HEADERS,
            "X-UniProt-Release": "2026_04",
            "X-UniProt-Release-Date": "15-October-2026",
            "X-Total-Results": "9512",
        }
    )
    diff = _compare_fingerprints(before, after)
    assert diff is not None
    assert set(diff["changed"]) == {"source_version", "extras"}


def test_release_date_alone_never_signals_drift():
    from hvantk.core.plugin.drift_runner import _compare_fingerprints

    a, _ = _probe()
    b, _ = _probe(
        headers={**_LIVE_HEADERS, "X-UniProt-Release-Date": "03-September-2026"}
    )
    assert _compare_fingerprints(a, b) is None


@pytest.mark.parametrize("missing", ["X-UniProt-Release", "X-Total-Results"])
def test_missing_release_or_count_header_fails_closed(missing):
    headers = {k: v for k, v in _LIVE_HEADERS.items() if k != missing}
    with pytest.raises(DriftProbeError, match=missing):
        _probe(headers=headers)


def test_blank_release_header_fails_closed():
    # Present but whitespace-only -- distinct from the missing-header case above:
    # `.strip()` must still reduce it to "no version signal" rather than recording
    # a blank string as a release tag.
    with pytest.raises(DriftProbeError, match="X-UniProt-Release"):
        _probe(headers={**_LIVE_HEADERS, "X-UniProt-Release": "  "})


def test_non_json_body_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(PROBE_URL, text="<html>", headers=_LIVE_HEADERS)
        with pytest.raises(DriftProbeError, match="Non-JSON"):
            fetch_fingerprint()


def test_non_integer_count_fails_closed():
    with pytest.raises(DriftProbeError, match="X-Total-Results"):
        _probe(headers={**_LIVE_HEADERS, "X-Total-Results": "many"})


def test_zero_results_fails_closed():
    with pytest.raises(DriftProbeError, match="zero results"):
        _probe(body={"results": []})


@pytest.mark.parametrize("total", ["0", "-1"])
def test_impossible_total_results_fails_closed(total):
    """X-Total-Results must be >= the number of results actually returned.

    A bad edge/proxy response can report a count lower than what it sends back --
    here 0 or -1 while `results` still carries the one entry the size=1 probe
    requested. Without this check `int(total_raw)` succeeds either way and the
    impossible count is recorded as a routine content change, so the drift bot
    would propose it as the new baseline instead of failing closed.
    """
    with pytest.raises(DriftProbeError, match="X-Total-Results"):
        _probe(headers={**_LIVE_HEADERS, "X-Total-Results": total})


def test_http_error_is_a_probe_error():
    # 503 is in RETRY_STATUSES (request_with_retry backs off across the probe's 3
    # attempts, ~6s of real sleep) -- 404 is not retried, so this stays a fast, offline
    # check of the same fail-closed path (HTTP failure -> DriftProbeError) without
    # patching time.sleep.
    with pytest.raises(DriftProbeError, match="HTTP"):
        _probe(status=404)


def test_an_empty_200_is_retried_rather_than_reported_as_drift(monkeypatch):
    """The call site of the #271 fix, not just the helper.

    `request_with_retry` grew `retry_on_empty_body` for this probe specifically, and
    `test_http_retry.py` covers the helper thoroughly -- but deleting the one
    `retry_on_empty_body=True` argument here passed the entire suite. The fix and its
    only consumer were tested separately, so the wire between them was not tested at
    all.

    Without it the empty 200 reaches `resp.json()` and the probe raises
    `Non-JSON response`, which the drift bot files as an issue for what is
    a transient upstream blip -- exactly what #271 recorded.
    """
    monkeypatch.setattr("hvantk.core.utils.http.time.sleep", lambda _s: None)

    with requests_mock.Mocker() as m:
        m.get(
            PROBE_URL,
            [
                {"status_code": 200, "text": ""},
                {"status_code": 200, "json": _BODY, "headers": _LIVE_HEADERS},
            ],
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "2026_03"


def test_worst_case_retry_budget_fits_under_the_runner_timeout():
    """The probe's retry budget must fit under drift_runner's SIGALRM with room to
    spare, or a probe that legitimately exhausts every retry gets killed mid-request
    instead of raising the ordinary DriftProbeError callers already handle.

    Checked against `run_drift_checks`'s default. The `hvantk drift` CLI passes its own
    `--timeout` through explicitly; `hvantk/tests/test_drift_cli.py` pins that default
    to this one, so the check covers both (a skills test may not import `tools/`).

    Note: the 51 s bound covers the exponential-backoff path; a Retry-After sleep is
    clamped to max_sleep_s and can exceed it, and the drift runner's SIGALRM is the
    backstop in that case.
    """
    import inspect

    from hvantk.core.plugin import drift_runner
    from hvantk.core.utils.http import DEFAULT_BACKOFF_S, DEFAULT_MAX_SLEEP_S, _backoff
    from hvantk.skills.uniprot_ptm.drift_probe import _ATTEMPTS, _TIMEOUT_S

    runner_timeout = (
        inspect.signature(drift_runner.run_drift_checks).parameters["timeout"].default
    )
    worst_case = sum(
        _backoff(DEFAULT_BACKOFF_S, a, DEFAULT_MAX_SLEEP_S) for a in range(1, _ATTEMPTS)
    ) + _ATTEMPTS * sum(_TIMEOUT_S)
    assert worst_case < runner_timeout
