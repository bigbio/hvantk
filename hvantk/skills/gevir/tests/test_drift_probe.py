"""GeVIR drift probe should fingerprint the supplementary object offline.

Runs OFFLINE via requests_mock, so CI never hits Springer. The probe replaced a
stub sentinel (issue #177) once the article's ESM object was confirmed to answer
a HEAD with both Content-Length and a content-hash ETag.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.gevir.drift_probe import GEVIR_SUPPLEMENTARY_URL, fetch_fingerprint

_FILENAME = "41588_2019_560_MOESM3_ESM.xlsx"
_RESPONSE_HEADERS = {
    "Content-Length": "10270511",
    "ETag": '"6423adf134a669acc357f619d1162009"',
    "Last-Modified": "Tue, 14 Nov 2023 18:23:50 GMT",
}


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_RESPONSE_HEADERS)
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 2
    assert fp["headers"][_FILENAME]["etag"] == "6423adf134a669acc357f619d1162009"
    assert fp["headers"][_FILENAME]["content_length"] == "10270511"


def test_content_signal_is_not_recorded_as_a_schema_checksum():
    """`drift_to_pr` reads `checksums` as a hash of the column-header row, so a raw
    ETag there tiers every routine content update as a SCHEMA change and tells the
    reviewer to check builder.py. This probe fetches no body and has no schema
    signal to offer, so `checksums` must stay empty."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_RESPONSE_HEADERS)
        fp = fetch_fingerprint()

    assert fp["checksums"] == {}
    assert fp["headers"], "signal must still be compared, just not as a schema hash"


def test_last_modified_is_demoted_out_of_the_compared_surface():
    """A byte-identical republish must not open a pull request (hgnc precedent)."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_RESPONSE_HEADERS)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert fp["informational"]["last_modified"] == "Tue, 14 Nov 2023 18:23:50 GMT"
    # The negative half: the timestamp must appear nowhere that drift compares.
    assert "last_modified" not in fp["headers"][_FILENAME]
    assert "last_modified" not in fp["checksums"]


def test_weak_etag_normalizes_to_the_same_tag_as_the_strong_form():
    """CDNs flip strong<->weak on a byte-identical object. `strip('"')` left the
    `W/` prefix attached, so the mangled value read as drift and was then committed
    as the new baseline."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_RESPONSE_HEADERS)
        strong = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.head(
            GEVIR_SUPPLEMENTARY_URL,
            headers={
                **_RESPONSE_HEADERS,
                "ETag": 'W/"6423adf134a669acc357f619d1162009"',
            },
        )
        weak = fetch_fingerprint()

    assert weak["headers"] == strong["headers"]


def test_no_content_signal_fails_closed():
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers={})
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_fingerprint()


def test_etag_only_response_fails_closed():
    """Accepting it would record `content_length: None`; the bot commits that as the
    new baseline and every later normal response then reads as drift forever."""
    with requests_mock.Mocker() as m:
        m.head(
            GEVIR_SUPPLEMENTARY_URL,
            headers={"ETag": '"6423adf134a669acc357f619d1162009"'},
        )
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_fingerprint()


def test_empty_etag_fails_closed():
    with requests_mock.Mocker() as m:
        m.head(
            GEVIR_SUPPLEMENTARY_URL,
            headers={"Content-Length": "10270511", "ETag": '""'},
        )
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_fingerprint()


def test_probe_requests_identity_encoding():
    """Without this a compressing server drops Content-Length and suffixes the ETag."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_RESPONSE_HEADERS)
        fetch_fingerprint()
        assert len(m.request_history) == 1
        assert m.request_history[0].headers.get("Accept-Encoding") == "identity"


def test_compressed_response_fails_closed():
    with requests_mock.Mocker() as m:
        m.head(
            GEVIR_SUPPLEMENTARY_URL,
            headers={
                "Content-Length": "123",
                "ETag": '"abc"',
                "Content-Encoding": "gzip",
            },
        )
        with pytest.raises(DriftProbeError, match="gzip-encoded"):
            fetch_fingerprint()
