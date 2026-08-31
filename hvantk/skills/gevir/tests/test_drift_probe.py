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
_HEADERS = {
    "Content-Length": "10270511",
    "ETag": '"6423adf134a669acc357f619d1162009"',
    "Last-Modified": "Tue, 14 Nov 2023 18:23:50 GMT",
}


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_HEADERS)
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 1
    # Springer serves an MD5 over the body, so the ETag is a true content digest
    # and belongs in the compared surface, unquoted.
    assert fp["checksums"][_FILENAME] == "6423adf134a669acc357f619d1162009"
    assert fp["extras"]["content_length"] == "10270511"
    assert "fetched_at" in fp


def test_last_modified_is_demoted_out_of_the_compared_surface():
    """A byte-identical republish must not open a pull request (hgnc precedent)."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers=_HEADERS)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert fp["informational"]["last_modified"] == "Tue, 14 Nov 2023 18:23:50 GMT"


def test_no_content_signal_fails_closed():
    """With neither validator the fingerprint would compare equal forever."""
    with requests_mock.Mocker() as m:
        m.head(GEVIR_SUPPLEMENTARY_URL, headers={})
        with pytest.raises(DriftProbeError, match="no content signal"):
            fetch_fingerprint()


def test_etag_alone_is_sufficient():
    """Content-Length may be absent if a real content digest is still available."""
    with requests_mock.Mocker() as m:
        m.head(
            GEVIR_SUPPLEMENTARY_URL,
            headers={"ETag": '"6423adf134a669acc357f619d1162009"'},
        )
        fp = fetch_fingerprint()

    assert fp["checksums"][_FILENAME] == "6423adf134a669acc357f619d1162009"
    assert fp["extras"]["content_length"] is None
