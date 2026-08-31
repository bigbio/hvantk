"""INSIDER drift probe should fingerprint both products with the documented shape.

Runs OFFLINE via requests_mock, so CI never hits the live portal. The probe
replaced a stub sentinel (issue #177) once the portal's direct file paths were
confirmed addressable; these tests pin the properties that make it a real
comparator rather than a false-green.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.insider.drift_probe import (
    INSIDER_BED_URL,
    INSIDER_INTERFACES_URL,
    fetch_fingerprint,
)

_BED = "Whole_Human_Interactome_Interface_hg38.bed"
_INTERFACES = "H_sapiens_interfacesALL.txt"


def _mock_both(m, *, bed_headers=None, interfaces_headers=None):
    m.head(
        INSIDER_BED_URL,
        headers=bed_headers
        if bed_headers is not None
        else {
            "Content-Length": "1171848841",
            "ETag": '"45d8fe89-566adb9ac3780"',
            "Last-Modified": "Mon, 05 Mar 2018 17:33:34 GMT",
        },
    )
    m.head(
        INSIDER_INTERFACES_URL,
        headers=interfaces_headers
        if interfaces_headers is not None
        else {
            "Content-Length": "49762981",
            "ETag": '"2f752a5-61883677f7fa2"',
            "Last-Modified": "Wed, 15 May 2024 19:48:36 GMT",
        },
    )


def test_fetch_fingerprint_covers_both_products():
    """One probe serves two datasets, so both files must appear in the fingerprint."""
    with requests_mock.Mocker() as m:
        _mock_both(m)
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 1
    assert fp["extras"][_BED] == "1171848841"
    assert fp["extras"][_INTERFACES] == "49762981"
    assert "fetched_at" in fp


def test_etag_and_last_modified_are_demoted_out_of_the_compared_surface():
    """nginx derives the ETag from (size, mtime), so it inherits Last-Modified's
    failure mode: a byte-identical re-upload would open a no-op pull request.
    Both belong in `informational`, which drift comparison ignores."""
    with requests_mock.Mocker() as m:
        _mock_both(m)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert fp["informational"][_INTERFACES]["etag"] == "2f752a5-61883677f7fa2"
    assert fp["informational"][_BED]["last_modified"] == "Mon, 05 Mar 2018 17:33:34 GMT"
    # The compared surface must carry no validator that moves on a re-upload.
    assert "etag" not in fp["extras"]


def test_missing_content_length_fails_closed():
    """Content-Length is the only compared signal; without it the fingerprint would
    compare equal to its baseline forever, which is the false-green the stub avoided.

    Regression guard for the real failure hit while seeding the baseline: nginx
    gzips text/plain on the fly and then omits Content-Length entirely.
    """
    with requests_mock.Mocker() as m:
        _mock_both(m, interfaces_headers={"Content-Encoding": "gzip"})
        with pytest.raises(DriftProbeError, match="omitted Content-Length"):
            fetch_fingerprint()


def test_probe_requests_identity_encoding():
    """The probe must opt out of compression, or the server drops Content-Length."""
    with requests_mock.Mocker() as m:
        _mock_both(m)
        fetch_fingerprint()
        assert all(
            r.headers.get("Accept-Encoding") == "identity" for r in m.request_history
        )
