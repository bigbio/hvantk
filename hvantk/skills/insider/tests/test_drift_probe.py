"""INSIDER drift probes should fingerprint each product independently.

Runs OFFLINE via requests_mock, so CI never hits the live portal. The probes
replaced a stub sentinel (issue #177) once the portal's direct file paths were
confirmed addressable; these tests pin the properties that make them real
comparators rather than false-greens.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.insider.drift_probe import (
    INSIDER_BED_FILENAME,
    INSIDER_BED_URL,
    fetch_fingerprint,
)
from hvantk.skills.insider.interfaces.drift_probe import (
    INSIDER_INTERFACES_FILENAME,
    INSIDER_INTERFACES_URL,
)
from hvantk.skills.insider.interfaces.drift_probe import (
    fetch_fingerprint as fetch_interfaces_fingerprint,
)

_BED_HEADERS = {
    "Content-Length": "1171848841",
    "ETag": '"45d8fe89-566adb9ac3780"',
    "Last-Modified": "Mon, 05 Mar 2018 17:33:34 GMT",
}
_INTERFACES_HEADERS = {
    "Content-Length": "49762981",
    "ETag": '"2f752a5-61883677f7fa2"',
    "Last-Modified": "Wed, 15 May 2024 19:48:36 GMT",
}


def test_bed_probe_fingerprints_only_its_own_file():
    """The two products are versioned independently, so a variants fingerprint
    must not carry the interfaces file -- otherwise an interfaces-only change
    rewrites the provenance of every BED-derived artifact."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_BED_URL, headers=_BED_HEADERS)
        fp = fetch_fingerprint()

    assert set(fp["headers"]) == {INSIDER_BED_FILENAME}
    assert INSIDER_INTERFACES_FILENAME not in fp["headers"]
    assert fp["headers"][INSIDER_BED_FILENAME]["content_length"] == "1171848841"


def test_interfaces_probe_fingerprints_only_its_own_file():
    with requests_mock.Mocker() as m:
        m.head(INSIDER_INTERFACES_URL, headers=_INTERFACES_HEADERS)
        fp = fetch_interfaces_fingerprint()

    assert set(fp["headers"]) == {INSIDER_INTERFACES_FILENAME}
    assert INSIDER_BED_FILENAME not in fp["headers"]


def test_a_bed_outage_cannot_mask_interfaces_drift():
    """Sharing one probe meant a fault on the frozen-since-2018 BED aborted before
    the interfaces file was ever reached. Split probes must be independent."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_BED_URL, status_code=503)
        m.head(INSIDER_INTERFACES_URL, headers=_INTERFACES_HEADERS)

        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
        # The interfaces probe is unaffected.
        fp = fetch_interfaces_fingerprint()

    assert fp["headers"][INSIDER_INTERFACES_FILENAME]["content_length"] == "49762981"


def test_etag_is_in_the_compared_surface():
    """The ETag is hex(size)-hex(mtime), a strict superset of Content-Length, so it
    catches an equal-size edit that Content-Length alone cannot. It must sit under
    `headers` (compared), not `informational` (ignored by drift comparison)."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_INTERFACES_URL, headers=_INTERFACES_HEADERS)
        fp = fetch_interfaces_fingerprint()

    assert fp["headers"][INSIDER_INTERFACES_FILENAME]["etag"] == "2f752a5-61883677f7fa2"
    assert "etag" not in fp.get("informational", {}).get(
        INSIDER_INTERFACES_FILENAME, {}
    )


def test_equal_size_edit_is_detected():
    """The regression this probe exists to catch: same length, changed content."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_INTERFACES_URL, headers=_INTERFACES_HEADERS)
        before = fetch_interfaces_fingerprint()

    with requests_mock.Mocker() as m:
        m.head(
            INSIDER_INTERFACES_URL,
            headers={**_INTERFACES_HEADERS, "ETag": '"2f752a5-99999999999999"'},
        )
        after = fetch_interfaces_fingerprint()

    assert before["headers"] != after["headers"]


def test_schema_signal_is_populated_so_drift_is_not_auto_batched():
    """`drift_to_pr.classify_risk` tiers a diff as "routine" unless `headers` or
    `checksums` moved. With both empty, a changed `track name=` format would be
    swept into a batch telling the reviewer the schema signal was unchanged."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_BED_URL, headers=_BED_HEADERS)
        fp = fetch_fingerprint()

    from hvantk.core.plugin.api import PROBE_FINGERPRINT_IGNORED_KEYS

    schema_keys = {"headers", "checksums"}
    assert not schema_keys & PROBE_FINGERPRINT_IGNORED_KEYS
    assert fp["headers"], "headers must carry signal or every change reads as routine"


def test_missing_content_length_fails_closed():
    """Regression guard for the real failure hit while seeding the baseline: the
    server gzips text/plain on the fly and then omits Content-Length entirely."""
    with requests_mock.Mocker() as m:
        m.head(INSIDER_INTERFACES_URL, headers={"ETag": '"abc"'})
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_interfaces_fingerprint()


def test_empty_etag_fails_closed():
    """`ETag: ""` would normalize to '' and be stored as an empty digest, which
    placeholder_baseline_reason then reads as a hand-seeded baseline -- trapping
    the dataset in a probe_failed loop that regenerating cannot clear."""
    with requests_mock.Mocker() as m:
        m.head(
            INSIDER_INTERFACES_URL,
            headers={"Content-Length": "49762981", "ETag": '""'},
        )
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_interfaces_fingerprint()


def test_compressed_response_fails_closed():
    """Identity is asserted request-side; a proxy may compress anyway, and the
    recorded length would then describe the compressed body."""
    with requests_mock.Mocker() as m:
        m.head(
            INSIDER_INTERFACES_URL,
            headers={
                "Content-Length": "123",
                "ETag": '"abc"',
                "Content-Encoding": "gzip",
            },
        )
        with pytest.raises(DriftProbeError, match="gzip-encoded"):
            fetch_interfaces_fingerprint()


@pytest.mark.parametrize(
    ("probe", "url", "headers"),
    [
        (fetch_fingerprint, INSIDER_BED_URL, _BED_HEADERS),
        (fetch_interfaces_fingerprint, INSIDER_INTERFACES_URL, _INTERFACES_HEADERS),
    ],
)
def test_probe_requests_identity_encoding(probe, url, headers):
    """Without this the server compresses and Content-Length disappears."""
    with requests_mock.Mocker() as m:
        m.head(url, headers=headers)
        probe()
        assert len(m.request_history) == 1
        assert m.request_history[0].headers.get("Accept-Encoding") == "identity"
