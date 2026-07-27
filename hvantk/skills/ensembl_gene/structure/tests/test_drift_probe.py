"""The ensembl-gene:structure probe must actually read the upstream server.

Runs OFFLINE: requests_mock stubs the HEAD so CI never touches ftp.ensembl.org.

The bug these tests pin is not a wrong value but an absent request. The probe used
to return ``{"release", "url"}`` -- two repo constants -- and therefore compared equal
to its own committed baseline no matter what Ensembl did, while its docstring claimed
it was a real fingerprint. A probe that cannot fail cannot detect anything.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.plugin.drift_runner import _compare_fingerprints
from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)
from hvantk.skills.ensembl_gene.structure.drift_probe import fetch_fingerprint

_LIVE_HEADERS = {
    "ETag": '"3d2b9d9-61ffc57e3eacb"',
    "Last-Modified": "Sun, 18 Aug 2024 22:02:07 GMT",
    "Content-Length": "64141785",
}


def test_fingerprint_records_live_validators():
    """Every field must come from the response, not from a repo constant.

    ``release`` is the one exception and is kept deliberately: bumping the pin without
    rebuilding should register as drift.
    """
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers=_LIVE_HEADERS)
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 2
    assert fp["source_version"] == "Sun, 18 Aug 2024 22:02:07 GMT"
    assert fp["headers"][ENSEMBL_GTF_FILENAME] == {
        "etag": '"3d2b9d9-61ffc57e3eacb"',
        "content_length": "64141785",
    }
    assert fp["extras"]["release"] == ENSEMBL_RELEASE
    assert "fetched_at" in fp


def test_probe_actually_calls_upstream():
    """The probe must issue a request, not synthesise a fingerprint offline."""
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers=_LIVE_HEADERS)
        fetch_fingerprint()

    assert m.call_count == 1
    assert m.request_history[0].method == "HEAD"
    assert m.request_history[0].url == ENSEMBL_GTF_URL


def test_upstream_failure_becomes_a_probe_error():
    """An unreachable upstream is probe_failed, never a silent clean or drifted."""
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, status_code=404)
        with pytest.raises(DriftProbeError, match="HTTP failure"):
            fetch_fingerprint()


def test_changed_upstream_file_is_detected():
    """A re-issued GTF (new ETag/size) must compare as drift, not clean."""
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers=_LIVE_HEADERS)
        baseline = fetch_fingerprint()

    reissued = {**_LIVE_HEADERS, "ETag": '"reissued"', "Content-Length": "64141999"}
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers=reissued)
        observed = fetch_fingerprint()

    assert _compare_fingerprints(baseline, observed) is not None


def test_unchanged_upstream_file_is_clean():
    """The other half of the contract: a real probe must not report spurious drift.

    Ensembl's archive serves a static file, so unlike ClinGen there is no request-time
    validator here to filter out -- two probes of an unchanged file must simply agree.
    """
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers=_LIVE_HEADERS)
        first = fetch_fingerprint()
        second = fetch_fingerprint()

    assert _compare_fingerprints(first, second) is None


def test_server_without_validators_fails_loudly():
    """Neither ETag nor Content-Length means nothing stable to compare."""
    with requests_mock.Mocker() as m:
        m.head(ENSEMBL_GTF_URL, headers={"Last-Modified": "x"})
        with pytest.raises(DriftProbeError, match="nothing stable"):
            fetch_fingerprint()
