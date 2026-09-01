"""gnomad-metrics drift probe should fingerprint every constraint object offline.

Runs OFFLINE via requests_mock. The probe replaced a stub sentinel (issue #177,
which rated the source "marginally feasible") once the public GCS objects were
confirmed to answer a HEAD with an MD5 ETag.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.gnomad_metrics.drift_probe import OBJECT_PATHS, fetch_fingerprint
from hvantk.skills.gnomad_metrics.shared.constants import GNOMAD_RELEASE_BASE_URL

_BY_GENE = "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"


def _url(path):
    return f"{GNOMAD_RELEASE_BASE_URL}/{path}"


def _headers(n):
    """Distinct per object, so a probe that mixed them up would be caught."""
    return {
        "ETag": f'"{n:032x}"',
        "Content-Length": str(1000 + n),
        "x-goog-generation": str(1_500_000_000_000_000 + n),
        "Last-Modified": f"Thu, 20 Aug 2020 15:19:{n:02d} GMT",
    }


def _mock_all(m, overrides=None):
    for n, path in enumerate(OBJECT_PATHS):
        m.head(_url(path), headers=(overrides or {}).get(path, _headers(n)))


def test_each_object_lands_under_its_own_key():
    """Distinct fixtures per object, so assigning one object's validators to all
    three would fail here."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert set(fp["headers"]) == set(OBJECT_PATHS)
    for n, path in enumerate(OBJECT_PATHS):
        assert fp["headers"][path]["content_length"] == str(1000 + n)
        assert fp["headers"][path]["etag"] == f"{n:032x}"


def test_signals_are_compared_not_stashed_in_checksums():
    """`_conventions` §12 defines `checksums` as a sha256 over the bytes used to
    derive `headers`. This probe fetches no body, so recording a raw validator
    there would misrepresent the contract."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert fp["checksums"] == {}
    assert fp["headers"], "signal must still be compared"


def test_generation_change_is_visible():
    """x-goog-generation moves on every rewrite even when bytes are identical."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        before = fetch_fingerprint()

    bumped = {OBJECT_PATHS[0]: {**_headers(0), "x-goog-generation": "9999"}}
    with requests_mock.Mocker() as m:
        _mock_all(m, overrides=bumped)
        after = fetch_fingerprint()

    assert before["headers"] != after["headers"]


def test_last_modified_is_demoted_out_of_the_compared_surface():
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert fp["informational"][_BY_GENE].startswith("Thu, 20 Aug 2020")
    assert "last_modified" not in fp["headers"][_BY_GENE]


@pytest.mark.parametrize(
    ("label", "bad_headers"),
    [
        ("no headers at all", {}),
        ("content-length only", {"Content-Length": "4609488"}),
        ("etag only", {"ETag": '"abc"'}),
        ("empty etag", {"Content-Length": "4609488", "ETag": '""'}),
    ],
)
def test_partial_validators_fail_closed(label, bad_headers):
    """Each discriminating case, not just the both-missing one. An ETag-only
    response used to record a null length; a Content-Length-only response used to
    drop the object out of the compared surface entirely."""
    with requests_mock.Mocker() as m:
        _mock_all(m, overrides={OBJECT_PATHS[0]: bad_headers})
        with pytest.raises(DriftProbeError, match="no usable content signal"):
            fetch_fingerprint()


def test_compressed_response_fails_closed():
    bad = {**_headers(0), "Content-Encoding": "gzip"}
    with requests_mock.Mocker() as m:
        _mock_all(m, overrides={OBJECT_PATHS[0]: bad})
        with pytest.raises(DriftProbeError, match="gzip-encoded"):
            fetch_fingerprint()


def test_one_object_failing_still_reports_the_others():
    """Raising on the first failure meant a fault on either frozen v2.1.1 object
    aborted before v4.0 -- the current release -- was ever probed. The error must
    name what failed rather than stopping at the first."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        m.head(_url(OBJECT_PATHS[0]), status_code=503)
        with pytest.raises(DriftProbeError) as excinfo:
            fetch_fingerprint()

    msg = str(excinfo.value)
    assert f"1 of {len(OBJECT_PATHS)}" in msg
    assert OBJECT_PATHS[0] in msg


def test_probe_requests_identity_encoding():
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fetch_fingerprint()
        assert len(m.request_history) == len(OBJECT_PATHS)
        assert all(
            r.headers.get("Accept-Encoding") == "identity" for r in m.request_history
        )
