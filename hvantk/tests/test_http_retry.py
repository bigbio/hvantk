"""Guard the drift-probe HTTP retry helper.

Regression: on 2026-08-04 GenCC answered the scheduled drift regeneration with HTTP 429.
No probe retried, so it raised, no PR could be opened for the drifted dataset, and the
whole workflow run failed -- with nothing wrong with either the data or the code.

The clamp is tested as carefully as the retry because it is the reason this helper exists
instead of ``urllib3.util.Retry``: that honours ``Retry-After`` via ``time.sleep``
with no upper bound, so a host answering ``Retry-After: 3600`` would park CI for an hour.
Deliberately non-Hail so it runs in the default suite.
"""
from __future__ import annotations

import pytest
import requests

from hvantk.core.utils import http as http_util

URL = "https://example.invalid/drift.tsv"


@pytest.fixture
def slept(monkeypatch):
    """Record sleep durations instead of actually sleeping."""
    calls: list[float] = []
    monkeypatch.setattr(http_util.time, "sleep", calls.append)
    return calls


def test_retries_past_a_429_and_returns_the_success(requests_mock, slept):
    """The GenCC failure, end to end: transient 429 then a good response."""
    requests_mock.get(URL, [{"status_code": 429}, {"status_code": 200, "text": "ok"}])

    resp = http_util.request_with_retry("GET", URL)

    assert resp.status_code == 200
    assert resp.text == "ok"
    assert len(slept) == 1, "should have backed off exactly once"


def test_retry_after_is_honoured_but_clamped(requests_mock, slept):
    """A long Retry-After must not be able to stall the job.

    This is the specific behaviour urllib3 does NOT give us.
    """
    requests_mock.get(
        URL,
        [
            {"status_code": 429, "headers": {"Retry-After": "3600"}},
            {"status_code": 200, "text": "ok"},
        ],
    )

    resp = http_util.request_with_retry("GET", URL, max_sleep_s=30.0)

    assert resp.status_code == 200
    assert slept == [30.0], "Retry-After should be clamped to max_sleep_s, not obeyed"


def test_exhausted_retries_return_the_last_response(requests_mock, slept):
    """Callers keep their existing error handling.

    The helper never calls raise_for_status itself, so a response that is still 429
    after the final attempt comes back as a 429 and raises exactly where it used to.
    """
    requests_mock.get(URL, status_code=429)

    resp = http_util.request_with_retry("GET", URL, attempts=3)

    assert resp.status_code == 429
    assert len(slept) == 2, "3 attempts means 2 gaps"
    with pytest.raises(requests.HTTPError):
        resp.raise_for_status()


def test_connection_errors_are_retried_but_other_request_errors_are_not(
    requests_mock, slept
):
    """A flaky handshake is the same class of problem as a 503; a bad URL is not."""
    requests_mock.get(URL, [{"exc": requests.ConnectionError}, {"text": "ok"}])
    assert http_util.request_with_retry("GET", URL).text == "ok"
    assert len(slept) == 1

    slept.clear()
    requests_mock.get(URL, exc=requests.URLRequired)
    with pytest.raises(requests.URLRequired):
        http_util.request_with_retry("GET", URL)
    assert slept == [], "non-transient errors must not be retried"


@pytest.mark.parametrize(
    "header,expected",
    [
        ("120", 120.0),
        ("  90  ", 90.0),
        ("not-a-date", None),
        ("", None),
        (None, None),
        ("Wed, 21 Oct 1990 07:28:00 GMT", 0.0),  # in the past -> clamped to 0, not negative
    ],
)
def test_parse_retry_after(header, expected):
    assert http_util.parse_retry_after(header) == expected
