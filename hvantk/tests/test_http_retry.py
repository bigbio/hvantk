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

import threading

import pytest
import requests

from hvantk.core.utils import http as http_util

URL = "https://example.invalid/drift.tsv"


@pytest.fixture
def slept(monkeypatch):
    """Record this thread's sleep durations instead of actually sleeping.

    ``monkeypatch.setattr(http_util.time, "sleep", ...)`` rebinds the attribute
    on the shared ``time`` module, so it intercepts *every* thread in the
    process, not just the retry loop under test. Whenever the full suite leaves
    a JVM/py4j polling thread alive, its ``time.sleep(1)`` calls land in this
    list and the length assertions fail for a reason unrelated to backoff --
    observed at 625,818 captured sleeps once and 31,085 another time, while the
    same tests pass in isolation. Filtering by thread keeps the assertions
    about the retry loop, which always runs on the calling thread.
    """
    calls: list[float] = []
    test_thread = threading.get_ident()
    real_sleep = http_util.time.sleep

    def record(seconds: float) -> None:
        if threading.get_ident() == test_thread:
            calls.append(seconds)
            return
        # Every other thread still really sleeps. Swallowing their sleeps
        # instead would turn a JVM/py4j poller's `while True: ...; sleep(1)`
        # into a hot spin on a core for the duration of these tests -- quiet,
        # but worse than the noisy failure this fixture is fixing.
        real_sleep(seconds)

    monkeypatch.setattr(http_util.time, "sleep", record)
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


# --- review findings on #268 ---------------------------------------------------------


@pytest.mark.parametrize("kwargs", [
    {"backoff_s": -1.0},
    {"max_sleep_s": -5.0},
])
def test_negative_timings_are_rejected_before_the_first_request(kwargs):
    """A negative value would otherwise surface only AFTER a transient failure, replacing
    the upstream error the caller was retrying with ValueError('sleep length must be
    non-negative')."""
    with pytest.raises(ValueError, match="must be >= 0"):
        http_util.request_with_retry("GET", URL, **kwargs)


def test_backoff_is_capped_at_max_sleep_for_any_attempt():
    """`backoff_s * 2 ** (attempt - 1)` is evaluated BEFORE min() clamps it, so a high
    attempt NUMBER raises `OverflowError: int too large to convert to float` -- not
    merely a slow computation, an outright crash.

    Asserted on the helper rather than by driving 1,200 real failures through the retry
    loop: the `slept` fixture patches `time.sleep` globally, so under a full-suite run
    that version captured sleeps from unrelated code (625,818 of them) and failed for a
    reason having nothing to do with backoff. This form is deterministic and still fails
    if the exponent cap is removed.
    """
    assert http_util._backoff(2.0, 1, 30.0) == 2.0
    assert http_util._backoff(2.0, 4, 30.0) == 16.0
    assert http_util._backoff(2.0, 99, 30.0) == 30.0
    assert http_util._backoff(2.0, 10**6, 30.0) == 30.0
