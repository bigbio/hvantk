"""Retrying HTTP request helper for drift probes and downloaders.

Thirteen of the twenty-two drift probes call ``requests`` directly and none of them
retried, so a single transient upstream response could fail the whole scheduled drift
run. That is exactly what happened on 2026-08-04: GenCC answered the regenerate request
with HTTP 429, the probe raised ``DriftProbeError``, no PR could be opened for the
drifted dataset, and the workflow exited non-zero -- while nothing was actually wrong
with either the data or the code.

Why not ``urllib3.util.Retry`` mounted on an ``HTTPAdapter``, which is the obvious
answer: it honours ``Retry-After`` through ``Retry.sleep_for_retry``, which calls
``time.sleep(retry_after)`` with **no upper bound**. ``backoff_max`` (120s by default)
caps only the exponential-backoff path, not this one -- verified against urllib3 2.6.3.
A rate-limited host answering ``Retry-After: 3600`` would therefore park a CI job for an
hour. This helper honours the header, because ignoring a server's explicit backoff
request is how you get rate-limited harder, but clamps it to ``max_sleep_s``.
"""
from __future__ import annotations

import logging
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from typing import Iterable

import requests

logger = logging.getLogger(__name__)

# 429 plus the transient 5xx family. 500 is included deliberately: several of the sources
# these probes hit answer overload with a bare 500 rather than a 503.
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})

DEFAULT_ATTEMPTS = 4
DEFAULT_BACKOFF_S = 2.0
DEFAULT_MAX_SLEEP_S = 30.0


def parse_retry_after(value: str | None) -> float | None:
    """``Retry-After`` as seconds, or ``None`` if absent/unparseable.

    RFC 9110 permits both forms -- delta-seconds (``120``) and an HTTP-date
    (``Wed, 21 Oct 2026 07:28:00 GMT``) -- and real servers send both. A date in the
    past clamps to 0 rather than going negative.
    """
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    try:
        return max(0.0, float(value))
    except ValueError:
        pass
    try:
        when = parsedate_to_datetime(value)
    except (TypeError, ValueError):
        return None
    if when is None:
        return None
    if when.tzinfo is None:  # HTTP-dates are GMT; naive means UTC here
        when = when.replace(tzinfo=timezone.utc)
    return max(0.0, (when - datetime.now(timezone.utc)).total_seconds())


def request_with_retry(
    method: str,
    url: str,
    *,
    attempts: int = DEFAULT_ATTEMPTS,
    backoff_s: float = DEFAULT_BACKOFF_S,
    max_sleep_s: float = DEFAULT_MAX_SLEEP_S,
    retry_statuses: Iterable[int] = RETRY_STATUSES,
    session: requests.Session | None = None,
    **kwargs,
) -> requests.Response:
    """``requests.request`` with retries on transient statuses and connection errors.

    Returns the final :class:`requests.Response` **without** calling
    ``raise_for_status``, so callers keep their existing error handling: a response that
    is still 429 after the last attempt comes back as a 429 and raises where it always
    did. Only the transient case changes behaviour.

    Connection errors and timeouts are retried too -- they are the same class of problem
    as a 503, and a probe that dies on one flaky TCP handshake is the bug this fixes.
    Other ``RequestException`` subclasses (a malformed URL, say) are not retried, since
    repeating them cannot help.
    """
    if attempts < 1:
        raise ValueError(f"attempts must be >= 1, got {attempts}")

    retry_statuses = frozenset(retry_statuses)
    caller = session if session is not None else requests

    for attempt in range(1, attempts + 1):
        is_last = attempt == attempts
        try:
            response = caller.request(method, url, **kwargs)
        except (requests.Timeout, requests.ConnectionError) as exc:
            if is_last:
                raise
            sleep_s = min(backoff_s * 2 ** (attempt - 1), max_sleep_s)
            logger.warning(
                "%s %s failed (%s); retrying in %.1fs (attempt %d/%d)",
                method.upper(), url, type(exc).__name__, sleep_s, attempt, attempts,
            )
            time.sleep(sleep_s)
            continue

        if is_last or response.status_code not in retry_statuses:
            return response

        # Retryable, with attempts left. Read Retry-After BEFORE closing: a streamed
        # response holds its connection until closed or consumed, and leaving discarded
        # attempts open would leak one connection per retry.
        retry_after = parse_retry_after(response.headers.get("Retry-After"))
        response.close()
        backoff = backoff_s * 2 ** (attempt - 1)
        sleep_s = min(retry_after if retry_after is not None else backoff, max_sleep_s)
        logger.warning(
            "%s %s returned %d; retrying in %.1fs (attempt %d/%d)",
            method.upper(), url, response.status_code, sleep_s, attempt, attempts,
        )
        time.sleep(sleep_s)

    raise AssertionError("unreachable: loop returns or raises on the last attempt")
