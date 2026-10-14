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

``tls_context()`` and ``stream_to_file()`` below are unrelated to retries.

``tls_context()`` is for the call sites that use ``urllib.request.urlopen`` directly
instead of ``requests``: clinvar's ``verify_md5``, the ``get_available_versions`` of
clingen, gencc and uniprot-ptm, and the ``urllib`` fallback in uniprot-ptm's
``_fetch_json``. On a python.org macOS build (this repo's base Python) the platform's
default trust store is not configured at all, so a bare ``urlopen`` HTTPS call fails
with ``CERTIFICATE_VERIFY_FAILED``, while ``requests`` (which bundles ``certifi``) is
unaffected. ``tls_context()`` adds certifi's bundle on top of whatever the platform
already trusts.

``stream_to_file()`` is the one atomic download primitive for every plugin
downloader that writes a single file to disk: ``hgnc``, ``gnomad-metrics``,
``ensembl-gene:structure``, ``gwas-catalog`` and ``insider:interfaces`` call it
directly, and ``hvantk.core.utils.file_utils.download_file`` (clinvar, clingen, gencc)
delegates to it. It goes through ``requests``, so it needs no ``tls_context()``.
"""

from __future__ import annotations

import logging
import os
import shutil
import ssl
import tempfile
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping

import certifi
import requests
from tqdm import tqdm

from hvantk.core.plugin.drift_runner import _current_umask

logger = logging.getLogger(__name__)

# 429 plus the transient 5xx family. 500 is included deliberately: several of the sources
# these probes hit answer overload with a bare 500 rather than a 503.
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})

#: The User-Agent hvantk identifies itself with -- the value the ClinGen, GenCC, ClinVar
#: and UniProt downloaders send. A probe that leaves requests' default
#: (``python-requests/<version>``) in place presents as a generic scraper to a site behind
#: bot protection (dbnsfp.org sits behind Cloudflare), and a challenge served to the CI
#: runner would fail that probe closed on every scheduled run (#408).
USER_AGENT = "hvantk/1.0"

DEFAULT_ATTEMPTS = 4
DEFAULT_BACKOFF_S = 2.0
DEFAULT_MAX_SLEEP_S = 30.0

# 2**32 s is ~136 years; any clamp is reached long before this.
_MAX_EXP = 32


@lru_cache(maxsize=1)
def tls_context() -> ssl.SSLContext:
    """Cached SSL context that adds certifi's CA bundle on top of the system store.

    ``ssl.create_default_context()`` with no arguments still calls
    ``load_default_certs()`` internally, so it keeps the platform trust store AND
    honours ``SSL_CERT_FILE``/``SSL_CERT_DIR`` -- unlike
    ``create_default_context(cafile=certifi.where())``, which uses *only* the given
    file and silently drops both. ``load_verify_locations()`` then layers certifi's
    bundle on top, so a host whose cert chains to a CA only certifi ships still
    verifies, without narrowing what a correctly configured system already trusted.
    This is deliberately not "the same way ``requests`` does it" -- ``requests``
    (via ``urllib3``) uses certifi exclusively, full stop. This is only a fallback
    for the platform gap ``requests`` happens not to have (python.org's macOS build
    ships no system trust store at all), so it adds rather than replaces. Cached
    because constructing an ``SSLContext`` re-parses both CA sources from disk on
    every call.
    """
    ctx = ssl.create_default_context()
    ctx.load_verify_locations(cafile=certifi.where())
    return ctx


_CHUNK_SIZE = 1024 * 1024  # 1 MiB


def stream_to_file(
    url: str,
    dest: str | Path,
    *,
    timeout: float,
    headers: Mapping[str, str] | None = None,
    progress_label: str | None = None,
) -> None:
    """Stream *url* to *dest* atomically, through a sibling temp file.

    Promoted from a private copy in ``hvantk/skills/gwas_catalog/cli.py`` (now
    deleted in favour of this one) to be the one download primitive every plugin
    downloader that writes a single file uses. A plain ``open(dest, "wb")`` would
    leave a truncated file at *dest* if the network drops mid-download; this writes
    to a temp file in *dest*'s directory first and only ``os.replace``s it into
    place once the transfer completes, so *dest* either does not exist or holds the
    complete download, never a partial one -- including on a short read (a response
    whose connection closes before its declared ``Content-Length`` is fully
    delivered). ``urllib3`` (which ``requests`` sits on) only enforces that itself
    since its 2.0 release (``enforce_content_length`` defaulting True); ``requests``
    allows ``urllib3>=1.21.1`` and this package pins no floor above that, so the
    check below is done explicitly rather than relied on from the dependency.
    Raises ``requests.exceptions.ChunkedEncodingError`` -- a ``RequestException``,
    so callers' existing "on download failure" messages stay right -- on a
    mismatch, before the temp file is renamed into place, so ``finally`` below
    still removes it. Uses ``requests`` rather than ``urllib.request``, so it is
    not subject to the trust-store gap ``tls_context()`` works around and needs no
    context of its own.

    Parameters
    ----------
    url : str
        URL to fetch.
    dest : str | Path
        Local path to write the response body to.
    timeout : float
        Socket timeout in seconds, applied to the connection and to each read -- not a
        cap on total transfer time, so a slow-but-steady transfer keeps succeeding and
        only a genuinely stalled one times out.
    headers : Mapping[str, str], optional
        Extra request headers (e.g. ``User-Agent``).
    progress_label : str, optional
        When given, show a ``tqdm`` progress bar with this label while streaming.
        Omitted by default (matching the behaviour this replaces at most call
        sites); ``hvantk.core.utils.file_utils.download_file`` passes its
        ``file_name`` here to keep its existing progress bar.

    Raises
    ------
    requests.exceptions.RequestException
        On a network, HTTP-level, or short-read failure. *dest* is left untouched
        (neither created nor modified) in every failure case.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(dest.parent), prefix=dest.name + ".")
    os.close(fd)
    done = False
    try:
        with requests.get(
            url, stream=True, timeout=timeout, headers=dict(headers or {})
        ) as response:
            response.raise_for_status()
            declared_length = response.headers.get("content-length")
            total = int(declared_length) if declared_length is not None else 0
            written = 0
            with open(tmp_name, "wb") as fh:
                if progress_label is not None:
                    with tqdm(
                        desc=progress_label,
                        total=total or None,
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                    ) as bar:
                        for chunk in response.iter_content(chunk_size=_CHUNK_SIZE):
                            if chunk:
                                size = fh.write(chunk)
                                written += size
                                bar.update(size)
                else:
                    for chunk in response.iter_content(chunk_size=_CHUNK_SIZE):
                        if chunk:
                            written += fh.write(chunk)

            # Short-read check (see the docstring for why it is not left to
            # urllib3). Only for an identity-encoded response: with a
            # Content-Encoding, Content-Length is the compressed size, while
            # `iter_content()` yields decoded bytes.
            if (
                declared_length is not None
                and "content-encoding" not in response.headers
                and written != total
            ):
                raise requests.exceptions.ChunkedEncodingError(
                    f"Response ended prematurely: {written} bytes read, "
                    f"Content-Length declared {total}"
                )
        # mkstemp always creates mode 0600; match what a normal file create (or the
        # existing file, on overwrite) would get, mirroring
        # drift_runner.write_fingerprint's atomic-publish idiom.
        if dest.exists():
            shutil.copymode(dest, tmp_name)
        else:
            os.chmod(tmp_name, 0o666 & ~_current_umask())
        os.replace(tmp_name, dest)
        done = True
    finally:
        if not done and os.path.exists(tmp_name):
            os.remove(tmp_name)


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


def _is_streamed(kwargs: dict, session: requests.Session | None) -> bool:
    """Whether this request is actually streamed, per-request kwarg OR session default.

    Checking only ``kwargs["stream"]`` reads a genuinely streamed request as buffered:
    ``Session.merge_environment_settings`` does ``stream = merge_setting(stream,
    self.stream)`` (requests 2.32.5), so ``session.stream = True`` with no per-request
    kwarg streams the response anyway. The empty-body check would then touch
    ``response.content`` and pull the whole body into memory -- measured at 84 MB peak
    for a 40 MB download that the per-request form keeps at 0.2 MB -- and set
    ``_content_consumed``, silently defeating the streaming the caller asked for on a
    path whose own comment says it avoids exactly that. No caller sets ``session.stream``
    today, but gnomad_metrics and gencc already pass ``session=``, so it is one line away.
    """
    if kwargs.get("stream"):
        return True
    return bool(session is not None and getattr(session, "stream", False))


def _backoff(backoff_s: float, attempt: int, max_sleep_s: float) -> float:
    """Exponential backoff for `attempt`, clamped to `max_sleep_s`.

    The exponent is capped before it is used. `backoff_s * 2 ** (attempt - 1)` looks
    harmless because min() clamps the result, but the multiplication happens first: with
    a large `attempts` it raises `OverflowError: int too large to convert to float`
    rather than merely being slow. _MAX_EXP is far past any reachable max_sleep_s, so
    capping it changes no real result.
    """
    return min(backoff_s * 2 ** min(attempt - 1, _MAX_EXP), max_sleep_s)


def request_with_retry(
    method: str,
    url: str,
    *,
    attempts: int = DEFAULT_ATTEMPTS,
    backoff_s: float = DEFAULT_BACKOFF_S,
    max_sleep_s: float = DEFAULT_MAX_SLEEP_S,
    retry_statuses: Iterable[int] = RETRY_STATUSES,
    retry_on_empty_body: bool = False,
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

    ``retry_on_empty_body`` additionally retries a ``200``/``203`` response (excluding
    HEAD/OPTIONS, which have no body by definition) whose body is empty or
    whitespace-only. That is opt-in because an empty 200 is legitimate for plenty of endpoints,
    and only the caller knows. It exists because the failure it covers is invisible to
    every status-based rule: on 2026-09-21 the medRxiv API served ``200`` with
    ``Content-Type: application/json`` and zero bytes, so the pqtl probe raised "returned
    non-JSON" and filed #352 -- a transient upstream fault wearing a success code, which
    no amount of 5xx retrying would have caught. Not applied to streamed requests, where
    touching ``.content`` would consume the body the caller asked to stream.

    If every attempt comes back empty the final response is still returned -- but with a
    warning, because unlike an exhausted 429 it will not raise for the caller.
    """
    if attempts < 1:
        raise ValueError(f"attempts must be >= 1, got {attempts}")
    # Rejected here rather than at the sleep: a negative value only surfaces AFTER a
    # transient failure, so the caller would see ValueError("sleep length must be
    # non-negative") in place of the upstream error it was retrying.
    if backoff_s < 0 or max_sleep_s < 0:
        raise ValueError(
            f"backoff_s and max_sleep_s must be >= 0, got {backoff_s} and {max_sleep_s}"
        )

    retry_statuses = frozenset(retry_statuses)
    caller = session if session is not None else requests

    for attempt in range(1, attempts + 1):
        is_last = attempt == attempts
        try:
            response = caller.request(method, url, **kwargs)
        except (requests.Timeout, requests.ConnectionError) as exc:
            if is_last:
                raise
            sleep_s = _backoff(backoff_s, attempt, max_sleep_s)
            logger.warning(
                "%s %s failed (%s); retrying in %.1fs (attempt %d/%d)",
                method.upper(),
                url,
                type(exc).__name__,
                sleep_s,
                attempt,
                attempts,
            )
            time.sleep(sleep_s)
            continue

        empty_body = (
            retry_on_empty_body
            and not _is_streamed(kwargs, session)
            # `response.ok` is merely status < 400, which makes 204 No Content, a 304
            # from a conditional GET, and every redirect look like an "empty success" --
            # and a HEAD response has no body by definition. Retrying any of those burns
            # the full backoff and returns the identical response. gnomad_metrics
            # already routes a HEAD through this helper (without opting in today), so
            # the guard is one caller away from mattering rather than hypothetical.
            and method.upper() not in {"HEAD", "OPTIONS"}
            # An allowlist, not `2xx and not 204`: the rest of 2xx is *expected* to be
            # bodiless or non-representational, so a blanket range retries four codes
            # that can never improve -- 201 Created and 202 Accepted routinely answer
            # empty (202 is the standard async-job-submitted reply), 205 Reset Content
            # MUST NOT carry content per RFC 9110 s15.3.6, and an empty 206 is a range
            # the server had nothing for. Each would burn the full backoff to return
            # the identical response: the exact waste the 204 exclusion exists to stop.
            # 203 is in because it is a 200 relayed by a transforming proxy.
            and response.status_code in {200, 203}
            and not response.content.strip()
        )
        if is_last or (response.status_code not in retry_statuses and not empty_body):
            # `is_last` short-circuits, so an exhausted empty-body retry would otherwise
            # return in silence. The docstring's defence -- "returns the final Response
            # without raise_for_status, so callers keep their existing error handling" --
            # holds for statuses (a persistent 429 still raises) but NOT here: a 200
            # passes raise_for_status, so the caller receives a success indistinguishable
            # from a legitimately-empty endpoint after every attempt failed. Say so, or
            # this is #352 moved four attempts later rather than fixed.
            if is_last and empty_body:
                logger.warning(
                    "%s %s still returned %d with an empty body after %d attempt(s); "
                    "returning it anyway -- the caller sees a success-shaped response",
                    method.upper(),
                    url,
                    response.status_code,
                    attempts,
                )
            return response

        # Retryable, with attempts left. Read Retry-After BEFORE closing: a streamed
        # response holds its connection until closed or consumed, and leaving discarded
        # attempts open would leak one connection per retry.
        retry_after = parse_retry_after(response.headers.get("Retry-After"))
        response.close()
        backoff = _backoff(backoff_s, attempt, max_sleep_s)
        sleep_s = min(retry_after, max_sleep_s) if retry_after is not None else backoff
        logger.warning(
            "%s %s returned %d%s; retrying in %.1fs (attempt %d/%d)",
            method.upper(),
            url,
            response.status_code,
            " with an empty body" if empty_body else "",
            sleep_s,
            attempt,
            attempts,
        )
        time.sleep(sleep_s)

    raise AssertionError("unreachable: loop returns or raises on the last attempt")
