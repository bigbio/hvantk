"""gnomad-metrics drift probe: HEAD against the public constraint tables.

The constraint tables live in the ``gcp-public-data--gnomad`` release bucket over
plain HTTPS with no auth, so every table the builder can consume is directly
addressable. Issue #177 rated this source "marginally feasible (fragile,
network-only)"; in practice GCS returns a stronger fingerprint than the hgnc
reference does, because for a simple (non-composite) object the ETag is an MD5
over the body rather than an mtime-derived validator.

All declared tables are probed rather than only the default, so the fingerprint
covers whichever release a given build pinned. Each is a HEAD, so nothing is
transferred (the v4.0 table alone is 86 MB), and all of them share one
``requests.Session`` so the three objects cost one TLS handshake.

Compared surface: the MD5 ETag and Content-Length per object, plus
``x-goog-generation`` -- a GCS counter that changes on every object rewrite even
if the bytes happen to be identical, which makes a silent republish visible.
These live under ``headers`` rather than ``checksums`` because ``_conventions``
§ 12 defines ``checksums`` as "sha256 over the bytes used to derive ``headers``";
this probe fetches no body and computes no such digest, so recording a raw
validator there would misrepresent the contract. ``Last-Modified`` is demoted to
``informational`` following the hgnc precedent, where 8 of 8 drift PRs moved only
the timestamp.

Note that ``headers`` is a drift-bot schema key, so any change here is tiered
"schema" rather than "routine". That is deliberate and follows
``classify_risk``'s own stated policy -- "Defaults to schema for anything it
cannot read. Misclassifying a schema change as routine would bury it in a batch;
the reverse just opens one extra PR." A HEAD-only probe cannot see column names,
so the conservative tier is the correct one.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError, normalize_etag
from hvantk.core.utils.http import request_with_retry
from hvantk.skills.gnomad_metrics.shared.constants import (
    GNOMAD_CONSTRAINT_TABLES,
    GNOMAD_RELEASE_BASE_URL,
)

PROBE_VERSION = 2

# Evaluated once: GNOMAD_CONSTRAINT_TABLES is a module constant and cannot change
# at runtime.
OBJECT_PATHS: tuple[str, ...] = tuple(
    sorted(
        path
        for tables in GNOMAD_CONSTRAINT_TABLES.values()
        for path in tables.values()
    )
)

# Sized against the runner's budget, not per request. drift_cli defaults
# --timeout to 60s and enforces it with a single SIGALRM around the whole probe,
# while requests applies its timeout separately to connect and read. At the
# previous 30s this loop's worst case was 3 x 60 = 180s, so a merely slow bucket
# reported probe_failed on a healthy source. A (connect, read) pair keeps the
# whole loop inside the alarm.
_TIMEOUT_S = (5.0, 10.0)

# A server that compresses on the fly omits Content-Length and mangles the ETag.
_HEADERS = {"Accept-Encoding": "identity"}


def _head_object(
    session: requests.Session, path: str
) -> tuple[dict, str | None]:
    """Return (compared signals, Last-Modified) for one object, or raise."""
    url = f"{GNOMAD_RELEASE_BASE_URL}/{path}"
    try:
        head = request_with_retry(
            "HEAD",
            url,
            session=session,
            timeout=_TIMEOUT_S,
            allow_redirects=True,
            headers=_HEADERS,
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure for {path}: {exc}") from exc

    content_length = head.headers.get("Content-Length")
    etag = normalize_etag(head.headers.get("ETag"))

    # Fail closed on emptiness, and require BOTH validators. Accepting one alone
    # was the hole: an ETag-only response recorded a null length, and a
    # Content-Length-only response silently dropped the object's key out of
    # `headers` entirely -- after which the bot regenerates, commits the degraded
    # baseline, and the MD5 signal for that table is gone for good while drift
    # keeps reporting clean.
    if not content_length or not etag:
        raise DriftProbeError(
            f"gnomAD response for {path} carried no usable content signal "
            f"(Content-Length={content_length!r}, ETag={etag!r}); refusing to "
            "record a fingerprint that would compare equal to its baseline "
            "forever."
        )

    # Asserted request-side above; verified response-side here, because a caching
    # proxy may compress regardless and the recorded length would then describe
    # the compressed body rather than the object.
    encoding = (head.headers.get("Content-Encoding") or "identity").lower()
    if encoding != "identity":
        raise DriftProbeError(
            f"gnomAD response for {path} was {encoding}-encoded despite an "
            "identity request; Content-Length would describe the compressed body."
        )

    compared = {
        "content_length": content_length,
        "etag": etag,
        "generation": head.headers.get("x-goog-generation"),
    }
    return compared, head.headers.get("Last-Modified")


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live gnomAD constraint tables (HEAD only)."""
    headers: dict[str, dict[str, str | None]] = {}
    informational: dict[str, str | None] = {}
    failures: list[str] = []

    with requests.Session() as session:
        for path in OBJECT_PATHS:
            try:
                compared, last_modified = _head_object(session, path)
            except DriftProbeError as exc:
                # Collected rather than raised immediately: raising on the first
                # failure meant a fault on either frozen-since-2020 v2.1.1 object
                # aborted before v4.0 -- the current release -- was ever probed,
                # so a genuine v4.0 republish stayed invisible for as long as the
                # unrelated fault persisted.
                failures.append(str(exc))
            else:
                headers[path] = compared
                informational[path] = last_modified

    if failures:
        raise DriftProbeError(
            f"{len(failures)} of {len(OBJECT_PATHS)} gnomAD objects could not be "
            "fingerprinted: " + "; ".join(failures)
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": headers,
        "checksums": {},
        "informational": informational,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
