"""gnomad-metrics drift probe: HEAD against the public constraint tables.

The constraint tables live in the ``gcp-public-data--gnomad`` release bucket over
plain HTTPS with no auth, so every table the builder can consume is directly
addressable. Issue #177 rated this source "marginally feasible (fragile,
network-only)"; in practice GCS returns a stronger fingerprint than the hgnc
reference does, because for a simple (non-composite) object the ETag is an MD5
over the body rather than an mtime-derived validator.

All three declared tables are probed rather than only the default, so the
fingerprint covers whichever release a given build pinned. Each is a HEAD, so
nothing is transferred (the v4.0 table alone is 86 MB).

Compared surface: the MD5 ETag and Content-Length per object, plus
``x-goog-generation`` -- a GCS counter that changes on every object rewrite even
if the bytes happen to be identical, which makes a silent republish visible.
``Last-Modified`` is demoted to ``informational`` following the hgnc precedent,
where 8 of 8 drift PRs moved only a timestamp.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.gnomad_metrics.shared.constants import (
    GNOMAD_CONSTRAINT_TABLES,
    GNOMAD_RELEASE_BASE_URL,
)

PROBE_VERSION = 1
_TIMEOUT_S = 30

# A server that compresses on the fly omits Content-Length and mangles the ETag.
# GCS does not today, but the fingerprint must not silently change meaning if
# that ever turns on.
_HEADERS = {"Accept-Encoding": "identity"}


def _object_paths() -> list[str]:
    """Every declared constraint object path, in a stable order."""
    return sorted(
        path
        for tables in GNOMAD_CONSTRAINT_TABLES.values()
        for path in tables.values()
    )


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live gnomAD constraint tables (HEAD only)."""
    checksums: dict[str, str] = {}
    extras: dict[str, dict[str, str | None]] = {}
    informational: dict[str, dict[str, str | None]] = {}

    for path in _object_paths():
        url = f"{GNOMAD_RELEASE_BASE_URL}/{path}"
        try:
            head = requests.head(
                url, timeout=_TIMEOUT_S, allow_redirects=True, headers=_HEADERS
            )
            head.raise_for_status()
        except requests.RequestException as exc:
            raise DriftProbeError(f"HTTP failure for {path}: {exc}") from exc

        etag = head.headers.get("ETag")
        content_length = head.headers.get("Content-Length")
        # Fail closed. Without a content signal the fingerprint would compare
        # equal to its baseline forever -- the false-green the stub avoided.
        if etag is None and content_length is None:
            raise DriftProbeError(
                f"gnomAD response for {path} carried neither ETag nor "
                "Content-Length; refusing to record a fingerprint with no "
                "content signal."
            )

        if etag is not None:
            checksums[path] = etag.strip('"')
        extras[path] = {
            "content_length": content_length,
            "generation": head.headers.get("x-goog-generation"),
        }
        informational[path] = {"last_modified": head.headers.get("Last-Modified")}

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {},
        "checksums": checksums,
        "extras": extras,
        "informational": informational,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
