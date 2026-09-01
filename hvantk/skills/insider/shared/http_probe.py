"""Shared HEAD-fingerprint helper for the two INSIDER products.

``insider:variants`` and ``insider:interfaces`` read two different files from one
portal, so they share the fetch mechanics but must NOT share a fingerprint: they
are versioned independently upstream (the BED has not moved since 2018-03-05
while the interfaces table moved 2024-05-15), and a shared baseline makes the
drift bot record a ledger row only for the group's anchor dataset, so an
interfaces-only change can never be reported against the dataset that actually
needs rebuilding. Each dataset therefore calls this with its own file and commits
its own baseline.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError, normalize_etag

PROBE_VERSION = 2

INSIDER_BASE_URL = "http://interactomeinsider.yulab.org"

# The portal serves no HTTPS listener (a TLS connection to the same host fails),
# so the probe is forced onto cleartext. That is a real exposure: an interposing
# proxy can define what this records. It is mitigated, not solved, by comparing
# two independent validators rather than one -- an interstitial would have to
# forge a consistent (Content-Length, ETag) pair to pass unnoticed -- and by the
# fact that any such forgery reads as drift rather than as clean.
_TIMEOUT_S = 30

# The portal compresses text/plain on the fly when the client offers it, and a
# compressed response omits Content-Length entirely while appending "-gzip" to
# the ETag. requests sends "gzip, deflate" by default, so the probe must opt out
# or the interfaces file's signals disappear.
_HEADERS = {"Accept-Encoding": "identity"}


def head_fingerprint(filename: str, url: str) -> dict:
    """Fingerprint one INSIDER product by HEAD, transferring no body.

    Compared surface: Content-Length **and** the ETag, both under ``headers``.

    Two corrections over the first version of this probe. The ETag is not
    redundant with Content-Length: decoding the live values shows the tag is
    ``hex(size)-hex(mtime)`` (``0x2f752a5`` is exactly the interfaces file's
    49,762,981-byte length), so it moves on size *or* mtime while Content-Length
    moves only on size. Demoting it therefore made every equal-size edit -- a
    swapped accession, a corrected residue index -- undetectable. The hgnc
    precedent for demoting validators does not transfer: hgnc republishes
    byte-identical content weekly, whereas these two objects are static archives
    that have not moved in years, so the no-op-PR risk is near-nil and the signal
    given up was the only one that catches an equal-size change.

    And the signals live under ``headers`` rather than ``extras`` because the
    drift bot reads ``headers``/``checksums`` as its schema signal and tiers
    anything else as "routine". With both of those empty, a re-release that
    changed the ``track name=`` format -- the case insider/SKILL.md § 8 says
    breaks the parser -- would have been swept into a batch whose body tells the
    reviewer the schema signal is unchanged. peptideatlas ships this same shape:
    validator metadata under ``headers``, ``checksums`` left empty.
    """
    try:
        head = requests.head(
            url, timeout=_TIMEOUT_S, allow_redirects=True, headers=_HEADERS
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure for {filename}: {exc}") from exc

    content_length = head.headers.get("Content-Length")
    etag = normalize_etag(head.headers.get("ETag"))

    # Fail closed on emptiness, not just absence: a proxy answering `ETag: ""`
    # would otherwise record an empty digest, which placeholder_baseline_reason
    # later reads as a hand-seeded baseline -- trapping the dataset in a
    # probe_failed loop that regenerating cannot clear.
    if not content_length or not etag:
        raise DriftProbeError(
            f"INSIDER response for {filename} carried no usable content signal "
            f"(Content-Length={content_length!r}, ETag={etag!r}); refusing to "
            "record a fingerprint that would compare equal to its baseline "
            "forever."
        )

    # Asserted request-side above; verified response-side here, because a
    # caching proxy or gzip_static may compress regardless and the recorded
    # length would then describe the compressed body rather than the object.
    encoding = (head.headers.get("Content-Encoding") or "identity").lower()
    if encoding != "identity":
        raise DriftProbeError(
            f"INSIDER response for {filename} was {encoding}-encoded despite an "
            "identity request; Content-Length would describe the compressed body."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {filename: {"content_length": content_length, "etag": etag}},
        "checksums": {},
        "informational": {
            filename: {"last_modified": head.headers.get("Last-Modified")}
        },
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
