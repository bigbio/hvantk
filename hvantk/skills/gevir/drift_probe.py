"""GeVIR drift probe: HEAD against the published supplementary object.

GeVIR metrics are distributed as supplementary data to the Nature Genetics paper
(PMID 31873297, DOI 10.1038/s41588-019-0560-2), served from Springer's
static-content CDN at a stable, direct URL. Issue #177 recorded this source as
having "no probeable URL at all (publication PMID 31873297 only)" and shipped a
stub sentinel; that conflated *publication-only distribution* with *no
addressable URL*. The ESM object below answers a HEAD with both Content-Length
and a content-hash ETag, so an hgnc-style probe is feasible after all.

The authors' code repository (github.com/gevirank/gevir) is NOT a usable source:
it ships the analysis code, and its ``tables/`` directory holds only a
placeholder file.

Of the six MOESM slots for this article only MOESM3 is public (the rest answer
403). It is an Excel workbook, not the bgzipped TSV the builder consumes -- sheet
``table_2`` has to be extracted and converted first -- so this probe watches the
*published upstream*, while the local build input is materialized externally.

Compared surface: the ETag and Content-Length, both under ``headers``. Springer
serves the ETag as an MD5 over the object body, so unlike a size+mtime validator
it is a true content digest. They live under ``headers`` rather than ``checksums``
because the drift bot reads ``checksums`` as a hash of the column-header row --
a schema signal -- and this probe never fetches a body, so it has no schema
signal to offer; recording one there tiered every routine content update as a
schema change. peptideatlas ships this same shape: validator metadata under
``headers``, ``checksums`` left empty.

``Last-Modified`` is demoted to ``informational`` following the hgnc precedent,
where 8 of 8 drift PRs moved only the timestamp.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError, normalize_etag

PROBE_VERSION = 2

# Supplementary Tables workbook for Abramovs et al., Nat Genet 2020.
GEVIR_SUPPLEMENTARY_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41588-019-0560-2/MediaObjects/41588_2019_560_MOESM3_ESM.xlsx"
)

_FILENAME = "41588_2019_560_MOESM3_ESM.xlsx"
_TIMEOUT_S = 30

# Pinned so Content-Length always describes the stored object. A server that
# compresses on the fly omits Content-Length and appends a transform suffix to
# the ETag; Springer does not today, but the fingerprint must not silently change
# meaning if that ever turns on.
_HEADERS = {"Accept-Encoding": "identity"}


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live GeVIR supplementary object.

    Issues a single HEAD; the 10 MB workbook is never transferred.
    """
    try:
        head = requests.head(
            GEVIR_SUPPLEMENTARY_URL,
            timeout=_TIMEOUT_S,
            allow_redirects=True,
            headers=_HEADERS,
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    content_length = head.headers.get("Content-Length")
    etag = normalize_etag(head.headers.get("ETag"))

    # Fail closed on emptiness, not just absence. Requiring BOTH validators is
    # deliberate: accepting an ETag-only response would record
    # `content_length: None`, and because the drift bot regenerates drifted
    # baselines automatically, that null bakes in and every later normal response
    # reads as drift. An empty `ETag: ""` is likewise rejected rather than stored,
    # since placeholder_baseline_reason would then read the baseline as
    # hand-seeded and trap the dataset in a probe_failed loop that regenerating
    # cannot clear.
    if not content_length or not etag:
        raise DriftProbeError(
            f"GeVIR response carried no usable content signal "
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
            f"GeVIR response was {encoding}-encoded despite an identity request; "
            "Content-Length would describe the compressed body."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {_FILENAME: {"content_length": content_length, "etag": etag}},
        "checksums": {},
        "informational": {"last_modified": head.headers.get("Last-Modified")},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
