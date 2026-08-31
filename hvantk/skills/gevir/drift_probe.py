"""GeVIR drift probe: HEAD against the published supplementary table.

GeVIR metrics are distributed as supplementary data to the Nature Genetics paper
(PMID 31873297, DOI 10.1038/s41588-019-0560-2), served from Springer's
static-content CDN at a stable, direct URL. Issue #177 recorded this source as
having "no probeable URL at all (publication PMID 31873297 only)" and shipped a
stub sentinel instead. That conflated *publication-only distribution* with *no
addressable URL*: the ESM object below answers a HEAD with both Content-Length
and a content-hash ETag, so an hgnc-style probe is feasible after all.

The authors' code repository (github.com/gevirank/gevir), which the skill's
SKILL.md cites as the distribution point, is NOT a usable source: it ships the
analysis code, and its ``tables/`` directory holds only a placeholder file.

Of the six MOESM slots for this article only MOESM3 is public (the rest answer
403); it carries sheets ``table_2`` and ``table_6`` with the
``gene_id`` / ``gevir_percentile`` / ``loeuf_percentile`` / ``virlof_percentile``
columns the builder consumes.

Compared surface: the ETag, which Springer serves as an MD5 over the object body
and is therefore a true content digest, plus Content-Length. ``Last-Modified`` is
recorded under ``informational`` so a byte-identical republish cannot open a
pull request, following the hgnc precedent where 8 of 8 drift PRs moved only the
timestamp.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1

# Supplementary Tables workbook for Abramovs et al., Nat Genet 2020.
GEVIR_SUPPLEMENTARY_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41588-019-0560-2/MediaObjects/41588_2019_560_MOESM3_ESM.xlsx"
)

_FILENAME = "41588_2019_560_MOESM3_ESM.xlsx"
_TIMEOUT_S = 30

# Pinned so Content-Length always describes the stored object. A server that
# compresses on the fly omits Content-Length and mangles the ETag (the INSIDER
# probe hit exactly that); Springer does not today, but the fingerprint should
# not silently change meaning if that ever turns on.
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
    etag = head.headers.get("ETag")
    last_modified = head.headers.get("Last-Modified")

    # Fail closed. With neither validator the fingerprint would carry no content
    # signal at all and would compare equal to its baseline forever -- the exact
    # false-green the stub was put in place to avoid. probe_failed is loud and
    # recoverable; a contentless "clean" is neither.
    if content_length is None and etag is None:
        raise DriftProbeError(
            "GeVIR response carried neither Content-Length nor ETag; refusing to "
            "record a fingerprint with no content signal."
        )

    checksums = {}
    if etag is not None:
        # Springer serves an MD5 over the body, quoted per RFC 7232.
        checksums[_FILENAME] = etag.strip('"')

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {},
        "checksums": checksums,
        "extras": {"content_length": content_length},
        "informational": {"last_modified": last_modified},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
