"""UCSC Cell Browser drift probe.

The UCSC Cell Browser ships data per collection, with each collection
versioned independently. There is no global "UCSC Cell Browser version"
string. To produce a stable drift signal we fingerprint two things:

* The local pinned-collection catalog
  (``hvantk/skills/ucsc_cellbrowser/data/cells_ucsc_datasets.json``) --
  this changes only when we deliberately add/remove tracked collections,
  and is the source of truth for reproducible builds.
* The upstream root ``https://cells.ucsc.edu/dataset.json`` HEAD --
  ``Last-Modified``/``ETag`` flip whenever ANY collection in the
  browser is updated. Surfaced as ``upstream_version`` and
  ``upstream_checksum`` so consumers can see when something changed
  upstream even if our local pin is unchanged.

This is a "catalog freshness" probe: it does not detect column drift on
a per-collection basis (each collection has its own ``cellbrowser.conf``
+ expression matrix). Per-collection probing should be done by the
builder when materializing each collection.
"""

from __future__ import annotations

import hashlib
import importlib.resources
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
UCSC_DATASET_INDEX_URL = "https://cells.ucsc.edu/dataset.json"
_LOCAL_CATALOG_PACKAGE = "hvantk.skills.ucsc_cellbrowser"
_LOCAL_CATALOG_SUBDIR = "data"
_LOCAL_CATALOG_NAME = "cells_ucsc_datasets.json"
_TIMEOUT_S = 30


def _local_catalog_bytes() -> bytes:
    """Read the pinned UCSC Cell Browser catalog shipped with hvantk."""
    return (
        importlib.resources.files(_LOCAL_CATALOG_PACKAGE)
        .joinpath(_LOCAL_CATALOG_SUBDIR, _LOCAL_CATALOG_NAME)
        .read_bytes()
    )


def fetch_fingerprint() -> dict:
    """Fingerprint of the pinned UCSC catalog + upstream dataset index HEAD."""
    local_bytes = _local_catalog_bytes()
    local_checksum = hashlib.sha256(local_bytes).hexdigest()

    try:
        head = requests.head(
            UCSC_DATASET_INDEX_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    upstream_last_modified = head.headers.get("Last-Modified")
    upstream_etag = head.headers.get("ETag", "")
    upstream_content_length = head.headers.get("Content-Length", "")
    upstream_blob = (
        f"{upstream_etag}|{upstream_content_length}".encode("utf-8")
    )
    upstream_checksum = hashlib.sha256(upstream_blob).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": upstream_last_modified,
        "headers": {
            _LOCAL_CATALOG_NAME: ["local-catalog-hash"],
            "dataset.json": ["ETag", "Content-Length", "Last-Modified"],
        },
        "checksums": {
            _LOCAL_CATALOG_NAME: local_checksum,
            "dataset.json": upstream_checksum,
        },
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
