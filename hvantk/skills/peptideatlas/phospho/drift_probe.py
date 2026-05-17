"""PeptideAtlas drift probe: HEAD probe against the pinned phospho build.

PeptideAtlas publishes versioned phospho-proteome builds under URLs of the
shape::

    {PEPTIDEATLAS_PHOSPHO_BASE_URL}/{build_date}/atlas_build_{build_id}.tsv.zip

A new build is released once or twice per year. The pinned ``(build_date,
build_id)`` lives in :mod:`hvantk.ptm.constants` and is exported by the
dataset class as ``PeptideAtlasPhosphoDataset.from_latest()``. A HEAD request
against that URL returns ``Last-Modified`` and ``Content-Length``, which
together form a lightweight fingerprint that flips when a new build is cut
or the existing archive is regenerated.

We deliberately do not download the (multi-gigabyte) zip body - schema drift
inside the archive is caught by the builder snapshot tests, not this probe.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin_api import DriftProbeError
from hvantk.skills.peptideatlas.phospho.shared.datasets import PeptideAtlasPhosphoDataset

PROBE_VERSION = 1
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the pinned PeptideAtlas phospho build."""
    dataset = PeptideAtlasPhosphoDataset.from_latest()
    url = dataset.zip_url
    filename = url.rsplit("/", 1)[-1]

    try:
        resp = requests.head(url, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    last_modified = resp.headers.get("Last-Modified")
    content_length = resp.headers.get("Content-Length")

    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {
            filename: {
                "content_length": content_length,
                "build_date": dataset.build_date,
                "build_id": dataset.build_id,
            }
        },
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
