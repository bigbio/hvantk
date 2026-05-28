"""Drift probes for `onek-genomes:variants` and `onek-genomes:samples`.

Both probes are release-identity probes (the upstream sources are immutable
per release). A new release ships at a new URL — bumping the constants here
is how a new release flips the fingerprint.
"""
from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1

# --- Variants probe constants ---

_VARIANTS_RELEASE_DIR = "1000G_2504_high_coverage"
_VARIANTS_ANCHOR_FILE = (
    "CCDG_14151_B01_GRM_WGS_2020-08-05_chr22"
    ".filtered.shapeit2-duohmm-phased.vcf.gz"
)
_VARIANTS_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections"
_VARIANTS_URL = (
    f"{_VARIANTS_BASE}/{_VARIANTS_RELEASE_DIR}"
    f"/working/20201028_3202_phased/{_VARIANTS_ANCHOR_FILE}"
)

# --- Samples probe constants ---

_SAMPLES_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
_SAMPLES_VERSION = "integrated_call_samples_v3.20130502"

_TIMEOUT_S = 30


def fetch_variants_fingerprint() -> dict:
    """Release identity for the 1KG high-coverage callset.

    HEADs chr22 (the smallest autosome in the release). Last-Modified +
    Content-Length form the per-file anchor; release directory name forms
    the human-readable source_version.
    """
    try:
        resp = requests.head(_VARIANTS_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure on {_VARIANTS_URL}: {exc}") from exc

    return {
        "probe_version": PROBE_VERSION,
        "source_version": _VARIANTS_RELEASE_DIR,
        "headers": {_VARIANTS_ANCHOR_FILE: {
            "last_modified": resp.headers.get("Last-Modified"),
            "content_length": resp.headers.get("Content-Length"),
        }},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }


def fetch_samples_fingerprint() -> dict:
    """Release identity for the 1KG canonical samples panel."""
    try:
        resp = requests.head(_SAMPLES_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure on {_SAMPLES_URL}: {exc}") from exc

    return {
        "probe_version": PROBE_VERSION,
        "source_version": _SAMPLES_VERSION,
        "headers": {"integrated_call_samples_v3.20130502.ALL.panel": {
            "last_modified": resp.headers.get("Last-Modified"),
            "content_length": resp.headers.get("Content-Length"),
        }},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
