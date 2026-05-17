"""ClinVar drift probe: HEAD probe against the latest GRCh38 VCF.

ClinVar publishes the current monthly release at a stable URL
(``clinvar.vcf.gz``) under the NCBI FTP service. A HEAD request returns
``Last-Modified`` and ``Content-Length`` headers, which form a lightweight
fingerprint that changes once per monthly release.

We deliberately do not stream the VCF body itself (it is ~500 MB) -- the
schema is parsed from the VCF header inside the builder, and snapshot-level
schema drift is caught by ``tests/snapshots/schema.json`` regenerated against
the fixture, not by this probe.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin_api import DriftProbeError
from hvantk.core.constants import CLINVAR_FTP_BASE

PROBE_VERSION = 1
_FILENAME = "clinvar.vcf.gz"
_URL = f"{CLINVAR_FTP_BASE}/{_FILENAME}"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live ClinVar GRCh38 VCF.

    Issues a single HEAD request and captures ``Last-Modified`` and
    ``Content-Length`` headers. Does not download the VCF body.
    """
    try:
        resp = requests.head(_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    last_modified = resp.headers.get("Last-Modified")
    content_length = resp.headers.get("Content-Length")

    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {_FILENAME: {"content_length": content_length}},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
