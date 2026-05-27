"""Shared drift-probe helpers for CPTAC datasets.

CPTAC data is fetched via the ``cptac`` Python package, not via direct
HTTP. A meaningful drift signal therefore needs to flip when the
upstream package publishes a new release. We probe two signals:

* The installed ``cptac`` package version (``importlib.metadata``). A
  bump here means the local environment is on a different schema.
* The latest GitHub release of ``PayneLab/cptac``
  (``api.github.com/repos/PayneLab/cptac/releases/latest``). A bump
  here means upstream has shipped a new release that the local
  install has not yet picked up.

Both signals are combined into a single fingerprint, parameterized by
the dataset name (``protein_expression`` / ``phosphoproteomics``) so
the two CPTAC plugins can share the implementation while still
emitting per-dataset filenames in the ``headers`` / ``checksums``
dicts.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from typing import Optional

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
GITHUB_LATEST_RELEASE_URL = (
    "https://api.github.com/repos/PayneLab/cptac/releases/latest"
)
_TIMEOUT_S = 30


def cptac_version() -> Optional[str]:
    """Return the installed ``cptac`` package version, or ``None``."""
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8
        return None
    try:
        return version("cptac")
    except PackageNotFoundError:
        return None


def _fetch_latest_release() -> dict:
    """Hit the GitHub releases-latest API for ``PayneLab/cptac``."""
    try:
        resp = requests.get(
            GITHUB_LATEST_RELEASE_URL,
            timeout=_TIMEOUT_S,
            headers={"Accept": "application/vnd.github+json"},
            allow_redirects=True,
        )
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc
    except ValueError as exc:
        raise DriftProbeError(f"GitHub release JSON parse failure: {exc}") from exc


def fetch_fingerprint(dataset_name: str) -> dict:
    """Lightweight drift fingerprint for a CPTAC dataset.

    Parameters
    ----------
    dataset_name : str
        Short identifier for the CPTAC dataset (e.g. ``protein_expression``).
        Used to key the per-dataset entries in ``headers`` and ``checksums``.
    """
    pkg_version = cptac_version()
    if pkg_version is None:
        raise DriftProbeError(
            "The 'cptac' Python package is not installed; cannot fingerprint. "
            "Install it with: pip install cptac"
        )

    release = _fetch_latest_release()
    tag_name = release.get("tag_name", "")
    published_at = release.get("published_at", "")

    fingerprint_blob = f"{pkg_version}|{tag_name}|{published_at}".encode("utf-8")
    checksum = hashlib.sha256(fingerprint_blob).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": tag_name or pkg_version,
        "headers": {
            dataset_name: [
                "installed_cptac_version",
                "latest_release_tag",
                "latest_release_published_at",
            ],
        },
        "checksums": {dataset_name: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "extras": {
            "installed_cptac_version": pkg_version,
            "latest_release_tag": tag_name,
            "latest_release_published_at": published_at,
        },
    }
