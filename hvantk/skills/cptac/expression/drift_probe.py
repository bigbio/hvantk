"""CPTAC protein expression drift probe.

CPTAC protein-expression data is fetched via the ``cptac`` Python package
(`pip install cptac`), not direct HTTP. The installed ``cptac`` package
version is the de-facto release pin -- a bump indicates new bundled
datasets or upstream-data updates fetched by the package at runtime.

This probe shares the package-version fingerprint pattern with the
sibling :mod:`hvantk.skills.cptac.phospho.drift_probe`.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from hvantk.core.plugin_api import DriftProbeError

PROBE_VERSION = 1


def _cptac_version() -> Optional[str]:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8
        return None
    try:
        return version("cptac")
    except PackageNotFoundError:
        return None


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the installed CPTAC package version."""
    pkg_version = _cptac_version()
    if pkg_version is None:
        raise DriftProbeError(
            "The 'cptac' Python package is not installed; cannot fingerprint. "
            "Install it with: pip install cptac"
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": pkg_version,
        "headers": {
            "cptac": {
                "package_version": pkg_version,
            }
        },
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
