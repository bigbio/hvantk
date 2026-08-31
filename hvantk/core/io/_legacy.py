"""Legacy raw-file shim.

When `core/io.load(path)` finds an artifact file with no sidecar manifest,
it falls back to wrapping the file with `Provenance.unknown(reason=...)`.
This shim is the bridge that lets Phase A consumers read pre-Phase-B
artifact files (raw .parquet, .h5ad, .ht) without rebuilding everything.
"""
from __future__ import annotations

from pathlib import Path

from hvantk.core.models.provenance import Provenance


def unknown_provenance_for(path: Path) -> Provenance:
    return Provenance.unknown(reason=f"legacy file, no manifest: {path.name}")
