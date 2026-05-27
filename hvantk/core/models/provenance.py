"""Provenance: required metadata attached to every artifact instance.

Every artifact on disk or in memory carries a Provenance so drift detection,
caching, and reproducibility can work end-to-end. Use Provenance.unknown(reason)
only for tests, the legacy-file shim, or scripts that legitimately have no
upstream lineage; production code must plumb real provenance.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone


@dataclass(frozen=True)
class Provenance:
    plugin: str
    dataset: str
    plugin_version: str
    source_fingerprint: str
    schema_id: str
    build_timestamp: datetime
    builder_commit: str | None
    parents: tuple["Provenance", ...] = field(default_factory=tuple)

    @classmethod
    def unknown(cls, *, reason: str) -> "Provenance":
        return cls(
            plugin="<unknown>",
            dataset="<unknown>",
            plugin_version="<unknown>",
            source_fingerprint=f"<unknown: {reason}>",
            schema_id="<unknown>",
            build_timestamp=datetime.now(timezone.utc),
            builder_commit=None,
        )
