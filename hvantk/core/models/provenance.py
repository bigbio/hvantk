"""Provenance: required metadata attached to every artifact instance.

Every artifact on disk or in memory carries a Provenance so drift detection,
caching, and reproducibility can work end-to-end. Use Provenance.unknown(reason)
only for tests, the legacy-file shim, or scripts that legitimately have no
upstream lineage; production code must plumb real provenance.

``build_parameters`` is optional: a builder may use it to record the options that
shaped the artifact's contents (AlphaGenome records its ``output_types`` and
``ontology_curies`` filters). Most builders record nothing, so an empty dict does
not mean the build took no options. The value must be JSON-serialisable
(``BuildContext.provenance`` stores a JSON copy and rejects anything else), and it
is empty in sidecars written before the field existed.
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
    # A dict is unhashable, so it stays out of __hash__; __eq__ still compares it.
    build_parameters: dict[str, object] = field(default_factory=dict, hash=False)

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
