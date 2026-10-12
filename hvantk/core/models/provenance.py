"""Provenance: required metadata attached to every artifact instance.

Every artifact on disk or in memory carries a Provenance so drift detection,
caching, and reproducibility can work end-to-end. Use Provenance.unknown(reason)
only for tests, the legacy-file shim, or scripts that legitimately have no
upstream lineage; production code must plumb real provenance.

``build_parameters`` is optional: a builder may use it to record the options that
shaped the artifact's contents (AlphaGenome records its ``output_types`` and
``ontology_curies`` filters). ``{}`` and ``None`` mean different things, and both
are distinct from "has entries":

  * a dict with entries -- the builder ran with those options;
  * ``{}`` -- the builder recorded that it ran with NO options;
  * ``None`` -- nothing was recorded at all. This is what a sidecar written before
    this field existed reads back as (``core.io._manifest`` reads a missing key as
    None, never defaulting it to ``{}``) -- otherwise a filtered AlphaGenome
    artifact built before this field existed would look identical to one explicitly
    built with no filters.

Every Provenance a builder stamps today gets a real dict, never None -- only a
value read back from a legacy sidecar can be None, so callers that read an
EXISTING artifact's ``build_parameters`` (rather than one they just stamped) must
handle it, e.g. by treating None the same as "unknown" rather than iterating over
it. The value must be JSON-serialisable (``BuildContext.provenance`` stores a JSON
copy and rejects anything else).
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
    # None only ever comes from reading a legacy sidecar (see the class docstring);
    # a freshly stamped Provenance always gets a real dict here.
    build_parameters: dict[str, object] | None = field(default_factory=dict, hash=False)

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
