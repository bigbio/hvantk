"""BuildContext: the platform's gift to a skill's build_fn.

The platform supplies plugin name, version, source fingerprint and builder
commit; builders may also record content-affecting parameters with provenance.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone

from hvantk.core.models.provenance import Provenance


@dataclass(frozen=True)
class BuildContext:
    plugin: str
    dataset: str
    plugin_version: str
    source_fingerprint: str
    builder_commit: str | None

    def provenance(
        self, *, schema_id: str, build_parameters: dict[str, object] | None = None
    ) -> Provenance:
        return Provenance(
            plugin=self.plugin,
            dataset=self.dataset,
            plugin_version=self.plugin_version,
            source_fingerprint=self.source_fingerprint,
            schema_id=schema_id,
            build_timestamp=datetime.now(timezone.utc),
            builder_commit=self.builder_commit,
            build_parameters=build_parameters or {},
        )
