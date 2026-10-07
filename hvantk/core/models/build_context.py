"""BuildContext: the platform's gift to a skill's build_fn.

Skill authors supply schema_id when calling ctx.provenance(...). They may also
pass build_parameters, an optional JSON-serialisable dict of the options that
shaped the artifact; few builders do (see Provenance). Everything else (plugin
name, version, source fingerprint, builder commit) is platform-computed and
plumbed through this dataclass.
"""

from __future__ import annotations

import json
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
        """Stamp this build's provenance, storing a JSON copy of build_parameters.

        The copy makes a non-JSON value (a Path, NaN) raise here, before save()
        writes any data; detaches the record from the caller's dict; and stores a
        tuple as the list the sidecar reloads, so a reloaded provenance equals
        the stamped one.
        """
        build_parameters = json.loads(
            json.dumps(build_parameters or {}, sort_keys=True, allow_nan=False)
        )
        return Provenance(
            plugin=self.plugin,
            dataset=self.dataset,
            plugin_version=self.plugin_version,
            source_fingerprint=self.source_fingerprint,
            schema_id=schema_id,
            build_timestamp=datetime.now(timezone.utc),
            builder_commit=self.builder_commit,
            build_parameters=build_parameters,
        )
