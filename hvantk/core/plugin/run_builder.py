"""Phase B orchestrator: parse → drift probe → BuildContext → build_fn → save.

Skills register a build_fn with signature `(parsed_input, ctx, **params) -> Artifact`.
This orchestrator wraps that contract: it computes the source fingerprint via the
plugin's drift_probe, constructs a BuildContext, invokes the builder, validates the
returned artifact's type matches the plugin.yaml `artifact_type`, then persists via
`artifact.save(output_path)`.

Plugins that have NOT yet migrated to this contract continue to be invoked through
_TABLE_BUILDERS / _MATRIX_BUILDERS — those facades stay in place for Phase B coexistence.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from hvantk.core.models import BuildContext, Provenance
from hvantk.core.plugin.api import DatasetSpec


class BuilderContractError(RuntimeError):
    """Raised when a plugin's build_fn return value or schema_id violates the contract."""


def run_builder_for_spec(
    spec: DatasetSpec,
    *,
    parsed_input: Any,
    output_path: Path,
    plugin_version: str,
    builder_commit: str | None = None,
    **params,
) -> Provenance:
    """Run a plugin's build_fn under the Phase B contract.

    Args:
        spec: The DatasetSpec resolved from the plugin manifest. Must have
            artifact_type populated (Phase B-migrated plugins only).
        parsed_input: Whatever spec.parse_fn returned, or None if the plugin has none.
        output_path: Where to write the produced artifact.
        plugin_version: The plugin's declared version (from plugin.yaml top-level `version:`).
        builder_commit: Optional hvantk commit/version string for traceability.
        **params: Additional kwargs forwarded to the builder.

    Returns:
        The Provenance stamped onto the saved artifact.

    Raises:
        BuilderContractError: if the builder returns a non-Artifact, the wrong
            Artifact subclass, or stamps a schema_id different from the manifest.
    """
    if spec.artifact_type is None:
        raise BuilderContractError(
            f"{spec.name}: plugin manifest has no artifact_type; "
            f"migrate to Phase B contract before calling run_builder_for_spec()"
        )

    probe_result = spec.drift_probe()
    fingerprint = probe_result.get("fingerprint", "<no-fingerprint>")

    ctx = BuildContext(
        plugin=spec.name.split(":", 1)[0],
        dataset=spec.name,
        plugin_version=plugin_version,
        source_fingerprint=fingerprint,
        builder_commit=builder_commit,
    )

    artifact = spec.builder(parsed_input, ctx, **params)

    if not isinstance(artifact, spec.artifact_type):
        raise BuilderContractError(
            f"{spec.name}: build_fn returned {type(artifact).__name__}, "
            f"expected {spec.artifact_type.__name__}"
        )

    if spec.schema_id and artifact.provenance.schema_id != spec.schema_id:
        raise BuilderContractError(
            f"{spec.name}: build_fn stamped schema_id "
            f"{artifact.provenance.schema_id!r}, manifest declares {spec.schema_id!r}"
        )

    artifact.save(output_path)
    return artifact.provenance
