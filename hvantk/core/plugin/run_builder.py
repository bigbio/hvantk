"""Build orchestrator: parse → drift probe → BuildContext → build_fn → save.

Skills register a build_fn with signature `(parsed_input, ctx, **params) -> Artifact`.
This orchestrator wraps that contract: it computes the source fingerprint via the
plugin's drift_probe, constructs a BuildContext, invokes the builder, validates the
returned artifact's type matches the plugin.yaml `artifact_type`, then persists via
`artifact.save(output_path)`.

All in-tree plugins use this plugin builder contract; the orchestrator is the sole
dispatch path for `hvantk reprocess` and any other consumer of plugin builders.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
from pathlib import Path
from typing import Any

from hvantk.core.models import BuildContext, Provenance
from hvantk.core.plugin.api import (
    DatasetSpec,
    DriftProbeError,
    PROBE_FINGERPRINT_IGNORED_KEYS,
    PROBE_UNAVAILABLE_TOKEN,
)

logger = logging.getLogger(__name__)

#: The directory of the hvantk package this process runs: a checkout, or the
#: container image's own copy.
_PACKAGE_DIR = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
#: containers/hvantk.def writes this into the image's own copy of the package; a
#: checkout never has one (it is gitignored).
_BUILD_STAMP = "_build_commit.txt"


def _usable_commit(value: str | None) -> str | None:
    """``value`` stripped, or None when it is empty, the def file's "unknown"
    default, or a ``{{ ... }}`` template the container tool never substituted."""
    value = (value or "").strip()
    if not value or value == "unknown" or "{{" in value:
        return None
    return value


def resolve_builder_commit() -> str:
    """Identify the hvantk code this process runs, for a dataset's builder_commit.

    The first usable answer wins:

    1. git, when the running package directory is hvantk's *own* checkout --
       confirmed by requiring the repository's toplevel to be the package
       directory's parent, so a wheel install dropped into some unrelated
       project's .venv does not record that project's HEAD. "-dirty" is appended
       when the package's tracked files have uncommitted changes, since a bare
       commit then misdescribes what ran.
    2. The build stamp in the running package directory, which
       containers/hvantk.def writes into the image's own copy: the commit the
       image was built from. Only the running package's own stamp counts, so a
       checkout run inside the image (through PYTHONPATH) never records the
       image's commit.
    3. The HVANTK_BUILDER_COMMIT environment variable. The image has no git, so
       the hpc-hail-job templates pass the commit of the checkout they run.
    4. ``hvantk <version>``.

    Any failure falls through to the next source: a traceability nice-to-have
    must never break a build.
    """
    package_dir = _PACKAGE_DIR
    try:
        toplevel = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=package_dir,
            capture_output=True,
            text=True,
            timeout=5,
        )
        expected_root = os.path.realpath(os.path.dirname(package_dir))
        actual_root = os.path.realpath(toplevel.stdout.strip())
        if toplevel.returncode == 0 and actual_root == expected_root:
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=package_dir,
                capture_output=True,
                text=True,
                timeout=5,
            )
            commit = result.stdout.strip()
            if result.returncode == 0 and commit:
                dirty = subprocess.run(
                    [
                        "git",
                        "status",
                        "--porcelain",
                        "--untracked-files=no",
                        "--",
                        package_dir,
                    ],
                    cwd=package_dir,
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                if dirty.returncode == 0 and dirty.stdout.strip():
                    commit += "-dirty"
                return commit
    except Exception:
        pass

    try:
        with open(os.path.join(package_dir, _BUILD_STAMP)) as fh:
            stamped = _usable_commit(fh.read())
    except OSError:
        stamped = None
    if stamped:
        return stamped

    env_commit = _usable_commit(os.environ.get("HVANTK_BUILDER_COMMIT"))
    if env_commit:
        return env_commit

    from importlib.metadata import PackageNotFoundError, version

    try:
        return f"hvantk {version('hvantk')}"
    except PackageNotFoundError:
        return "hvantk <unknown>"


class BuilderContractError(RuntimeError):
    """Raised when a plugin's build_fn return value or schema_id violates the contract."""


def _coerce_fingerprint(probe_result: Any, dataset_name: str) -> str:
    """Build a stable fingerprint string from a drift probe's return dict.

    Probes return a dict like:
        {"probe_version": ..., "source_version": ..., "headers": {...},
         "checksums": {...}, "fetched_at": ...}

    The canonical fingerprint is a sha256 of the probe dict, EXCLUDING
    ``fetched_at`` (which is timestamp noise). If the probe explicitly
    returns a top-level ``fingerprint`` key, that wins (lets a probe
    publish its own canonical fingerprint).

    Raises BuilderContractError if the probe returned something other than a
    dict.
    """
    if probe_result is None:
        raise BuilderContractError(
            f"{dataset_name}: drift probe returned None; expected dict"
        )
    if not isinstance(probe_result, dict):
        raise BuilderContractError(
            f"{dataset_name}: drift probe returned {type(probe_result).__name__}, "
            f"expected dict"
        )
    if "fingerprint" in probe_result:
        return str(probe_result["fingerprint"])

    canonical = {
        k: v for k, v in probe_result.items() if k not in PROBE_FINGERPRINT_IGNORED_KEYS
    }
    payload = json.dumps(canonical, sort_keys=True, default=str)
    return f"sha256:{hashlib.sha256(payload.encode()).hexdigest()}"


def run_builder_for_spec(
    spec: DatasetSpec,
    *,
    parsed_input: Any,
    output_path: Path,
    plugin_version: str,
    builder_commit: str | None = None,
    **params,
) -> Provenance:
    """Run a plugin's build_fn under the builder contract.

    Args:
        spec: The DatasetSpec resolved from the plugin manifest. Must have
            artifact_type populated (plugins on the builder contract only).
        parsed_input: Whatever spec.parse_fn returned, or None if the plugin has none.
        output_path: Where to write the produced artifact.
        plugin_version: The plugin's declared version (from plugin.yaml top-level `version:`).
        builder_commit: The hvantk commit or version to record; None (the default)
            records what :func:`resolve_builder_commit` finds.
        **params: Additional kwargs forwarded to the builder.

    Returns:
        The Provenance stamped onto the saved artifact.

    Raises:
        BuilderContractError: if the builder returns a non-Artifact, the wrong
            Artifact subclass, or stamps a schema_id not accepted by the manifest.
    """
    if spec.artifact_type is None:
        raise BuilderContractError(
            f"{spec.name}: plugin manifest has no artifact_type; "
            f"migrate to the plugin builder contract before calling run_builder_for_spec()"
        )

    # A probe failure must not destroy the build. The probe supplies provenance
    # metadata, not build input, and every manual-acquisition plugin is built
    # from a file staged by hand -- often on a compute node with no egress. Before
    # the documentation-only plugins gained live probes their probes were pure
    # in-process calls, so those builds needed no network at all; letting a
    # DriftProbeError propagate here would have made `hvantk reprocess
    # --skip-download --no-check-drift` fail offline, and `--no-check-drift` gates
    # only the post-build check, not this call.
    #
    # The fallback is a self-describing token rather than a synthesized digest, so
    # the provenance record never implies a probe ran (same reasoning as
    # STUB_FINGERPRINT_TOKEN).
    try:
        probe_result = spec.drift_probe()
        fingerprint = _coerce_fingerprint(probe_result, spec.name)
    except DriftProbeError as exc:
        logger.warning(
            "%s: drift probe could not reach its source (%s); stamping provenance "
            "with %r. The artifact is built normally; run `hvantk drift %s` once "
            "connectivity is available to record a real fingerprint.",
            spec.name,
            exc,
            PROBE_UNAVAILABLE_TOKEN,
            spec.name,
        )
        fingerprint = PROBE_UNAVAILABLE_TOKEN

    ctx = BuildContext(
        plugin=spec.name.split(":", 1)[0],
        dataset=spec.name,
        plugin_version=plugin_version,
        source_fingerprint=fingerprint,
        builder_commit=(
            builder_commit if builder_commit is not None else resolve_builder_commit()
        ),
    )

    artifact = spec.builder(parsed_input, ctx, **params)

    if not isinstance(artifact, spec.artifact_type):
        raise BuilderContractError(
            f"{spec.name}: build_fn returned {type(artifact).__name__}, "
            f"expected {spec.artifact_type.__name__}"
        )

    schema_ids = getattr(spec, "schema_ids", ())
    if not isinstance(schema_ids, (tuple, list, set, frozenset)):
        schema_ids = ()
    accepted_schema_ids = {spec.schema_id, *schema_ids} - {None}
    if accepted_schema_ids and artifact.provenance.schema_id not in accepted_schema_ids:
        raise BuilderContractError(
            f"{spec.name}: build_fn stamped schema_id "
            f"{artifact.provenance.schema_id!r}, manifest accepts "
            f"{sorted(accepted_schema_ids)!r}"
        )

    artifact.save(output_path)
    return artifact.provenance
