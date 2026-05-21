"""Sidecar .provenance.json manifest read/write.

Manifest path is the artifact path with `.provenance.json` appended:
  rows.parquet                -> rows.parquet.provenance.json
  expr.h5ad                   -> expr.h5ad.provenance.json
  brca.geneset.json           -> brca.geneset.json.provenance.json
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from hvantk.core.models.provenance import Provenance


def manifest_path(artifact_path: Path) -> Path:
    return artifact_path.with_name(artifact_path.name + ".provenance.json")


def write_manifest(provenance: Provenance, artifact_path: Path) -> None:
    payload = _to_dict(provenance)
    manifest_path(artifact_path).write_text(json.dumps(payload, indent=2))


def read_manifest(artifact_path: Path) -> Provenance | None:
    mp = manifest_path(artifact_path)
    if not mp.exists():
        return None
    try:
        payload = json.loads(mp.read_text())
        return _from_dict(payload)
    except (json.JSONDecodeError, KeyError, ValueError) as e:
        import logging
        logging.getLogger(__name__).warning(
            "corrupt provenance manifest at %s (%s); falling back to legacy shim",
            mp, e,
        )
        return None


def _to_dict(p: Provenance) -> dict:
    return {
        "plugin": p.plugin,
        "dataset": p.dataset,
        "plugin_version": p.plugin_version,
        "source_fingerprint": p.source_fingerprint,
        "schema_id": p.schema_id,
        "build_timestamp": p.build_timestamp.isoformat(),
        "builder_commit": p.builder_commit,
        "parents": [_to_dict(parent) for parent in p.parents],
    }


def _from_dict(d: dict) -> Provenance:
    return Provenance(
        plugin=d["plugin"],
        dataset=d["dataset"],
        plugin_version=d["plugin_version"],
        source_fingerprint=d["source_fingerprint"],
        schema_id=d["schema_id"],
        build_timestamp=datetime.fromisoformat(d["build_timestamp"]),
        builder_commit=d.get("builder_commit"),
        parents=tuple(_from_dict(p) for p in d.get("parents", [])),
    )
