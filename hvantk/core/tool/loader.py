"""Discover hvantk tool manifests under hvantk/tools/**/<basename>.tool.yaml.

Tool manifests are purely metadata for agent/human discoverability; they do
not drive command registration (that still happens via the existing
``cli.add_command(...)`` wiring in ``hvantk/hvantk.py``). The manifest is
OPTIONAL: tools without a ``<basename>.tool.yaml`` continue to work, they
just do not show up in ``hvantk tools list``.

The module-level registry is built lazily on first access via
``get_registry()``; tests reset it via ``reset_registry_for_tests()``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import jsonschema
import yaml

from .api import (
    Subcommand,
    ToolLoadError,
    ToolRequirements,
    ToolSpec,
)

logger = logging.getLogger(__name__)

_SCHEMA_PATH = Path(__file__).parent / "manifest.schema.json"
_TOOLS_ROOT = Path(__file__).resolve().parent.parent.parent / "tools"


def _load_schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


class ToolRegistry:
    """Registry of all loaded tool manifests.

    Construction is empty; populate via ``load_from_tools_root``.
    """

    def __init__(self) -> None:
        self._tools: dict[str, ToolSpec] = {}
        self._load_errors: list[tuple[str, Exception]] = []
        self._schema = _load_schema()

    # --- Public lookup API ---

    def get_tool(self, name: str) -> ToolSpec:
        return self._tools[name]

    def list_tools(self, *, domain: str | None = None) -> list[ToolSpec]:
        out = list(self._tools.values())
        if domain is not None:
            out = [t for t in out if t.domain == domain]
        return out

    def load_errors(self) -> list[tuple[str, Exception]]:
        return list(self._load_errors)

    # --- Loading entry points ---

    def load_from_tools_root(self, root: Path | None = None) -> None:
        """Scan ``hvantk/tools/**/<basename>.tool.yaml`` (or the given root)."""
        root = root or _TOOLS_ROOT
        if not root.is_dir():
            return
        for path in sorted(root.glob("**/*.tool.yaml")):
            self._load_manifest(path)

    # --- Internal helpers ---

    def _load_manifest(self, manifest_path: Path) -> None:
        try:
            content = yaml.safe_load(manifest_path.read_text())
            jsonschema.validate(content, self._schema)
            spec = self._build_spec(content, manifest_path)
        except (yaml.YAMLError, jsonschema.ValidationError) as exc:
            err = ToolLoadError(f"{manifest_path}: {exc}")
            err.__cause__ = exc
            self._load_errors.append((str(manifest_path), err))
            return
        except Exception as exc:  # noqa: BLE001
            err = ToolLoadError(f"{manifest_path}: {exc}")
            err.__cause__ = exc
            self._load_errors.append((str(manifest_path), err))
            return

        if spec.name in self._tools:
            err = ToolLoadError(
                f"tool name collision: '{spec.name}' already registered "
                f"(new attempt from {manifest_path})"
            )
            self._load_errors.append((str(manifest_path), err))
            return
        self._tools[spec.name] = spec

    def _build_spec(self, manifest: dict, manifest_path: Path) -> ToolSpec:
        subs = tuple(
            Subcommand(
                name=s["name"],
                description=s["description"],
                inputs=tuple(s.get("inputs", [])),
                outputs=tuple(s.get("outputs", [])),
            )
            for s in manifest.get("subcommands", [])
        )
        reqs_dict = manifest.get("requires", {}) or {}
        reqs = ToolRequirements(
            hail=reqs_dict.get("hail", False),
            network=reqs_dict.get("network", False),
            extras=tuple(reqs_dict.get("extras", [])),
        )
        return ToolSpec(
            name=manifest["name"],
            domain=manifest["domain"],
            type=manifest["type"],
            description=manifest["description"],
            cli_module=manifest["cli"]["module"],
            cli_callable=manifest["cli"]["function"],
            purpose_short=manifest["purpose"]["short"],
            purpose_long=manifest["purpose"].get("long", ""),
            subcommands=subs,
            requires=reqs,
            manifest_path=str(manifest_path.resolve()),
        )


# --- Module-level singleton (lazy) ---

_REGISTRY: ToolRegistry | None = None


def get_registry() -> ToolRegistry:
    """Return the module-level registry, building it on first access."""
    global _REGISTRY
    if _REGISTRY is None:
        reg = ToolRegistry()
        reg.load_from_tools_root()
        _REGISTRY = reg
    return _REGISTRY


def reset_registry_for_tests() -> None:
    """Test-only: drop the cached registry so the next get_registry() rebuilds it."""
    global _REGISTRY
    _REGISTRY = None
