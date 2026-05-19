"""Plugin discovery and registration for hvantk providers.

Scans hvantk/skills/<provider>/plugin.yaml and the hvantk.providers entry-point
group, validates each manifest against the schema, resolves the declared
builder/probe callables, and constructs Provider records.

The module-level REGISTRY is built lazily on first access via get_registry().
"""

from __future__ import annotations

import importlib
import json
import logging
from importlib.metadata import entry_points
from pathlib import Path
from typing import Any, Callable, Iterable

import jsonschema
import yaml

from .api import (
    DatasetSpec,
    PluginLoadError,
    PluginNameCollision,
    Provider,
    TestPaths,
)

logger = logging.getLogger(__name__)

_SCHEMA_PATH = Path(__file__).parent / "manifest.schema.json"
_SKILLS_ROOT = Path(__file__).resolve().parent.parent.parent / "skills"
_ENTRY_POINT_GROUP = "hvantk.providers"


def _load_schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


class PluginRegistry:
    """Registry of all loaded provider plugins.

    Construction is empty; populate via load_from_directory,
    load_from_skills_root, and/or load_from_entry_points.
    """

    def __init__(self) -> None:
        self._providers: dict[str, Provider] = {}
        self._datasets: dict[str, DatasetSpec] = {}
        self._load_errors: list[tuple[str, Exception]] = []
        self._loaded_dirs: set[Path] = set()
        self._schema = _load_schema()

    # --- Public lookup API ---

    def get_provider(self, name: str) -> Provider:
        return self._providers[name]

    def get_dataset(self, name: str) -> DatasetSpec:
        return self._datasets[name]

    def list_providers(self) -> list[Provider]:
        return list(self._providers.values())

    def list_datasets(
        self, *, domain: str | None = None, backend: str | None = None
    ) -> list[DatasetSpec]:
        out = list(self._datasets.values())
        if domain is not None:
            out = [d for d in out if d.domain == domain]
        if backend is not None:
            out = [d for d in out if d.backend == backend]
        return out

    def load_errors(self) -> list[tuple[str, Exception]]:
        return list(self._load_errors)

    # --- Loading entry points ---

    def load_from_skills_root(self, root: Path | None = None) -> None:
        """Scan hvantk/skills/<provider>/plugin.yaml (or the given root)."""
        root = root or _SKILLS_ROOT
        if not root.is_dir():
            return
        for child in sorted(root.iterdir()):
            if child.is_dir() and not child.name.startswith("_"):
                if (child / "plugin.yaml").is_file():
                    self.load_from_directory(child)

    def load_from_directory(self, plugin_dir: Path) -> None:
        """Load a single plugin from its directory.

        Idempotent: loading the same resolved directory twice is a no-op. This
        prevents double-registration (and a spurious PluginNameCollision) when
        a plugin is discovered via both load_from_skills_root and
        load_from_entry_points after `poetry install` exposes the entry point.
        """
        plugin_dir = Path(plugin_dir).resolve()
        if plugin_dir in self._loaded_dirs:
            return
        plugin_id = str(plugin_dir)
        try:
            manifest = self._read_and_validate_manifest(plugin_dir / "plugin.yaml")
            provider = self._build_provider(manifest, plugin_dir)
            self._register(provider, plugin_id)
            self._loaded_dirs.add(plugin_dir)
        except PluginNameCollision:
            # Hard error: silent shadowing is the worst failure mode.
            raise
        except PluginLoadError as exc:
            self._load_errors.append((plugin_id, exc))
        except Exception as exc:  # noqa: BLE001
            err = PluginLoadError(str(exc))
            err.__cause__ = exc
            self._load_errors.append((plugin_id, err))

    def load_from_entry_points(self) -> None:
        """Iterate hvantk.providers entry points and load each."""
        try:
            eps: Iterable = entry_points(group=_ENTRY_POINT_GROUP)
        except TypeError:
            eps = entry_points().get(_ENTRY_POINT_GROUP, [])
        for ep in eps:
            try:
                module = ep.load()
                module_path = Path(module.__file__).resolve().parent
                self.load_from_directory(module_path)
            except PluginNameCollision:
                raise
            except Exception as exc:  # noqa: BLE001
                err = PluginLoadError(str(exc))
                err.__cause__ = exc
                self._load_errors.append((f"entry-point:{ep.name}", err))

    # --- Internal helpers ---

    def _read_and_validate_manifest(self, manifest_path: Path) -> dict:
        if not manifest_path.is_file():
            raise PluginLoadError(f"missing plugin.yaml at {manifest_path}")
        try:
            content = yaml.safe_load(manifest_path.read_text())
        except yaml.YAMLError as exc:
            raise PluginLoadError(f"invalid YAML: {exc}") from exc
        try:
            jsonschema.validate(content, self._schema)
        except jsonschema.ValidationError as exc:
            raise PluginLoadError(f"schema validation failed: {exc.message}") from exc
        return content

    def _build_provider(self, manifest: dict, plugin_dir: Path) -> Provider:
        datasets: list[DatasetSpec] = []
        for ds_manifest in manifest["datasets"]:
            try:
                ds = self._build_dataset_spec(
                    provider_name=manifest["name"],
                    ds_manifest=ds_manifest,
                    plugin_dir=plugin_dir,
                )
                datasets.append(ds)
            except PluginLoadError as exc:
                self._load_errors.append(
                    (f"{manifest['name']}:{ds_manifest.get('name', '?')}", exc)
                )
        catalog_rel = manifest.get("catalog")
        catalog_path = (
            str((plugin_dir / catalog_rel).resolve()) if catalog_rel else None
        )
        # Derive the manifest's primary domain independently of builder
        # imports so catalog routing still works when a builder cannot
        # be imported (e.g. hail missing in a non-hail dev env).
        manifest_domains: dict[str, int] = {}
        for ds in manifest.get("datasets", []):
            d = ds.get("domain")
            if d:
                manifest_domains[d] = manifest_domains.get(d, 0) + 1
        primary_domain = (
            max(manifest_domains, key=manifest_domains.get) if manifest_domains else None
        )
        return Provider(
            name=manifest["name"],
            version=manifest["version"],
            datasets=tuple(datasets),
            catalog_path=catalog_path,
            primary_domain=primary_domain,
        )

    def _build_dataset_spec(
        self, *, provider_name: str, ds_manifest: dict, plugin_dir: Path
    ) -> DatasetSpec:
        compound = f"{provider_name}:{ds_manifest['name']}"
        builder = self._resolve_callable(
            ds_manifest["builder"]["module"], ds_manifest["builder"]["function"]
        )
        probe = self._resolve_callable(
            ds_manifest["drift_probe"]["module"],
            ds_manifest["drift_probe"]["function"],
        )
        tests = ds_manifest["tests"]
        test_paths = TestPaths(
            command=tests["command"],
            fixture=str((plugin_dir / tests["fixture"]).resolve()),
            schema_snapshot=str((plugin_dir / tests["schema_snapshot"]).resolve()),
            row_snapshot=str((plugin_dir / tests["row_snapshot"]).resolve()),
            drift_fingerprint=str((plugin_dir / tests["drift_fingerprint"]).resolve()),
        )
        skill_path = str((plugin_dir / ds_manifest["skill"]).resolve())

        lifecycle = ds_manifest.get("lifecycle") or {}
        download_fn = None
        parse_fn = None
        if "download" in lifecycle:
            download_fn = self._resolve_callable(
                lifecycle["download"]["module"],
                lifecycle["download"]["function"],
            )
        if "parse" in lifecycle:
            parse_fn = self._resolve_callable(
                lifecycle["parse"]["module"],
                lifecycle["parse"]["function"],
            )

        return DatasetSpec(
            name=compound,
            domain=ds_manifest["domain"],
            backend=ds_manifest["backend"],
            builder=builder,
            drift_probe=probe,
            skill_path=skill_path,
            test_paths=test_paths,
            download_fn=download_fn,
            parse_fn=parse_fn,
        )

    def _resolve_callable(self, module_path: str, func_name: str) -> Callable[..., Any]:
        try:
            module = importlib.import_module(module_path)
        except ImportError as exc:
            raise PluginLoadError(
                f"cannot import module '{module_path}': {exc}"
            ) from exc
        try:
            func = getattr(module, func_name)
        except AttributeError as exc:
            raise PluginLoadError(
                f"module '{module_path}' has no function '{func_name}'"
            ) from exc
        if not callable(func):
            raise PluginLoadError(
                f"'{module_path}.{func_name}' is not callable"
            )
        return func

    def _register(self, provider: Provider, plugin_id: str) -> None:
        if provider.name in self._providers:
            raise PluginNameCollision(
                f"provider name collision: '{provider.name}' already registered "
                f"(new attempt from {plugin_id})"
            )
        self._providers[provider.name] = provider
        for ds in provider.datasets:
            if ds.name in self._datasets:
                # Defensive: with compound keys this is unreachable when the
                # provider-name check above passes. Kept as a backstop.
                raise PluginNameCollision(
                    f"dataset name collision: '{ds.name}'"
                )
            self._datasets[ds.name] = ds


# --- Module-level singleton (lazy) ---

_REGISTRY: PluginRegistry | None = None


def get_registry() -> PluginRegistry:
    """Return the module-level registry, building it on first access."""
    global _REGISTRY
    if _REGISTRY is None:
        reg = PluginRegistry()
        reg.load_from_skills_root()
        reg.load_from_entry_points()
        _REGISTRY = reg
    return _REGISTRY


def reset_registry_for_tests() -> None:
    """Test-only: drop the cached registry so the next get_registry() rebuilds it."""
    global _REGISTRY
    _REGISTRY = None
