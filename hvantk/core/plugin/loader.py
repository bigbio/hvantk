"""Plugin discovery and registration for hvantk providers.

Scans hvantk/skills/<provider>/plugin.yaml and the hvantk.providers entry-point
group, validates each manifest against the schema, resolves the declared
builder/probe callables, and constructs Provider records.

Discovery is split into two passes:

  Pass 1 (descriptive, eager) — reads YAML, validates against the schema, and
  populates ``DatasetManifest`` objects.  No callables are imported.  The result
  is always available via ``registry.list_manifests()``.

  Pass 2 (executable, lazy) — on first ``get_dataset(name)`` call, resolves the
  callables declared in the manifest via importlib and caches the resulting
  ``DatasetSpec``.  If an import fails the manifest remains visible but
  ``get_dataset`` raises.

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
    DatasetManifest,
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
        # Pass-1 index: all manifests (descriptive only, no callables).
        self._manifests: dict[str, DatasetManifest] = {}
        # Pass-2 cache: successfully-resolved executable specs.
        self._datasets: dict[str, DatasetSpec] = {}
        self._load_errors: list[tuple[str, Exception]] = []
        self._loaded_dirs: set[Path] = set()
        self._schema = _load_schema()
        # CLI entries per plugin, keyed by plugin name.
        self._cli_entries: dict[str, list[dict]] = {}

    # --- Public lookup API ---

    def get_provider(self, name: str) -> Provider:
        return self._providers[name]

    def get_dataset(self, name: str) -> DatasetSpec:
        """Return the executable spec, resolving callables lazily on first call."""
        if name in self._datasets:
            return self._datasets[name]
        if name not in self._manifests:
            raise KeyError(name)
        spec = self._resolve_spec(self._manifests[name])
        self._datasets[name] = spec
        return spec

    def list_providers(self) -> list[Provider]:
        return list(self._providers.values())

    def list_datasets(
        self, *, domain: str | None = None, backend: str | None = None
    ) -> list[DatasetSpec]:
        """Return successfully-bound executable specs (lazy resolution per entry).

        Only specs whose callables can be imported are returned.  Use
        ``list_manifests()`` to get a full descriptive view that survives
        missing optional runtimes.
        """
        out = []
        for name, dm in self._manifests.items():
            try:
                spec = self.get_dataset(name)
            except Exception:  # noqa: BLE001
                continue
            out.append(spec)
        if domain is not None:
            out = [d for d in out if d.domain == domain]
        if backend is not None:
            out = [d for d in out if d.backend == backend]
        return out

    def list_manifests(self) -> list[DatasetManifest]:
        """Return ALL dataset manifests (descriptive), regardless of bind status."""
        return list(self._manifests.values())

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

    # --- CLI wiring ---

    def apply_plugin_downloaders(self, click_group: Any) -> None:
        """Wire each plugin's downloader CLI into the given click group.

        Reads the ``cli:`` block from each plugin manifest (cached at load
        time).  For entries whose ``command`` ends in ``-download``, strips
        the suffix to derive the subcommand name under the given group,
        lazily resolves the function, and calls
        ``click_group.add_command(fn, name=short_name)``.

        Entries that don't end in ``-download`` are skipped — they're
        reserved for future top-level commands.

        Idempotent: skips entries whose short name is already registered.
        """
        for plugin_name, entries in self._cli_entries.items():
            for entry in entries:
                command = entry["command"]
                if not command.endswith("-download"):
                    continue
                short_name = command[: -len("-download")]
                if short_name in click_group.commands:
                    continue
                try:
                    fn = self._resolve_callable(entry["module"], entry["function"])
                except PluginLoadError as exc:
                    logger.warning(
                        "skipping downloader %r from plugin %r: %s",
                        short_name,
                        plugin_name,
                        exc,
                    )
                    continue
                click_group.add_command(fn, name=short_name)

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

    def _build_dataset_manifest(
        self,
        *,
        provider_name: str,
        ds_manifest: dict,
        plugin_dir: Path,
        plugin_version: str | None = None,
    ) -> DatasetManifest:
        """Pass 1: build descriptive manifest — no callable imports."""
        compound = f"{provider_name}:{ds_manifest['name']}"
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
        return DatasetManifest(
            name=compound,
            domain=ds_manifest["domain"],
            backend=ds_manifest["backend"],
            skill_path=skill_path,
            test_paths=test_paths,
            plugin_name=provider_name,
            plugin_version=plugin_version,
            artifact_type_name=ds_manifest.get("artifact_type"),
            schema_id=ds_manifest.get("schema_id"),
            builder_ref=(
                ds_manifest["builder"]["module"],
                ds_manifest["builder"]["function"],
            ),
            drift_probe_ref=(
                ds_manifest["drift_probe"]["module"],
                ds_manifest["drift_probe"]["function"],
            ),
            download_ref=(
                (lifecycle["download"]["module"], lifecycle["download"]["function"])
                if "download" in lifecycle
                else None
            ),
            parse_ref=(
                (lifecycle["parse"]["module"], lifecycle["parse"]["function"])
                if "parse" in lifecycle
                else None
            ),
            has_download_fn="download" in lifecycle,
            has_parse_fn="parse" in lifecycle,
        )

    def _resolve_spec(self, dm: DatasetManifest) -> DatasetSpec:
        """Pass 2: import callables and return an executable DatasetSpec."""
        builder = self._resolve_callable(*dm.builder_ref)
        drift_probe = self._resolve_callable(*dm.drift_probe_ref)
        download_fn = (
            self._resolve_callable(*dm.download_ref) if dm.download_ref else None
        )
        parse_fn = (
            self._resolve_callable(*dm.parse_ref) if dm.parse_ref else None
        )
        artifact_type = None
        if dm.artifact_type_name:
            from hvantk.core import models as _models

            artifact_type = getattr(_models, dm.artifact_type_name, None)
            if artifact_type is None:
                raise PluginLoadError(
                    f"{dm.name}: artifact_type {dm.artifact_type_name!r} not found "
                    f"in hvantk.core.models"
                )
        return DatasetSpec(
            name=dm.name,
            domain=dm.domain,
            backend=dm.backend,
            builder=builder,
            drift_probe=drift_probe,
            skill_path=dm.skill_path,
            test_paths=dm.test_paths,
            download_fn=download_fn,
            parse_fn=parse_fn,
            plugin_version=dm.plugin_version,
            artifact_type=artifact_type,
            schema_id=dm.schema_id,
        )

    def _build_provider(self, manifest: dict, plugin_dir: Path) -> Provider:
        """Build a Provider via two-pass discovery.

        Pass 1 (always runs): build DatasetManifest for every declared dataset
        and store in self._manifests.  No callable imports.

        Pass 2 (attempted eagerly here, cached on demand): try to resolve
        specs for the Provider.datasets tuple so callers that never call
        get_dataset() still get populated Provider objects.  Failures are
        recorded in _load_errors; the manifest entry always survives.
        """
        dm_list: list[DatasetManifest] = []
        for ds_manifest in manifest["datasets"]:
            try:
                dm = self._build_dataset_manifest(
                    provider_name=manifest["name"],
                    ds_manifest=ds_manifest,
                    plugin_dir=plugin_dir,
                    plugin_version=manifest.get("version"),
                )
                dm_list.append(dm)
                self._manifests[dm.name] = dm
            except PluginLoadError as exc:
                self._load_errors.append(
                    (f"{manifest['name']}:{ds_manifest.get('name', '?')}", exc)
                )

        # Eagerly attempt to bind specs (best-effort; failures are soft).
        datasets: list[DatasetSpec] = []
        for dm in dm_list:
            try:
                spec = self._resolve_spec(dm)
                self._datasets[dm.name] = spec
                datasets.append(spec)
            except PluginLoadError as exc:
                self._load_errors.append((dm.name, exc))
            except Exception as exc:  # noqa: BLE001
                err = PluginLoadError(str(exc))
                err.__cause__ = exc
                self._load_errors.append((dm.name, err))

        # Cache CLI entries for downloader wiring.
        cli_entries = manifest.get("cli", [])
        if cli_entries:
            self._cli_entries[manifest["name"]] = cli_entries

        catalog_rel = manifest.get("catalog")
        catalog_path = (
            str((plugin_dir / catalog_rel).resolve()) if catalog_rel else None
        )
        # Derive primary domain from YAML (not from resolved specs) so that
        # catalog routing works even when builder imports fail.
        manifest_domains: dict[str, int] = {}
        for ds in manifest.get("datasets", []):
            d = ds.get("domain")
            if d:
                manifest_domains[d] = manifest_domains.get(d, 0) + 1
        primary_domain = (
            max(manifest_domains, key=manifest_domains.get)  # type: ignore[arg-type]
            if manifest_domains
            else None
        )
        return Provider(
            name=manifest["name"],
            version=manifest["version"],
            datasets=tuple(datasets),
            catalog_path=catalog_path,
            primary_domain=primary_domain,
            manifests=tuple(dm_list),
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
        # Register manifests (all datasets, bind or not).
        for dm in provider.manifests:
            if dm.name in self._manifests:
                # Defensive: already set during _build_provider; this is a no-op.
                pass
        # Register successfully-bound specs.
        for ds in provider.datasets:
            if ds.name in self._datasets:
                # Defensive: already cached; no action needed.
                pass


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
