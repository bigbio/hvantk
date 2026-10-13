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
import re
from importlib.metadata import entry_points
from pathlib import Path
from typing import Any, Callable, Iterable, get_args

import jsonschema
import yaml

from .api import (
    Acquisition,
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

#: Same pattern manifest.schema.json enforces on a provider `name`. Kept here so the
#: error path can judge a name WITHOUT running the validation that just failed.
_PROVIDER_NAME_RE = re.compile(r"^[a-z][a-z0-9-]*$")

#: Load-error unit prefix for a manifest `cli:` entry whose callable will not resolve
#: (`cli:<command>`). Such a row is about the command alone -- a dataset whose own
#: callables fail with it is recorded separately, under `provider:dataset` -- so
#: `hvantk drift` does not count it.
CLI_UNIT_PREFIX = "cli:"


def is_cli_unit(unit: str) -> bool:
    """True for a load-error unit recorded for a manifest `cli:` entry."""
    return unit.startswith(CLI_UNIT_PREFIX)


def _load_schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


def _provider_id_hint(plugin_dir: Path) -> str:
    """Name a failed plugin the way its datasets are named, not by filesystem path.

    A load error is reported to humans and to machines: `drift --all --json` turns it
    into a row whose `dataset_name` goes verbatim into the "Drift probes failing" issue
    body, and `drift_to_pr.py` splits that value on ``:`` to look up maintainers. An
    absolute path leaks the build machine's layout into a public issue
    (``/home/runner/work/hvantk/...``) and gives the maintainer lookup a directory where
    it expects a provider.

    The manifest's own ``name`` is the right label, but the common reason for landing
    here is that the manifest would not parse -- so read it best-effort, and accept it
    only if it is shaped like a provider id. A manifest broken enough to reach this
    function may well be broken *in its name* (the `broken-manifest` fixture declares
    ``BROKEN UPPERCASE``), and echoing an invalid name into an issue body and a
    maintainer lookup helps nobody.

    The fallback is the directory's basename, which is a real identifier even though it
    is not always the provider name -- ``skills/gwas_catalog/`` declares ``gwas-catalog``
    -- so the manifest value is preferred when it is usable.

    Never raises: this runs on the error path, where a second failure would mask the
    first.
    """
    try:
        content = yaml.safe_load((plugin_dir / "plugin.yaml").read_text())
        name = content.get("name") if isinstance(content, dict) else None
        if isinstance(name, str) and _PROVIDER_NAME_RE.fullmatch(name):
            return name
    except Exception:  # noqa: BLE001 - best effort; the fallback is always usable
        pass
    return plugin_dir.name


def _acquisition_of(ds_manifest: dict, lifecycle: dict, compound: str) -> Acquisition:
    """Build the Acquisition record, rejecting the one incoherent combination.

    ``acquisition.mode: byo`` and ``lifecycle.download`` contradict each other: the
    first says the plugin cannot fetch its inputs, the second hands it a function that
    does. Whichever the runner honoured, the manifest would be lying about the other,
    so this is a load error rather than a precedence rule.

    Absent block means ``mode: download``, which keeps every existing manifest valid.
    """
    raw = ds_manifest.get("acquisition") or {}
    acq = Acquisition(
        mode=raw.get("mode", "download"),
        reason=raw.get("reason"),
        instructions=raw.get("instructions"),
    )
    if acq.is_byo and "download" in lifecycle:
        raise PluginLoadError(
            f"{compound}: acquisition.mode is 'byo' but lifecycle.download is declared; "
            "a dataset either fetches its own inputs or it does not."
        )
    return acq


def _dataset_label(content: Any, path: list) -> str | None:
    """The ``name:`` of the dataset at ``content["datasets"][path[1]]``, if ``path``
    (a jsonschema ``ValidationError.absolute_path``, already converted to a list)
    points there, else None.

    A manifest declares a LIST of datasets, so a schema failure inside one of them is
    reported by ``jsonschema`` only as an array index -- the ``PluginLoadError`` this
    raises is then filed under the PROVIDER's name (``_record_load_error`` has no
    dataset to attribute it to), with nothing saying which of a multi-dataset plugin's
    entries actually broke. Defensive: this augments an error message and must never
    raise one of its own, so any shape mismatch (a non-dict entry, an out-of-range
    index, a manifest that is not even a dict) falls through to None rather than
    masking the real ValidationError.
    """
    if len(path) < 2 or path[0] != "datasets" or not isinstance(path[1], int):
        return None
    try:
        dataset = content["datasets"][path[1]]
    except (LookupError, TypeError):
        return None
    name = dataset.get("name") if isinstance(dataset, dict) else None
    return name if isinstance(name, str) and name else f"#{path[1]}"


def _explain(exc: jsonschema.ValidationError, content: Any = None) -> str:
    """Render a schema error, translating the ones whose raw text explains nothing,
    and naming the dataset when the failing path points inside one.

    A conditional rule reports as `... should not be valid under {'required': [...]}`,
    which states the mechanism and not the problem. This is the message a plugin
    author sees -- at ``hvantk plugins validate`` time, and since #408 at load time
    too -- so the two ``if/then`` rules on ``acquisition``/``lifecycle`` (#360) get a
    sentence instead.

    ``content`` is the parsed manifest ``exc`` came from. Passing it lets a failure
    inside ``datasets[i]`` be attributed to that dataset's own ``name:`` instead of
    the bare provider, which is all a load error was recorded under before #408 --
    with no way to tell which entry of a multi-dataset plugin actually broke.
    ``hvantk plugins validate`` checks one manifest at a time and already names it on
    the command line, so it omits ``content`` and keeps today's wording.
    """
    path = list(exc.absolute_path)
    if exc.validator == "not" and "lifecycle" in path:
        message = (
            "acquisition.mode is 'byo' but lifecycle.download is declared. A dataset "
            "either fetches its own inputs or it cannot -- drop the downloader, or set "
            "acquisition.mode to 'download'."
        )
    elif exc.validator == "required" and (
        exc.validator_value == ["lifecycle"]
        or (exc.validator_value == ["download"] and path and path[-1] == "lifecycle")
    ):
        message = (
            "acquisition.mode is 'download' but no lifecycle.download is declared. "
            "Declare the downloader, or set acquisition.mode to 'byo' with a reason; "
            "omit the acquisition block while a downloader is still to be written."
        )
    else:
        where = " -> ".join(str(x) for x in path)
        message = f"{exc.message}{f' (at {where})' if where else ''}"
    dataset = _dataset_label(content, path) if content is not None else None
    return f"dataset {dataset!r}: {message}" if dataset is not None else message


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
        # Pass-2 negative cache: spec resolutions that failed eagerly.
        # Stored so get_dataset() returns a consistent error within a
        # session — without this, a transient resolution failure during
        # Pass 2 could appear to "self-heal" on a later get_dataset() call
        # (same registry, two different answers). Cleared on reset only.
        self._failed_datasets: dict[str, PluginLoadError] = {}
        self._load_errors: list[tuple[str, Exception]] = []
        self._loaded_dirs: set[Path] = set()
        self._schema = _load_schema()
        # CLI entries per plugin, keyed by plugin name.
        self._cli_entries: dict[str, list[dict]] = {}
        # Resolved `cli:` callables and cached failures, keyed by (plugin, command).
        # Filled on demand by `resolve_cli_entries` / `apply_plugin_downloaders`, never
        # at load time (see `resolve_cli_entries`).
        self._cli_callables: dict[tuple[str, str], Callable[..., Any]] = {}
        self._failed_cli: dict[tuple[str, str], PluginLoadError] = {}

    # --- Public lookup API ---

    def get_provider(self, name: str) -> Provider:
        return self._providers[name]

    def get_dataset(self, name: str) -> DatasetSpec:
        """Return the executable spec, resolving callables lazily on first call.

        If the dataset's spec failed to resolve during eager Pass 2, raises
        the cached :class:`PluginLoadError` rather than re-attempting the
        import. This keeps the registry deterministic for the lifetime of
        a session even when an underlying import is flaky.
        """
        if name in self._datasets:
            return self._datasets[name]
        if name in self._failed_datasets:
            raise self._failed_datasets[name]
        if name not in self._manifests:
            raise KeyError(name)
        try:
            spec = self._resolve_spec(self._manifests[name])
        except PluginLoadError as exc:
            self._failed_datasets[name] = exc
            self._record_load_error(name, exc)
            raise
        except Exception as exc:  # noqa: BLE001
            # Same wrapping as the eager pass in _build_provider: a module that raises
            # something other than ImportError at import time must produce the same cached
            # PluginLoadError on every lookup, not a different exception on the first.
            err = PluginLoadError(str(exc))
            err.__cause__ = exc
            self._failed_datasets[name] = err
            self._record_load_error(name, err)
            raise err from exc
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
            except PluginLoadError:
                # Already recorded by whichever path cached it -- the eager pass in
                # _build_provider or the lazy path in get_dataset -- so recording here
                # would duplicate the row. get_dataset raises nothing else for a name
                # taken from self._manifests.
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
                else:
                    # A manifest that fails to LOAD is recorded; one that fails to
                    # EXIST was not, so the whole provider left the registry with no
                    # row, no warning and exit 0 -- verbatim the #351 failure. The
                    # loader accepts only `plugin.yaml`, so a `.yml` typo, a bad
                    # rebase, or a packaging glob that stops shipping the manifest
                    # all land here. Underscore-prefixed directories are skipped
                    # above (`_conventions`, `_hooks`), so anything reaching this
                    # branch is shaped like a provider and claims to be one.
                    #
                    # Mark the directory as attempted too (#374): 13 of
                    # 21 in-tree providers are ALSO declared as entry points pointing
                    # at this same directory, and get_registry() runs
                    # load_from_entry_points() right after this method. Without this,
                    # that second pass calls load_from_directory() on the identical
                    # path, which fails again -- with a differently-worded message --
                    # and doubles this one broken unit's row.
                    self._loaded_dirs.add(child.resolve())
                    self._record_load_error(
                        child.name, PluginLoadError(f"no plugin.yaml in {child}")
                    )

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
        # Mark as attempted regardless of outcome (#374): this used to
        # happen only on success, so a directory whose manifest failed to load was
        # never marked, and get_registry()'s second pass (load_from_entry_points, which
        # resolves 13 of 21 in-tree providers back to this same directory) retried it
        # and recorded a second, differently-worded error for the same unit.
        self._loaded_dirs.add(plugin_dir)
        plugin_id = str(plugin_dir)
        try:
            manifest = self._read_and_validate_manifest(plugin_dir / "plugin.yaml")
            provider = self._build_provider(manifest, plugin_dir)
            self._register(provider, plugin_id)
        except PluginNameCollision:
            # Hard error: silent shadowing is the worst failure mode.
            raise
        except PluginLoadError as exc:
            self._record_load_error(_provider_id_hint(plugin_dir), exc)
        except Exception as exc:  # noqa: BLE001
            err = PluginLoadError(str(exc))
            err.__cause__ = exc
            self._record_load_error(_provider_id_hint(plugin_dir), err)

    def _record_load_error(self, unit: str, exc: Exception) -> None:
        """Collect a load failure AND say so.

        ``unit`` is whatever disappeared from the registry -- a plugin directory, an
        entry point, or a single ``provider:dataset`` key -- or a ``cli:<command>``
        entry whose callable did not resolve (`resolve_cli_entries`), which costs that
        command alone; a dataset that fails with it is recorded under its own key.

        Collecting silently is the problem this exists to stop. `plugins list` prints
        the survivors with no "N failed" line, and `drift-health.yml` greps the report
        for `probe_failed`, so anything that never registered emits no row at all. A
        dataset can go from monitored to unmonitored with every signal still green
        (#351).

        Every path that drops something from the registry must come through here.
        `get_dataset`'s lazy bind records on first failure and re-raises the cached
        error afterwards; `list_datasets` merely skips those cached failures, which is
        why it is the one `except` in this class that does not record. When
        this landed it covered only the two directory-load paths, leaving the entry-point
        path and all three dataset-binding paths silent -- and the dataset ones are the
        likelier failure, because they fire on a callable that will not import, which no
        amount of manifest validation can catch. (Manifest-level breakage takes the
        OTHER route: `jsonschema.validate` runs in `_read_and_validate_manifest`, so a
        malformed entry fails the whole document and is recorded at the provider level
        by `load_from_directory` -- one of the two paths already covered.)
        """
        self._load_errors.append((unit, exc))
        if is_cli_unit(unit):
            # "NOT registered and will not be drift-checked" would be false here (#408).
            logger.warning(
                "%r failed to load; that command is NOT available: %s", unit, exc
            )
            return
        logger.warning(
            "%r failed to load; it is NOT registered and will not be drift-checked: %s",
            unit,
            exc,
        )

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
                self._record_load_error(f"entry-point:{ep.name}", err)

    # --- CLI wiring ---

    def resolve_cli_entries(self) -> None:
        """Resolve every manifest ``cli:`` entry, recording each failure as a load error.

        A failure is recorded under ``cli:<command>``, not under the provider: the
        command is lost, while the provider's datasets bind (or fail) on their own, so
        `hvantk drift` ignores these rows. Before this, a downloader module that would
        not import was only logged by `apply_plugin_downloaders`, and its subcommand
        silently vanished from `hvantk download`.

        Called by `hvantk plugins errors`, so the check does not depend on `hvantk
        download` having run in the same process, and through `apply_plugin_downloaders`
        by the download CLI, which needs the callables anyway. Never at load time:
        a `cli:` module need not be one the datasets import, and it may pull heavy
        dependencies that every other registry consumer would then pay for.

        Idempotent: each entry is resolved once and its outcome cached, so calling this
        and `apply_plugin_downloaders` in one process records a failure once.
        """
        for plugin_name, entries in self._cli_entries.items():
            for entry in entries:
                self._resolve_cli_entry(plugin_name, entry)

    def _resolve_cli_entry(
        self, plugin_name: str, entry: dict
    ) -> Callable[..., Any] | None:
        """The callable of one ``cli:`` entry, or None once its failure is recorded."""
        key = (plugin_name, entry["command"])
        if key in self._cli_callables:
            return self._cli_callables[key]
        if key in self._failed_cli:
            return None
        try:
            fn = self._resolve_callable(entry["module"], entry["function"])
        except PluginLoadError as exc:
            err = exc
        except Exception as exc:  # noqa: BLE001
            # Same wrapping as the eager dataset pass: a module that raises something
            # other than ImportError at import time is a load failure too, not a crash
            # of whichever command happened to resolve it.
            err = PluginLoadError(str(exc))
            err.__cause__ = exc
        else:
            self._cli_callables[key] = fn
            return fn
        self._failed_cli[key] = err
        self._record_load_error(f"{CLI_UNIT_PREFIX}{entry['command']}", err)
        return None

    def apply_plugin_downloaders(self, click_group: Any) -> None:
        """Wire each plugin's downloader CLI into the given click group.

        Reads the ``cli:`` block from each plugin manifest (cached at load
        time).  For entries whose ``command`` ends in ``-download``, strips
        the suffix to derive the subcommand name under the given group,
        lazily resolves the function, and calls
        ``click_group.add_command(fn, name=short_name)``.

        Entries that don't end in ``-download`` are skipped — they're
        reserved for future top-level commands.

        An entry whose callable does not resolve is skipped, and recorded as a
        ``cli:<command>`` load error (see `resolve_cli_entries`).

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
                fn = self._resolve_cli_entry(plugin_name, entry)
                if fn is None:
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
            raise PluginLoadError(
                f"schema validation failed: {_explain(exc, content)}"
            ) from exc
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
            schema_ids=tuple(ds_manifest.get("schema_ids", ())),
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
            acquisition=_acquisition_of(ds_manifest, lifecycle, compound),
        )

    def _resolve_spec(self, dm: DatasetManifest) -> DatasetSpec:
        """Pass 2: import callables and return an executable DatasetSpec."""
        builder = self._resolve_callable(*dm.builder_ref)
        drift_probe = self._resolve_callable(*dm.drift_probe_ref)
        download_fn = (
            self._resolve_callable(*dm.download_ref) if dm.download_ref else None
        )
        parse_fn = self._resolve_callable(*dm.parse_ref) if dm.parse_ref else None
        artifact_type = None
        if dm.artifact_type_name:
            from hvantk.core import models as _models

            artifact_type = getattr(_models, dm.artifact_type_name, None)
            if artifact_type is None:
                raise PluginLoadError(
                    f"{dm.name}: artifact_type {dm.artifact_type_name!r} not found "
                    f"in hvantk.core.models"
                )
            # Must be one of the four artifact classes. An `isinstance(_, type)`
            # check is not enough: it rejects the `Artifact` union and module
            # objects, but still admits any other class the package exports
            # (`Provenance`, `BuildContext`), which then loads clean and fails
            # deep inside run_builder_for_spec instead of here. Checking
            # membership makes the guard say what its message promises.
            from hvantk.core.models.artifact import Artifact

            if artifact_type not in get_args(Artifact):
                allowed = ", ".join(t.__name__ for t in get_args(Artifact))
                raise PluginLoadError(
                    f"{dm.name}: artifact_type {dm.artifact_type_name!r} must name a "
                    f"concrete artifact class ({allowed}), not {artifact_type!r}"
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
            schema_ids=dm.schema_ids,
            acquisition=dm.acquisition,
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
                self._record_load_error(
                    f"{manifest['name']}:{ds_manifest.get('name', '?')}", exc
                )

        # Eagerly attempt to bind specs (best-effort; failures are soft).
        # Failed bindings are cached in ``_failed_datasets`` so that a
        # later get_dataset() call re-raises the same error rather than
        # silently re-attempting and possibly succeeding (which would
        # leave the same registry returning two different answers in one
        # session).
        datasets: list[DatasetSpec] = []
        for dm in dm_list:
            try:
                spec = self._resolve_spec(dm)
                self._datasets[dm.name] = spec
                datasets.append(spec)
            except PluginLoadError as exc:
                self._failed_datasets[dm.name] = exc
                self._record_load_error(dm.name, exc)
            except Exception as exc:  # noqa: BLE001
                err = PluginLoadError(str(exc))
                err.__cause__ = exc
                self._failed_datasets[dm.name] = err
                self._record_load_error(dm.name, err)

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
            raise PluginLoadError(f"'{module_path}.{func_name}' is not callable")
        return func

    def _register(self, provider: Provider, plugin_id: str) -> None:
        if provider.name in self._providers:
            raise PluginNameCollision(
                f"provider name collision: '{provider.name}' already registered "
                f"(new attempt from {plugin_id})"
            )
        self._providers[provider.name] = provider


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
