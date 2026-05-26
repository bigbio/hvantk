"""Stable contracts for hvantk provider plugins.

In-tree and out-of-tree plugin authors import their types from this module.
The Provider dataclass is constructed by the loader (hvantk.core.plugin.loader)
from a plugin.yaml manifest plus resolved callables - it is not subclassed or
instantiated by plugin authors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Literal, Mapping, Protocol

Domain = Literal["genomics", "transcriptomics", "proteomics", "epigenomics", "mapping"]
Backend = Literal["hail", "anndata", "pandas"]

# Keys excluded from both the canonical fingerprint hash (run_builder) and the
# drift-comparison diff (drift_runner). `fetched_at` is timestamp noise.
# `probe_version` is orthogonal probe-implementation metadata; bumping it should
# not flip drift status or invalidate stored artifact fingerprints.
PROBE_FINGERPRINT_IGNORED_KEYS = frozenset({"fetched_at", "probe_version"})


class Builder(Protocol):
    """Phase B builder contract used by ``DatasetSpec.builder``.

    Plugin authors implement this signature; the platform's
    ``run_builder_for_spec`` invokes it after computing the source fingerprint
    from the drift probe and constructing a ``BuildContext``.

    ``parsed_input`` is whatever ``DatasetSpec.parse_fn`` returned (typically a
    path or a dict of paths for multi-input builders). Returns an Artifact
    subclass — ``AnnotationTable``, ``ExpressionMatrix``, or ``GeneSet`` —
    matching the manifest's ``artifact_type`` declaration.
    """

    def __call__(
        self,
        parsed_input: Any,
        ctx: "BuildContext",
        **params: Any,
    ) -> Any: ...


class DownloadFn(Protocol):
    """Optional download stage, declared via the ``lifecycle.download`` block
    in ``plugin.yaml`` (api_version >= 2). Fetches upstream data and writes
    raw files under ``raw_dir``. Extra kwargs come from ``--plugin-arg
    KEY=VALUE`` passthrough; the function ignores ones it doesn't consume.
    """

    def __call__(
        self,
        *,
        raw_dir: Path | str,
        **params: Any,
    ) -> None: ...


class ParseFn(Protocol):
    """Optional parse stage, declared via the ``lifecycle.parse`` block in
    ``plugin.yaml``. Reads the raw files ``download_fn`` produced from
    ``raw_dir`` and writes the intermediate representation that ``Builder``
    consumes to ``output_path``.
    """

    def __call__(
        self,
        *,
        raw_dir: Path | str,
        output_path: Path | str,
        **params: Any,
    ) -> None: ...


class PluginLoadError(Exception):
    """Raised when a plugin fails the protocol or its runtime requirements."""


class PluginNameCollision(PluginLoadError):
    """Raised when a plugin's name conflicts with an already-registered plugin.

    Subclass of PluginLoadError so callers that catch PluginLoadError still
    see collisions, but loaders that want to distinguish 'hard-fail collision'
    from 'soft-fail per-plugin failure' can catch this specifically.
    """


class DriftProbeError(Exception):
    """Raised by a drift probe on transient failure (network, timeout, parse).

    Distinct from drift itself: the runner classifies a DriftProbeError as
    status=probe_failed, separate from a status=drifted diff.
    """


@dataclass(frozen=True)
class TestPaths:
    """Validation artifacts for one dataset, as paths resolved from plugin.yaml.

    All paths are absolute (resolved by the loader) so callers do not need to
    know where the plugin folder lives.
    """

    __test__ = False  # not a pytest test class; pytest sees the leading "Test" and tries to collect

    command: str
    fixture: str
    schema_snapshot: str
    row_snapshot: str
    drift_fingerprint: str


@dataclass(frozen=True)
class DatasetManifest:
    """Descriptive view of a dataset declaration — no resolved callables.

    Populated by the loader's first pass (pure YAML reads). Catalog,
    ``hvantk plugins list``, and any UI that needs to enumerate datasets
    use this; it survives missing optional runtimes (e.g. hail) that
    would otherwise prevent the corresponding DatasetSpec from binding.

    To get an executable spec with resolved callables, call
    ``registry.get_dataset(name)`` — the loader resolves callables on
    demand and caches the resulting ``DatasetSpec``.
    """

    __test__ = False  # not a pytest test class

    name: str  # compound, e.g. "clinvar:variants"
    domain: Domain
    backend: Backend
    skill_path: str
    test_paths: TestPaths
    plugin_name: str
    plugin_version: str | None = None
    artifact_type_name: str | None = None  # e.g. "AnnotationTable"; resolved lazily
    schema_id: str | None = None
    has_download_fn: bool = False
    has_parse_fn: bool = False
    # Module/function references kept as strings so we can resolve lazily
    builder_ref: tuple[str, str] = field(default=("", ""))
    drift_probe_ref: tuple[str, str] = field(default=("", ""))
    download_ref: tuple[str, str] | None = None
    parse_ref: tuple[str, str] | None = None


@dataclass(frozen=True)
class DatasetSpec:
    """One dataset shipped by a provider plugin.

    `name` is the compound key (e.g. "hgnc:lookup"), composed by the loader
    from plugin.name + ":" + dataset.name. Plugin authors only write the
    bare dataset name in the manifest.

    Optional lifecycle callables (`download_fn`, `parse_fn`) are wired in by
    the loader when a manifest declares an api_version >= 2 `lifecycle:` block.
    Plugin contract for keyword arguments:

      download_fn(raw_dir=<path>)
          Fetch upstream data and write raw files under `raw_dir`.
      parse_fn(raw_dir=<path>, output_path=<path>)
          Read the raw files produced by `download_fn` from `raw_dir` and
          write the intermediate representation that the builder consumes
          to `output_path`.

    Both default to None for backward compatibility with api_version 1
    manifests; callers must check for None before invoking.
    """

    name: str
    domain: Domain
    backend: Backend
    builder: Builder
    drift_probe: Callable[[], Mapping[str, Any]]
    skill_path: str
    test_paths: TestPaths
    download_fn: DownloadFn | None = None
    parse_fn: ParseFn | None = None
    artifact_type: type | None = None
    schema_id: str | None = None
    plugin_version: str | None = None


@dataclass(frozen=True)
class Provider:
    """Registry record for one provider plugin.

    Constructed by the loader; not subclassed or instantiated by plugin
    authors. Equality is structural (frozen dataclass), so two Provider
    objects with the same fields compare equal - useful in tests.

    `catalog_path` is the absolute path to the plugin's per-plugin
    `catalog/datasets.json` (or None if the plugin does not ship one).
    Consumers load it lazily via JSON to avoid forcing a parse cost on
    plugin discovery.

    `primary_domain` records the manifest-declared dataset domain for the
    provider (the most common `datasets[].domain` value in plugin.yaml).
    It is populated independently of whether each builder import
    succeeded, so catalog-consumers can route entries even when a
    plugin's runtime dependencies (e.g. hail) are unavailable.

    `manifests` holds the full descriptive view for all declared datasets,
    populated by the loader's first (pure-YAML) pass. This is always
    populated regardless of whether callable imports succeeded.
    `datasets` holds only the successfully-bound executable specs.
    """

    name: str
    version: str
    datasets: tuple[DatasetSpec, ...]
    catalog_path: str | None = None
    primary_domain: str | None = None
    manifests: tuple[DatasetManifest, ...] = field(default=())


from hvantk.core.models.build_context import BuildContext  # noqa: F401
