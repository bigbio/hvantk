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

# Sentinel value for ``probe_status`` marking a drift probe as an intentional
# stub: a documentation-only / license-gated / publication-only source with no
# stable, programmatically-probeable direct URL. ``drift_runner`` reports these
# as status="stub" (a visible WARNING) instead of a silent false-green "clean"
# or a misleading "probe_failed". See issue #177.
PROBE_STATUS_STUB = "stub"

# Honest provenance token recorded by ``run_builder._coerce_fingerprint`` for a
# stubbed source (the ``fingerprint`` key wins there). Self-describing rather
# than a fake ``sha256:...`` hash, so provenance never implies a real probe ran.
STUB_FINGERPRINT_TOKEN = "stub:no-programmatic-source"

# Value a hand-seeded ``drift_fingerprint.json`` carries where a real checksum
# belongs. A baseline holding it was written by hand, never captured from a live
# probe, so it cannot equal any observed fingerprint.
PLACEHOLDER_CHECKSUM = "placeholder-regenerate-from-live-source"

# ``fetched_at`` value meaning "no probe has ever run": a hand-written baseline
# stamped at the Unix epoch rather than at a real fetch time.
_EPOCH_PREFIX = "1970-01-01"


def placeholder_baseline_reason(expected: Mapping[str, Any]) -> str | None:
    """Explain why ``expected`` is a hand-seeded baseline, or None if it is real.

    A committed ``drift_fingerprint.json`` is supposed to be the output of a real
    probe run. Six of them were instead seeded by hand -- placeholder checksum
    strings, ``fetched_at`` at the Unix epoch -- and a hand-seeded baseline can
    never equal a live observation. Comparing it produced a permanent, and
    therefore meaningless, ``drifted`` verdict: the scheduled drift bot opened the
    same no-op pull requests every night, which is how six dead comparators went
    unnoticed while looking maximally alive.

    ``drift_runner`` calls this before diffing and reports ``probe_failed`` -- the
    baseline is missing in substance even though the file exists -- so the
    condition surfaces as the configuration error it is, and is never mistaken
    for the upstream having moved.

    Three markers, none of which can occur in genuine probe output: the placeholder
    sentinel string, an empty checksum value, and ``fetched_at`` at the Unix epoch.

    Detection is deliberately narrow. An empty ``checksums`` *map* is not a marker,
    because probes that fingerprint HTTP validators instead of bodies legitimately
    ship one -- peptideatlas does. The distinction between an empty map and an empty
    value inside it is load-bearing, and both cases are pinned by tests.
    """
    checksums = expected.get("checksums")
    if isinstance(checksums, Mapping):
        for name, value in checksums.items():
            if isinstance(value, str) and value == PLACEHOLDER_CHECKSUM:
                return f"checksum for {name!r} is the placeholder sentinel"
            # The other seeding tell: cptac's two baselines carried an empty string
            # where a digest belongs. Those also had an epoch ``fetched_at``, so the
            # check below covered the real cases -- but a baseline hand-written with
            # a genuine timestamp would slip past. A probe that ran always produces
            # a digest.
            if isinstance(value, str) and not value.strip():
                return f"checksum for {name!r} is empty; no digest was ever computed"

    fetched_at = expected.get("fetched_at")
    if isinstance(fetched_at, str) and fetched_at.startswith(_EPOCH_PREFIX):
        return f"fetched_at is the Unix epoch ({fetched_at!r}); no probe ever ran"

    return None


def stub_fingerprint(reason: str) -> dict:
    """Build the structured sentinel a documentation-only drift probe returns.

    Use this from ``fetch_fingerprint()`` when the source cannot be fingerprinted
    programmatically (manual/gated/publication-only acquisition). ``reason``
    explains why (surfaced in the ``hvantk drift`` WARNING). Returning this marks
    the dataset as status="stub" rather than producing a false-green drift result.
    """
    return {
        "probe_status": PROBE_STATUS_STUB,
        "reason": reason,
        "fingerprint": STUB_FINGERPRINT_TOKEN,
    }


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

    #: Fields naming an on-disk validation artifact. ``command`` is excluded -- it is a
    #: shell string, not a path.
    ARTIFACT_FIELDS = ("fixture", "schema_snapshot", "row_snapshot", "drift_fingerprint")

    def missing_artifacts(self) -> tuple[tuple[str, str], ...]:
        """Return ``(field, path)`` for every declared artifact absent from disk.

        The loader resolves these paths but deliberately does not require them, so a
        manifest can declare a snapshot it does not ship and still load -- which is how
        the tree came to hold 25 dataset declarations against 10 snapshot files without
        anything failing. Callers that want the contract enforced (``hvantk plugins
        validate``, the coverage ratchet in the test suite) ask for the gap explicitly
        rather than each re-deriving the path layout.

        Returns an empty tuple when every declared artifact exists.
        """
        missing = []
        for field in self.ARTIFACT_FIELDS:
            path = getattr(self, field)
            if not Path(path).exists():
                missing.append((field, path))
        return tuple(missing)


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
