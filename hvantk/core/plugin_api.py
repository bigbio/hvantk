"""Stable contracts for hvantk provider plugins.

In-tree and out-of-tree plugin authors import their types from this module.
The Provider dataclass is constructed by the loader (hvantk.core.plugin_loader)
from a plugin.yaml manifest plus resolved callables - it is not subclassed or
instantiated by plugin authors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Literal, Mapping

Domain = Literal["genomics", "transcriptomics", "proteomics", "epigenomics", "mapping"]
Backend = Literal["hail", "anndata", "pandas"]


class PluginLoadError(Exception):
    """Raised when a plugin fails the protocol or its runtime requirements."""


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
class DatasetSpec:
    """One dataset shipped by a provider plugin.

    `name` is the compound key (e.g. "hgnc:lookup"), composed by the loader
    from plugin.name + ":" + dataset.name. Plugin authors only write the
    bare dataset name in the manifest.
    """

    name: str
    domain: Domain
    backend: Backend
    builder: Callable[..., Any]
    drift_probe: Callable[[], Mapping[str, Any]]
    skill_path: str
    test_paths: TestPaths


@dataclass(frozen=True)
class Provider:
    """Registry record for one provider plugin.

    Constructed by the loader; not subclassed or instantiated by plugin
    authors. Equality is structural (frozen dataclass), so two Provider
    objects with the same fields compare equal - useful in tests.
    """

    name: str
    version: str
    datasets: tuple[DatasetSpec, ...]
