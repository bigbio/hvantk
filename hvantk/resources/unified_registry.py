"""
Unified registry interface for hvantk.

Sources data from two places:

1. Per-plugin catalogs declared in each plugin's ``plugin.yaml`` via the
   ``catalog:`` field (resolved by ``hvantk.core.plugin_loader``).
2. Legacy per-domain JSON files under ``hvantk/resources/registry/<omics>/``.
   After the per-plugin migration, only ``genomics/datasets.json`` remains
   in-tree (it still holds orphan entries without an owning plugin).

**Collision precedence:** per-plugin catalogs win. If a legacy registry
entry duplicates an ``(omics_type, accession)`` already contributed by a
plugin catalog, the legacy entry is skipped and a WARNING is logged. There
are zero overlapping accessions today; the rule exists to make the
behavior deterministic the day someone moves an orphan entry into a plugin
catalog without remembering to delete the legacy copy.

The public API (``list_*_datasets``, ``get_dataset``, ``search``,
``get_stats``) is unchanged.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# Map a plugin's primary `domain` (declared on its DatasetSpec entries) to the
# omics bucket. Domains and omics types currently use the same vocabulary, but
# we go through this dict so the mapping is explicit and easy to override.
_DOMAIN_TO_OMICS = {
    "transcriptomics": "transcriptomics",
    "proteomics": "proteomics",
    "genomics": "genomics",
    "epigenomics": "epigenomics",
    # mapping-type plugins (e.g. hgnc lookup) carry no per-omics data
    # entries and so contribute nothing to the registry buckets.
}


def _provider_primary_domain(provider) -> Optional[str]:
    """Return the dominant `domain` for this provider.

    Prefer the loader's manifest-derived `primary_domain` (populated even
    when a builder import fails). Fall back to inspecting the loaded
    DatasetSpec list for compatibility with Provider instances built by
    older loader versions.
    """
    primary = getattr(provider, "primary_domain", None)
    if primary:
        return primary
    if not provider.datasets:
        return None
    counts: Dict[str, int] = {}
    for ds in provider.datasets:
        counts[ds.domain] = counts.get(ds.domain, 0) + 1
    return max(counts, key=counts.get)


# Per-entry fields that unambiguously identify a proteomics record. Used to
# route entries that share a catalog with transcriptomics records (e.g. the
# Human Protein Atlas E-PROT-* entries inside the expression-atlas catalog).
_PROTEOMICS_FIELDS = {
    "abundance_unit",
    "protein_database",
    "protein_inference_method",
    "quantification_method",
}


def _infer_omics_for_entry(entry: Dict[str, Any], fallback: Optional[str]) -> Optional[str]:
    """Pick an omics bucket for a single catalog entry.

    Most entries adopt the owning plugin's primary domain (``fallback``).
    A small number of provider catalogs mix domains - the expression-atlas
    catalog, for example, contains both bulk-RNA-seq entries and Human
    Protein Atlas E-PROT-* proteomics entries. For those we look at the
    entry payload directly.
    """
    accession = (entry.get("accession") or "").upper()
    if accession.startswith("E-PROT") or any(k in entry for k in _PROTEOMICS_FIELDS):
        return "proteomics"
    data_level = (entry.get("data_level") or "").lower()
    if data_level == "protein":
        return "proteomics"
    return fallback


class HvantkRegistry:
    """Main interface for accessing all omics datasets in hvantk."""

    omics_types = ["transcriptomics", "proteomics", "genomics", "epigenomics"]

    def __init__(self, registry_root: Optional[Path] = None):
        """Initialize the unified registry."""
        self.registry_root = registry_root or Path(__file__).parent / "registry"
        self._cache: Dict[str, List[Dict[str, Any]]] = {
            t: [] for t in self.omics_types
        }
        self._load_registry()

    def _load_registry(self) -> None:
        """Aggregate per-plugin catalogs and the remaining legacy registry files.

        Plugin catalogs are loaded first; the legacy genomics registry is the
        fallback. Any legacy entry whose ``(omics_type, accession)`` is already
        covered by a plugin catalog is skipped with a WARNING — per-plugin wins.
        """
        # 1. Per-plugin catalogs (declared via plugin.yaml `catalog:`)
        try:
            # Import locally to avoid a hard dependency cycle when this module
            # is imported before the plugin loader is needed.
            from hvantk.core.plugin import loader as plugin_loader

            registry = plugin_loader.get_registry()
        except Exception as exc:  # noqa: BLE001
            logger.warning("plugin loader unavailable; per-plugin catalogs skipped: %s", exc)
            registry = None

        if registry is not None:
            for provider in registry.list_providers():
                catalog_path = getattr(provider, "catalog_path", None)
                if not catalog_path:
                    continue
                try:
                    with open(catalog_path, "r") as fh:
                        entries = json.load(fh)
                except (OSError, json.JSONDecodeError) as exc:
                    logger.warning(
                        "failed to load catalog for %s (%s): %s",
                        provider.name,
                        catalog_path,
                        exc,
                    )
                    continue

                primary = _DOMAIN_TO_OMICS.get(
                    _provider_primary_domain(provider) or "", None
                )
                if primary is None:
                    # Mapping-type or unknown-domain plugins do not feed
                    # the omics buckets.
                    continue
                for entry in entries:
                    omics = _infer_omics_for_entry(entry, primary)
                    if omics is None:
                        continue
                    self._cache.setdefault(omics, []).append(entry)

        # Snapshot the accessions contributed by the plugin layer; the legacy
        # pass below uses this set to enforce the "per-plugin wins" precedence.
        plugin_owned: Dict[str, set] = {
            otype: {
                e.get("accession") for e in self._cache.get(otype, []) if e.get("accession")
            }
            for otype in self.omics_types
        }

        # 2. Legacy per-omics registry files (only genomics remains today)
        for omics_type in self.omics_types:
            datasets_file = self.registry_root / omics_type / "datasets.json"
            if not datasets_file.is_file():
                continue
            try:
                with open(datasets_file, "r") as fh:
                    legacy_entries = json.load(fh)
            except (OSError, json.JSONDecodeError) as exc:
                logger.error("Error loading %s: %s", omics_type, exc)
                continue
            for entry in legacy_entries:
                accession = entry.get("accession")
                if accession and accession in plugin_owned.get(omics_type, set()):
                    logger.warning(
                        "legacy %s entry %r duplicates a per-plugin catalog "
                        "accession; skipping legacy copy (per-plugin wins per "
                        "unified_registry precedence rule)",
                        omics_type,
                        accession,
                    )
                    continue
                self._cache[omics_type].append(entry)

    def list_transcriptomics_datasets(self) -> List[Dict]:
        """Get all transcriptomics datasets."""
        return self._cache.get("transcriptomics", [])

    def list_proteomics_datasets(self) -> List[Dict]:
        """Get all proteomics datasets."""
        return self._cache.get("proteomics", [])

    def list_genomics_datasets(self) -> List[Dict]:
        """Get all genomics datasets."""
        return self._cache.get("genomics", [])

    def list_epigenomics_datasets(self) -> List[Dict]:
        """Get all epigenomics datasets."""
        return self._cache.get("epigenomics", [])

    def get_dataset(self, accession: str) -> Optional[Dict]:
        """Get a specific dataset by accession across all omics types."""
        for omics_type in self.omics_types:
            for dataset in self._cache.get(omics_type, []):
                if dataset.get("accession") == accession:
                    result = dataset.copy()
                    result["_omics_type"] = omics_type
                    return result
        return None

    def search(
        self,
        query: str = "",
        omics_type: Optional[str] = None,
        organism: Optional[str] = None,
        data_source: Optional[str] = None,
    ) -> List[Dict]:
        """Search datasets with various filters.

        All filters are case-insensitive substring matches. ``query`` matches
        against title + description + accession; ``organism`` and ``data_source``
        match against their named fields. Substring semantics keep the CLI's
        ``--organism Homo`` ergonomic without forcing callers to know the exact
        canonical form (e.g. ``Homo sapiens``).
        """
        results = []
        search_types = [omics_type] if omics_type else self.omics_types
        organism_lc = organism.lower() if organism else None
        data_source_lc = data_source.lower() if data_source else None
        query_lc = query.lower() if query else None

        for otype in search_types:
            for dataset in self._cache.get(otype, []):
                if query_lc:
                    searchable = (
                        f"{dataset.get('title', '')} "
                        f"{dataset.get('description', '')} "
                        f"{dataset.get('accession', '')}"
                    ).lower()
                    if query_lc not in searchable:
                        continue

                if organism_lc and organism_lc not in (
                    dataset.get("organism") or ""
                ).lower():
                    continue

                if data_source_lc and data_source_lc not in (
                    dataset.get("data_source") or ""
                ).lower():
                    continue

                result = dataset.copy()
                result["_omics_type"] = otype
                results.append(result)

        return results

    def get_stats(self) -> Dict:
        """Get registry statistics."""
        stats = {
            "total_datasets": sum(len(datasets) for datasets in self._cache.values()),
            "by_omics_type": {
                otype: len(datasets) for otype, datasets in self._cache.items()
            },
            "organisms": {},
            "data_sources": {},
        }

        # Count organisms and data sources
        for omics_type, datasets in self._cache.items():
            for dataset in datasets:
                organism = dataset.get("organism", "Unknown")
                data_source = dataset.get("data_source", "Unknown")

                stats["organisms"][organism] = stats["organisms"].get(organism, 0) + 1
                stats["data_sources"][data_source] = (
                    stats["data_sources"].get(data_source, 0) + 1
                )

        return stats

    def reload(self):
        """Reload registry from files."""
        self._cache = {t: [] for t in self.omics_types}
        self._load_registry()


# Backward compatibility - maintain existing function signatures
def load_transcriptomics_registry() -> List[Dict]:
    """Backward compatibility function for transcriptomics data."""
    registry = HvantkRegistry()
    return registry.list_transcriptomics_datasets()


def load_expression_atlas_datasets() -> List[Dict]:
    """Backward compatibility function for Expression Atlas data."""
    registry = HvantkRegistry()
    return [
        d
        for d in registry.list_transcriptomics_datasets()
        if d.get("data_source") == "Expression_Atlas"
    ]


def load_ucsc_datasets() -> List[Dict]:
    """Backward compatibility function for UCSC data."""
    registry = HvantkRegistry()
    return [
        d
        for d in registry.list_transcriptomics_datasets()
        if d.get("data_source") == "UCSC"
    ]
