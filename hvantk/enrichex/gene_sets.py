"""
Gene set data structures and I/O for EnrichEx module.

This module provides the core data structures for representing gene sets
and collections of gene sets, along with functions for loading and saving
them in various formats.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Union

logger = logging.getLogger(__name__)


@dataclass
class GeneSet:
    """A named set of genes with optional metadata.

    Attributes
    ----------
    name : str
        Name of the gene set (e.g., "T_cell", "microglia")
    genes : Set[str]
        Set of gene identifiers (typically gene symbols or Ensembl IDs)
    source : str
        Source description (e.g., "cluster:T_cell", "pathway:apoptosis")
    metadata : Dict[str, Any]
        Additional metadata about the gene set
    """

    name: str
    genes: Set[str]
    source: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_genes(self) -> int:
        """Return the number of genes in this set."""
        return len(self.genes)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns
        -------
        Dict[str, Any]
            Dictionary representation with sorted gene list
        """
        return {
            "name": self.name,
            "genes": sorted(self.genes),  # Sorted for reproducibility
            "source": self.source,
            "n_genes": self.n_genes,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GeneSet":
        """Create GeneSet from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with 'name', 'genes', optional 'source' and 'metadata'

        Returns
        -------
        GeneSet
            GeneSet instance
        """
        return cls(
            name=data["name"],
            genes=set(data["genes"]),
            source=data.get("source", ""),
            metadata=data.get("metadata", {}),
        )


@dataclass
class GeneSetCollection:
    """Collection of gene sets with a common background universe.

    This class manages multiple gene sets and provides methods for
    iteration, access, and serialization.

    Attributes
    ----------
    gene_sets : Dict[str, GeneSet]
        Dictionary mapping gene set names to GeneSet objects
    background_genes : Set[str]
        Background universe of all possible genes
    source_description : str
        Description of the source of these gene sets
    metadata : Dict[str, Any]
        Additional metadata about the collection
    """

    gene_sets: Dict[str, GeneSet]
    background_genes: Set[str]
    source_description: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        """Return number of gene sets in the collection."""
        return len(self.gene_sets)

    def __iter__(self):
        """Iterate over GeneSet objects."""
        return iter(self.gene_sets.values())

    def get(self, name: str) -> Optional[GeneSet]:
        """Get a gene set by name.

        Parameters
        ----------
        name : str
            Name of the gene set

        Returns
        -------
        Optional[GeneSet]
            GeneSet if found, None otherwise
        """
        return self.gene_sets.get(name)

    def names(self) -> List[str]:
        """Return list of gene set names.

        Returns
        -------
        List[str]
            Sorted list of gene set names
        """
        return sorted(self.gene_sets.keys())

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns
        -------
        Dict[str, Any]
            Dictionary representation with all gene sets
        """
        return {
            "gene_sets": {name: gs.to_dict() for name, gs in self.gene_sets.items()},
            "background_genes": sorted(self.background_genes),
            "n_gene_sets": len(self),
            "n_background": len(self.background_genes),
            "source_description": self.source_description,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GeneSetCollection":
        """Create GeneSetCollection from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with 'gene_sets', 'background_genes', etc.

        Returns
        -------
        GeneSetCollection
            GeneSetCollection instance
        """
        gene_sets = {}
        for name, gs_data in data.get("gene_sets", {}).items():
            gene_sets[name] = GeneSet.from_dict(gs_data)

        return cls(
            gene_sets=gene_sets,
            background_genes=set(data.get("background_genes", [])),
            source_description=data.get("source_description", ""),
            metadata=data.get("metadata", {}),
        )

    def save(self, path: Union[str, Path]) -> None:
        """Save to JSON file.

        Parameters
        ----------
        path : Union[str, Path]
            Output file path
        """
        path = Path(path)
        logger.info(f"Saving gene set collection to {path}")

        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

        logger.info(f"Saved {len(self)} gene sets to {path}")

    @classmethod
    def load(cls, path: Union[str, Path]) -> "GeneSetCollection":
        """Load from JSON file.

        Parameters
        ----------
        path : Union[str, Path]
            Input file path

        Returns
        -------
        GeneSetCollection
            Loaded gene set collection
        """
        path = Path(path)
        logger.info(f"Loading gene set collection from {path}")

        with open(path) as f:
            data = json.load(f)

        collection = cls.from_dict(data)
        logger.info(
            f"Loaded {len(collection)} gene sets with {len(collection.background_genes)} background genes"
        )

        return collection


def load_gene_sets_from_dict(
    gene_sets_dict: Dict[str, List[str]],
    background_genes: Optional[Set[str]] = None,
    source: str = "user_provided",
) -> GeneSetCollection:
    """Create GeneSetCollection from simple {name: [genes]} dictionary.

    This is a convenience function for creating gene set collections
    from simple Python dictionaries.

    Parameters
    ----------
    gene_sets_dict : Dict[str, List[str]]
        Dictionary mapping gene set names to gene lists
    background_genes : Set[str], optional
        Background gene universe. If None, uses union of all genes.
    source : str
        Source description

    Returns
    -------
    GeneSetCollection
        Gene set collection

    Examples
    --------
    >>> gene_sets = {
    ...     "T_cell": ["CD3D", "CD3E", "CD3G"],
    ...     "B_cell": ["CD19", "CD79A", "CD79B"]
    ... }
    >>> collection = load_gene_sets_from_dict(gene_sets)
    >>> print(len(collection))
    2
    """
    gene_sets = {}
    all_genes = set()

    for name, genes in gene_sets_dict.items():
        genes_set = set(genes)
        all_genes.update(genes_set)
        gene_sets[name] = GeneSet(name=name, genes=genes_set, source=source)

    if background_genes is None:
        background_genes = all_genes

    return GeneSetCollection(
        gene_sets=gene_sets,
        background_genes=background_genes,
        source_description=source,
    )


def load_marker_genes(
    marker_file: Union[str, Path],
    cluster_column: str = "cluster",
    gene_column: str = "gene",
    background_genes: Optional[Set[str]] = None,
) -> GeneSetCollection:
    """Load pre-computed marker genes from TSV file.

    Expected format: TSV with at least cluster and gene columns.

    Parameters
    ----------
    marker_file : Union[str, Path]
        Path to marker gene file (TSV format)
    cluster_column : str
        Column name for cluster identifiers
    gene_column : str
        Column name for gene identifiers
    background_genes : Set[str], optional
        Background gene universe. If None, uses all genes in file.

    Returns
    -------
    GeneSetCollection
        Gene set collection with one set per cluster

    Examples
    --------
    >>> collection = load_marker_genes("markers.tsv", cluster_column="cell_type")
    """
    import pandas as pd

    marker_file = Path(marker_file)
    logger.info(f"Loading marker genes from {marker_file}")

    df = pd.read_csv(marker_file, sep="\t")

    if cluster_column not in df.columns:
        raise ValueError(f"Column '{cluster_column}' not found in {marker_file}")
    if gene_column not in df.columns:
        raise ValueError(f"Column '{gene_column}' not found in {marker_file}")

    gene_sets = {}
    all_genes = set()

    for cluster, group in df.groupby(cluster_column):
        genes = set(group[gene_column].tolist())
        all_genes.update(genes)
        gene_sets[str(cluster)] = GeneSet(
            name=str(cluster),
            genes=genes,
            source=f"markers:{cluster}",
        )

    if background_genes is None:
        background_genes = all_genes

    logger.info(f"Loaded {len(gene_sets)} gene sets from {marker_file}")

    return GeneSetCollection(
        gene_sets=gene_sets,
        background_genes=background_genes,
        source_description=f"Marker genes from {marker_file.name}",
    )
