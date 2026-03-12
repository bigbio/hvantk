"""
Gene set data structures and I/O utilities.

This module provides:
- GeneSet: Named set of genes with optional metadata
- GeneSetCollection: Collection of gene sets with background universe
- Loading functions for various formats (JSON, GMT, TSV, plain text)
- Simple gene set loading utilities

Supported formats:
- JSON: Full metadata support (canonical format)
- GMT: Tab-separated interchange format for GSEA/MSigDB compatibility
- TSV: Marker gene files with cluster/gene columns
- Plain text: One gene per line
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, TYPE_CHECKING, Union

if TYPE_CHECKING:
    import hail as hl
    import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    # Core data structures
    "GeneSet",
    "GeneSetCollection",
    # Loading functions
    "load_gene_sets_from_dict",
    "load_marker_genes",
    "load_gene_sets",
    # Marker extraction
    "extract_marker_gene_sets",
    # Simple utilities
    "load_gene_set",
    "load_sample_chd_gene_set",
]


# =============================================================================
# Core Data Structures
# =============================================================================


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
    iteration, access, and serialization to JSON and GMT formats.

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

    # -------------------------------------------------------------------------
    # JSON I/O
    # -------------------------------------------------------------------------

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
            f"Loaded {len(collection)} gene sets with "
            f"{len(collection.background_genes)} background genes"
        )

        return collection

    # -------------------------------------------------------------------------
    # GMT I/O
    # -------------------------------------------------------------------------

    def save_gmt(
        self,
        path: Union[str, Path],
        description_field: Optional[str] = None,
    ) -> None:
        """Save to GMT (Gene Matrix Transposed) file.

        GMT format: <set_name>\\t<description>\\t<gene1>\\t<gene2>...

        Note: GMT format does not preserve full metadata. Use JSON for
        complete serialization.

        Parameters
        ----------
        path : Union[str, Path]
            Output file path
        description_field : str, optional
            Metadata field to use as description column. If None, uses
            the source field. If the field is missing, uses empty string.
        """
        path = Path(path)
        logger.info(f"Saving gene set collection to GMT: {path}")

        with open(path, "w") as f:
            for name in sorted(self.gene_sets.keys()):
                gs = self.gene_sets[name]

                # Determine description
                if description_field and description_field in gs.metadata:
                    description = str(gs.metadata[description_field])
                else:
                    description = gs.source or ""

                # Write GMT line: name, description, genes...
                genes_str = "\t".join(sorted(gs.genes))
                f.write(f"{gs.name}\t{description}\t{genes_str}\n")

        logger.info(f"Saved {len(self)} gene sets to GMT: {path}")

    @classmethod
    def load_gmt(
        cls,
        path: Union[str, Path],
        background_genes: Optional[Set[str]] = None,
        background_strategy: str = "union",
    ) -> "GeneSetCollection":
        """Load from GMT (Gene Matrix Transposed) file.

        GMT format: <set_name>\\t<description>\\t<gene1>\\t<gene2>...

        Parameters
        ----------
        path : Union[str, Path]
            Input GMT file path
        background_genes : Set[str], optional
            Explicit background gene set. If None, uses background_strategy.
        background_strategy : str
            Strategy for computing background when background_genes is None:
            - "union": Union of all genes across all sets (default)
            - "none": Empty background (caller must provide later)

        Returns
        -------
        GeneSetCollection
            Loaded gene set collection

        Raises
        ------
        ValueError
            If GMT file has invalid format or background_strategy is unknown
        """
        path = Path(path)
        logger.info(f"Loading gene set collection from GMT: {path}")

        gene_sets = {}
        all_genes: Set[str] = set()

        with open(path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 3:
                    logger.warning(
                        f"GMT line {line_num}: expected at least 3 columns, "
                        f"got {len(parts)}. Skipping."
                    )
                    continue

                name = parts[0]
                description = parts[1]
                genes = set(parts[2:])

                # Filter empty gene entries
                genes = {g for g in genes if g}

                if not genes:
                    logger.warning(
                        f"GMT line {line_num}: gene set '{name}' has no genes. "
                        "Skipping."
                    )
                    continue

                all_genes.update(genes)

                gene_sets[name] = GeneSet(
                    name=name,
                    genes=genes,
                    source=description,
                    metadata={"gmt_description": description},
                )

        # Determine background
        if background_genes is not None:
            bg = background_genes
        elif background_strategy == "union":
            bg = all_genes
        elif background_strategy == "none":
            bg = set()
        else:
            raise ValueError(
                f"Unknown background_strategy: {background_strategy}. "
                "Use 'union' or 'none'."
            )

        collection = cls(
            gene_sets=gene_sets,
            background_genes=bg,
            source_description=f"GMT file: {path.name}",
        )

        logger.info(
            f"Loaded {len(collection)} gene sets with "
            f"{len(collection.background_genes)} background genes from GMT"
        )

        return collection


# =============================================================================
# Loading Functions
# =============================================================================


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
    all_genes: Set[str] = set()

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
    all_genes: Set[str] = set()

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


def load_gene_sets(
    path: Union[str, Path],
    background_genes: Optional[Set[str]] = None,
    **kwargs,
) -> GeneSetCollection:
    """Load gene sets from file, auto-detecting format.

    Supported formats:
    - .json: JSON format with full metadata
    - .gmt: GMT (Gene Matrix Transposed) format

    Parameters
    ----------
    path : Union[str, Path]
        Path to gene set file
    background_genes : Set[str], optional
        Background gene universe. For GMT files, defaults to union of all genes.
    **kwargs
        Additional arguments passed to format-specific loaders

    Returns
    -------
    GeneSetCollection
        Loaded gene set collection

    Raises
    ------
    ValueError
        If file format cannot be determined
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".json":
        collection = GeneSetCollection.load(path)
        if background_genes is not None:
            collection.background_genes = background_genes
        return collection

    elif suffix == ".gmt":
        return GeneSetCollection.load_gmt(
            path,
            background_genes=background_genes,
            **kwargs,
        )

    else:
        raise ValueError(
            f"Unknown gene set file format: {suffix}. " "Supported formats: .json, .gmt"
        )


# =============================================================================
# Marker Gene Extraction
# =============================================================================


def extract_marker_gene_sets(
    summary: Union["hl.Table", "pd.DataFrame", str],
    n_markers: int = 200,
    min_fold_change: float = 1.5,
    min_fraction_expressed: float = 0.1,
    method: str = "fold_change",
    gene_id_field: str = "gene_id",
    gene_name_field: Optional[str] = "gene_name",
) -> GeneSetCollection:
    """Extract top marker genes per group from an expression summary Table.

    Operates on the summary Table produced by
    :func:`~hvantk.utils.matrix_utils.summarize_expression`.  The summary
    Table is small (one row per gene) so this function works entirely in
    pandas — no Hail needed after the initial conversion.

    Parameters
    ----------
    summary : hl.Table, pd.DataFrame, or str
        Summary Table (or path to ``.ht``) with a ``stats`` dict mapping
        group labels to ``struct{mean, fraction_expressed, n_cells}``.
        Also accepts an already-converted pandas DataFrame.
    n_markers : int
        Maximum markers per group.
    min_fold_change : float
        Minimum fold change to qualify as a marker.
    min_fraction_expressed : float
        Minimum fraction of cells expressing a gene in the group.
    method : str
        ``"fold_change"`` (mean_in_group / mean_in_other_groups) or
        ``"specificity"`` (mean_in_group / mean_across_all).
    gene_id_field : str
        Column name for gene IDs in the summary.
    gene_name_field : str or None
        Column name for gene names.  None to omit.

    Returns
    -------
    GeneSetCollection
        One GeneSet per group, with per-gene fold-change scores in metadata.
    """
    import pandas as pd

    if method not in ("fold_change", "specificity"):
        raise ValueError(
            f"Unknown method '{method}'. Use 'fold_change' or 'specificity'."
        )

    df = _summary_to_dataframe(summary, gene_id_field, gene_name_field)

    # Discover groups from the stats columns (prefixed by the group label)
    stat_cols = [c for c in df.columns if c.endswith("_mean")]
    groups = [c.rsplit("_mean", 1)[0] for c in stat_cols]

    if not groups:
        raise ValueError("No groups found in summary table.")

    # Build wide arrays for vectorised computation
    means = pd.DataFrame({g: df[f"{g}_mean"] for g in groups}, index=df.index)
    fracs = pd.DataFrame(
        {g: df[f"{g}_fraction_expressed"] for g in groups}, index=df.index
    )

    all_genes = set(df[gene_id_field].tolist())
    gene_sets: Dict[str, GeneSet] = {}

    for group in groups:
        group_mean = means[group]
        group_frac = fracs[group]

        # Compute fold change
        if method == "fold_change":
            other_cols = [g for g in groups if g != group]
            if other_cols:
                other_mean = means[other_cols].mean(axis=1)
            else:
                other_mean = group_mean  # single group edge case
            # Avoid division by zero
            denom = other_mean.replace(0, 1e-10)
            fc = group_mean / denom
        else:  # specificity
            global_mean = means.mean(axis=1)
            denom = global_mean.replace(0, 1e-10)
            fc = group_mean / denom

        # Apply filters
        mask = (group_frac >= min_fraction_expressed) & (fc >= min_fold_change)
        candidates = df.loc[mask].copy()
        candidates["_fc"] = fc[mask]

        # Rank and take top N
        candidates = candidates.sort_values("_fc", ascending=False).head(n_markers)

        if candidates.empty:
            logger.warning("Group '%s': no markers passed filters.", group)
            continue

        # Use gene names if available, else gene IDs
        use_names = gene_name_field and gene_name_field in candidates.columns
        gene_col = gene_name_field if use_names else gene_id_field
        genes = set(candidates[gene_col].tolist())

        # Store per-gene scores in metadata
        scores = {
            row[gene_col]: round(row["_fc"], 4) for _, row in candidates.iterrows()
        }

        gene_sets[group] = GeneSet(
            name=group,
            genes=genes,
            source=f"marker:{group}",
            metadata={"fold_changes": scores, "method": method},
        )

    logger.info(
        "Extracted markers for %d/%d groups (method=%s, top_n=%d)",
        len(gene_sets),
        len(groups),
        method,
        n_markers,
    )

    return GeneSetCollection(
        gene_sets=gene_sets,
        background_genes=all_genes,
        source_description=(
            f"Marker genes extracted via {method} "
            f"(n={n_markers}, min_fc={min_fold_change}, "
            f"min_frac={min_fraction_expressed})"
        ),
    )


def _summary_to_dataframe(
    summary,
    gene_id_field: str = "gene_id",
    gene_name_field: Optional[str] = "gene_name",
) -> "pd.DataFrame":
    """Convert a summary Table/path/DataFrame to a wide pandas DataFrame.

    The ``stats`` dict column is exploded into ``{group}_mean``,
    ``{group}_fraction_expressed``, ``{group}_n_cells`` columns.
    """
    import pandas as pd

    if isinstance(summary, str):
        import hail as hl

        summary = hl.read_table(summary)

    if not isinstance(summary, pd.DataFrame):
        # Assume Hail Table
        summary = summary.to_pandas()

    if "stats" not in summary.columns:
        # Already a wide DataFrame — return as-is
        return summary

    # Explode the stats dict into wide columns
    rows = []
    for _, row in summary.iterrows():
        flat = {gene_id_field: row.get(gene_id_field)}
        if gene_name_field and gene_name_field in row.index:
            flat[gene_name_field] = row.get(gene_name_field)
        stats = row["stats"]
        if isinstance(stats, dict):
            for group, s in stats.items():
                if hasattr(s, "mean"):
                    flat[f"{group}_mean"] = s.mean
                    flat[f"{group}_fraction_expressed"] = s.fraction_expressed
                    flat[f"{group}_n_cells"] = s.n_cells
                elif isinstance(s, dict):
                    flat[f"{group}_mean"] = s.get("mean", 0)
                    flat[f"{group}_fraction_expressed"] = s.get("fraction_expressed", 0)
                    flat[f"{group}_n_cells"] = s.get("n_cells", 0)
        rows.append(flat)

    return pd.DataFrame(rows)


# =============================================================================
# Simple Gene Set Utilities
# =============================================================================


def load_gene_set(
    path: Optional[Union[str, Path]] = None,
    genes: Optional[Iterable[str]] = None,
    comment_prefix: str = "#",
    delimiter: Optional[str] = None,
    column: int = 0,
    strip_version: bool = False,
) -> Set[str]:
    """Load a gene set from a file and/or an iterable of gene symbols.

    This is a simple utility for loading a single set of genes from a
    plain text file. For structured gene set collections, use
    GeneSetCollection.load() or load_gene_sets().

    Parameters
    ----------
    path : Optional[Union[str, Path]]
        Path to a text file containing gene identifiers. One gene per line
        (default) or delimited lines where the gene column is selected via
        `column`.
    genes : Optional[Iterable[str]]
        Additional genes to include (e.g. passed in code). Can be used
        without `path`.
    comment_prefix : str
        Lines beginning with this prefix are ignored.
    delimiter : Optional[str]
        Delimiter for splitting lines. If None, lines are treated as
        single tokens.
    column : int
        Column index (0-based) to extract when `delimiter` is provided.
    strip_version : bool
        If True, strip transcript / gene version suffix after '.'
        (e.g. ENSG0001.5 -> ENSG0001).

    Returns
    -------
    Set[str]
        Set of unique gene identifiers.

    Raises
    ------
    ValueError
        If neither path nor genes are provided.
    """
    gene_set: Set[str] = set()

    if path is None and genes is None:
        raise ValueError("Provide at least one of 'path' or 'genes'.")

    if path is not None:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Gene set file not found: {p}")
        with p.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(comment_prefix):
                    continue
                if delimiter is not None:
                    parts = line.split(delimiter)
                    if column >= len(parts):
                        continue
                    token = parts[column].strip()
                else:
                    token = line
                if strip_version and "." in token:
                    token = token.split(".", 1)[0]
                if token:
                    gene_set.add(token)

    if genes is not None:
        for g in genes:
            token = g.strip()
            if strip_version and "." in token:
                token = token.split(".", 1)[0]
            if token:
                gene_set.add(token)

    return gene_set


def load_sample_chd_gene_set() -> Set[str]:
    """Return the legacy sample CHD-associated gene set (for backward compatibility)."""
    return {
        "GATA4",
        "NKX2-5",
        "TBX5",
        "NOTCH1",
        "CHD7",
        "TBX1",
        "MYH6",
        "ACTC1",
        "MYH7",
        "TNNT2",
    }
