"""Build metadata structs for annotating Hail Tables and MatrixTables globals."""

from datetime import datetime

import hail as hl

SOURCE_DESCRIPTIONS = {
    "ClinVar": "Clinically relevant variant-disease associations from NCBI",
    "ClinGen": "Gene-disease validity curations from ClinGen",
    "dbNSFP": "Functional prediction and conservation scores for coding variants",
    "Ensembl": "Gene annotations from Ensembl BioMart",
    "GeVIR": "Gene-level intolerance metrics (VIRLoF, LOEUF)",
    "gnomAD": "Population allele frequency and constraint metrics from gnomAD",
    "CCR": "Constrained coding regions from Havrilla et al.",
    "Interactome": "Protein-protein interaction interface intervals from INSIDER",
    "HGNC": "HUGO Gene Nomenclature Committee gene identifiers and mappings",
    "UCSC": "Single-cell expression data from UCSC Cell Browser",
    "ExpressionAtlas": "Bulk/single-cell RNA-seq from EMBL-EBI Expression Atlas",
    "CPTAC": "Proteomics expression data from the Clinical Proteomic Tumor Analysis Consortium",
}


def _normalize_source_name(source_name: str) -> str:
    """Normalize source names for case/format-insensitive description lookup.

    The normalization lowercases the input and removes non-alphanumeric
    characters so aliases that differ by case, spaces, hyphens, or punctuation
    resolve to the same metadata key.

    Parameters
    ----------
    source_name : str
        Raw source name to normalize.

    Returns
    -------
    str
        Lowercase alphanumeric normalization of ``source_name``.
    """
    return "".join(ch for ch in source_name.lower() if ch.isalnum())


_NORMALIZED_SOURCE_DESCRIPTIONS = {
    _normalize_source_name(k): v for k, v in SOURCE_DESCRIPTIONS.items()
}


def _get_hvantk_version() -> str:
    """Get hvantk version from package metadata."""
    try:
        from importlib.metadata import version

        return version("hvantk")
    except Exception:
        return "unknown"


def build_table_metadata(
    source_name: str,
    input_path: str,
    ht: hl.Table,
) -> hl.struct:
    """Build an hvantk_metadata struct for a Hail Table.

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source.
    input_path : str
        Path to the raw input file.
    ht : hl.Table
        The Hail Table (used to extract schema info).

    Returns
    -------
    hl.struct
        Metadata struct suitable for ``ht.annotate_globals()``.
    """
    return hl.struct(
        hvantk_version=_get_hvantk_version(),
        source_name=source_name,
        source_description=_NORMALIZED_SOURCE_DESCRIPTIONS.get(
            _normalize_source_name(source_name), ""
        ),
        raw_input_path=input_path,
        build_date=datetime.now().isoformat(),
        reference_genome=str(ht.locus.dtype.reference_genome)
        if "locus" in ht.row
        else "NA",
        row_schema=str(ht.row.dtype),
        key_schema=str(ht.key.dtype),
        key_fields=list(ht.key),
        n_fields=len(ht.row),
    )


def build_matrix_metadata(
    source_name: str,
    input_path: str,
    mt: hl.MatrixTable,
    include_n_cols: bool = False,
) -> hl.struct:
    """Build an hvantk_metadata struct for a Hail MatrixTable.

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source.
    input_path : str
        Path to the raw input file.
    mt : hl.MatrixTable
        The Hail MatrixTable (used to extract schema info).
    include_n_cols : bool, optional
        Whether to materialize and include the number of columns via ``count_cols``.
        Defaults to False to avoid triggering an expensive action; when False,
        ``n_cols`` is set to a missing int64 value.

    Returns
    -------
    hl.struct
        Metadata struct suitable for ``mt.annotate_globals()``.
    """
    return hl.struct(
        hvantk_version=_get_hvantk_version(),
        source_name=source_name,
        source_description=_NORMALIZED_SOURCE_DESCRIPTIONS.get(
            _normalize_source_name(source_name), ""
        ),
        raw_input_path=input_path,
        build_date=datetime.now().isoformat(),
        row_schema=str(mt.row.dtype),
        col_schema=str(mt.col.dtype),
        entry_schema=str(mt.entry.dtype),
        key_schema=str(mt.row_key.dtype),
        key_fields=list(mt.row_key),
        n_row_fields=len(mt.row),
        n_col_fields=len(mt.col),
        n_cols=mt.count_cols() if include_n_cols else hl.missing(hl.tint64),
    )
