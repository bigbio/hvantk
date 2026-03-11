"""
Gene alias expansion utility using HGNC data.

Expands a gene set to include all known HGNC aliases and previous symbols,
ensuring that ClinVar gene filtering does not silently miss variants due
to symbol mismatches.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


def _load_hgnc_symbol_maps(
    hgnc_path: str,
) -> Tuple[Set[str], Dict[str, str], Dict[str, List[str]]]:
    """Load HGNC data and build symbol lookup maps.

    Supports both Hail Table (.ht) and TSV formats.

    Returns
    -------
    Tuple of:
    - canonical_symbols: Set of all current approved symbols
    - alias_to_canonical: Mapping from alias/prev_symbol → canonical symbol
    - canonical_to_aliases: Mapping from canonical symbol → list of aliases
    """
    path = Path(hgnc_path)

    if path.suffix == ".ht" or (path.is_dir() and (path / "_SUCCESS").exists()):
        return _load_from_hail_table(hgnc_path)
    else:
        return _load_from_tsv(hgnc_path)


def _load_from_hail_table(
    hgnc_path: str,
) -> Tuple[Set[str], Dict[str, str], Dict[str, List[str]]]:
    """Load symbol maps from an HGNC Hail Table."""
    import hail as hl

    ht = hl.read_table(hgnc_path)
    row_fields = set(ht.row)

    # Collect gene_symbol and alias/prev fields
    fields_to_select = ["gene_symbol"]
    if "alias_symbols" in row_fields:
        fields_to_select.append("alias_symbols")
    if "prev_symbols" in row_fields:
        fields_to_select.append("prev_symbols")

    data = ht.select(*fields_to_select).collect()

    canonical_symbols: Set[str] = set()
    alias_to_canonical: Dict[str, str] = {}
    canonical_to_aliases: Dict[str, List[str]] = {}

    # Pass 1: collect all canonical symbols
    for row in data:
        canonical_symbols.add(row.gene_symbol)

    # Pass 2: build alias mappings (canonical set is complete)
    for row in data:
        symbol = row.gene_symbol
        all_aliases: List[str] = []

        if "alias_symbols" in fields_to_select and row.alias_symbols:
            for alias in row.alias_symbols:
                if alias and alias not in canonical_symbols:
                    if alias not in alias_to_canonical:
                        alias_to_canonical[alias] = symbol
                        all_aliases.append(alias)

        if "prev_symbols" in fields_to_select and row.prev_symbols:
            for prev in row.prev_symbols:
                if prev and prev not in canonical_symbols:
                    if prev not in alias_to_canonical:
                        alias_to_canonical[prev] = symbol
                        all_aliases.append(prev)

        if all_aliases:
            canonical_to_aliases[symbol] = all_aliases

    return canonical_symbols, alias_to_canonical, canonical_to_aliases


def _load_from_tsv(
    hgnc_path: str,
) -> Tuple[Set[str], Dict[str, str], Dict[str, List[str]]]:
    """Load symbol maps from an HGNC TSV file."""
    import pandas as pd

    df = pd.read_csv(hgnc_path, sep="\t", dtype=str, na_values=[""])

    # Normalize column names (HGNC uses 'symbol' or 'Approved symbol')
    col_map = {}
    for col in df.columns:
        lower = col.lower().replace(" ", "_")
        if lower in ("symbol", "approved_symbol"):
            col_map[col] = "gene_symbol"
        elif lower in ("alias_symbol", "alias_symbols", "alias_name"):
            col_map[col] = "alias_symbols"
        elif lower in ("prev_symbol", "prev_symbols", "previous_symbols"):
            col_map[col] = "prev_symbols"
    df = df.rename(columns=col_map)

    if "gene_symbol" not in df.columns:
        raise ValueError(
            f"HGNC TSV missing required column 'symbol' or 'Approved symbol'. "
            f"Found columns: {list(df.columns)}"
        )

    canonical_symbols: Set[str] = set()
    alias_to_canonical: Dict[str, str] = {}
    canonical_to_aliases: Dict[str, List[str]] = {}

    # Pass 1: collect all canonical symbols
    for _, row in df.iterrows():
        symbol = row.get("gene_symbol")
        if not symbol or pd.isna(symbol):
            continue
        canonical_symbols.add(str(symbol).strip())

    # Pass 2: build alias mappings (canonical set is complete)
    for _, row in df.iterrows():
        symbol = row.get("gene_symbol")
        if not symbol or pd.isna(symbol):
            continue
        symbol = str(symbol).strip()

        all_aliases: List[str] = []

        # Parse pipe-separated alias_symbols
        alias_str = row.get("alias_symbols")
        if alias_str and not pd.isna(alias_str):
            for alias in str(alias_str).split("|"):
                alias = alias.strip()
                if alias and alias not in canonical_symbols:
                    if alias not in alias_to_canonical:
                        alias_to_canonical[alias] = symbol
                        all_aliases.append(alias)

        # Parse pipe-separated prev_symbols
        prev_str = row.get("prev_symbols")
        if prev_str and not pd.isna(prev_str):
            for prev in str(prev_str).split("|"):
                prev = prev.strip()
                if prev and prev not in canonical_symbols:
                    if prev not in alias_to_canonical:
                        alias_to_canonical[prev] = symbol
                        all_aliases.append(prev)

        if all_aliases:
            canonical_to_aliases[symbol] = all_aliases

    return canonical_symbols, alias_to_canonical, canonical_to_aliases


def expand_gene_set_with_aliases(
    genes: List[str],
    hgnc_path: Optional[str] = None,
) -> Tuple[Set[str], Dict[str, str]]:
    """Expand a gene set to include all known HGNC aliases.

    For each user-provided gene symbol:
    - If it is a canonical symbol, also include its known aliases/prev_symbols
    - If it is an alias/prev_symbol, also include the canonical symbol

    This ensures that ClinVar's GENEINFO field (which may use any of these
    forms) is matched regardless of which symbol the user provides.

    Parameters
    ----------
    genes : list of str
        User-provided gene symbols.
    hgnc_path : str, optional
        Path to HGNC TSV or Hail Table. If None, returns the original set
        unchanged (no expansion).

    Returns
    -------
    Tuple of:
    - Expanded gene set (original + aliases + canonical symbols)
    - Mapping of resolved aliases: alias → canonical symbol (for logging).
      Only includes entries where a symbol was resolved to a different symbol.
    """
    if not hgnc_path:
        return set(genes), {}

    (
        canonical_symbols,
        alias_to_canonical,
        canonical_to_aliases,
    ) = _load_hgnc_symbol_maps(hgnc_path)

    expanded: Set[str] = set(genes)
    alias_map: Dict[str, str] = {}

    for gene in genes:
        if gene in canonical_symbols:
            # User provided a canonical symbol → add its aliases
            aliases = canonical_to_aliases.get(gene, [])
            for alias in aliases:
                if alias not in expanded:
                    expanded.add(alias)
                    alias_map[alias] = gene
        elif gene in alias_to_canonical:
            # User provided an alias → add the canonical symbol
            canonical = alias_to_canonical[gene]
            if canonical not in expanded:
                expanded.add(canonical)
                alias_map[gene] = canonical
            # Also add other aliases of the same canonical symbol
            for alias in canonical_to_aliases.get(canonical, []):
                if alias not in expanded:
                    expanded.add(alias)
                    alias_map[alias] = canonical

    return expanded, alias_map
