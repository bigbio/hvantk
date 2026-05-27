"""
Gene ID mapping utility using HGNC data.

This module provides the GeneMapper class for bidirectional gene ID mapping
using HGNC (HUGO Gene Nomenclature Committee) data as the central hub.
"""

import logging
from typing import Dict, List, Literal, Optional

import hail as hl

logger = logging.getLogger(__name__)

# Supported ID types for mapping
ID_TYPE = Literal[
    "hgnc_id",
    "gene_symbol",
    "ensembl_gene_id",
    "entrez_id",
    "uniprot_id",
]

# Fields that can be used as mapping sources (single-valued or first element of array)
MAPPABLE_FIELDS = {
    "hgnc_id": "hgnc_id",
    "gene_symbol": "gene_symbol",
    "ensembl_gene_id": "ensembl_gene_id",
    "entrez_id": "entrez_id",
    "uniprot_id": "uniprot_ids",  # Note: maps to array field
}


class GeneMapper:
    """
    Standalone utility for bidirectional gene ID mapping using HGNC data.

    This class provides methods to:
    - Map between different gene identifier types
    - Resolve gene symbol aliases to current approved symbols
    - Annotate Hail Tables with cross-reference IDs

    The mapper is designed as a standalone utility that can be composed
    with any hvantk component (streamers, tables, etc.).

    Example
    -------
    >>> hgnc_ht = hl.read_table("/tables/hgnc.ht")
    >>> mapper = GeneMapper(hgnc_ht)
    >>>
    >>> # Map Ensembl IDs to HGNC
    >>> mapper.map_to_hgnc(["ENSG00000012048"], source_type="ensembl_gene_id")
    {'ENSG00000012048': 'HGNC:1100'}
    >>>
    >>> # Resolve alias to current symbol
    >>> mapper.resolve_symbol("FANCD1")
    'BRCA2'
    """

    def __init__(self, hgnc_ht: hl.Table):
        """
        Initialize the mapper with an HGNC Hail Table.

        Parameters
        ----------
        hgnc_ht : hl.Table
            HGNC Hail Table built via ``hvantk reprocess hgnc:lookup`` and
            loaded via ``hl.read_table()``. Must be keyed by ``hgnc_id``.
        """
        self._ht = hgnc_ht
        self._validate_schema()

        # Lazy-loaded lookup dictionaries
        self._symbol_to_hgnc: Optional[Dict[str, str]] = None
        self._alias_to_hgnc: Optional[Dict[str, str]] = None
        self._ensembl_to_hgnc: Optional[Dict[str, str]] = None
        self._entrez_to_hgnc: Optional[Dict[str, str]] = None
        self._uniprot_to_hgnc: Optional[Dict[str, str]] = None
        self._hgnc_to_symbol: Optional[Dict[str, str]] = None

    def _validate_schema(self) -> None:
        """Validate that the HGNC table has required fields."""
        row_fields = set(self._ht.row)
        required = {"hgnc_id", "gene_symbol"}
        missing = required - row_fields
        if missing:
            raise ValueError(f"HGNC table missing required fields: {missing}")

        # Check key
        if list(self._ht.key) != ["hgnc_id"]:
            raise ValueError("HGNC table must be keyed by hgnc_id")

    def _build_symbol_lookup(self) -> None:
        """Build symbol -> hgnc_id lookup dictionary."""
        if self._symbol_to_hgnc is not None:
            return

        logger.info("Building symbol to HGNC ID lookup")
        data = self._ht.select("gene_symbol").collect()
        self._symbol_to_hgnc = {row.gene_symbol: row.hgnc_id for row in data}
        self._hgnc_to_symbol = {row.hgnc_id: row.gene_symbol for row in data}

    def _build_alias_lookup(self) -> None:
        """Build alias/prev_symbol -> hgnc_id lookup dictionary."""
        if self._alias_to_hgnc is not None:
            return

        logger.info("Building alias to HGNC ID lookup")
        row_fields = set(self._ht.row)

        # Collect alias and previous symbols
        self._alias_to_hgnc = {}

        if "alias_symbols" in row_fields:
            alias_ht = self._ht.select("alias_symbols")
            data = alias_ht.order_by(alias_ht.hgnc_id).collect()
            for row in data:
                if row.alias_symbols:
                    for alias in row.alias_symbols:
                        existing_hgnc = self._alias_to_hgnc.get(alias)
                        # Keep first mapping, but log collisions for visibility.
                        if existing_hgnc is None:
                            self._alias_to_hgnc[alias] = row.hgnc_id
                        elif existing_hgnc != row.hgnc_id:
                            logger.warning(
                                "Alias collision for %s: existing=%s new=%s",
                                alias,
                                existing_hgnc,
                                row.hgnc_id,
                            )

        if "prev_symbols" in row_fields:
            prev_ht = self._ht.select("prev_symbols")
            data = prev_ht.order_by(prev_ht.hgnc_id).collect()
            for row in data:
                if row.prev_symbols:
                    for prev in row.prev_symbols:
                        existing_hgnc = self._alias_to_hgnc.get(prev)
                        if existing_hgnc is None:
                            self._alias_to_hgnc[prev] = row.hgnc_id
                        elif existing_hgnc != row.hgnc_id:
                            logger.warning(
                                "Previous symbol collision for %s: existing=%s new=%s",
                                prev,
                                existing_hgnc,
                                row.hgnc_id,
                            )

    def _build_ensembl_lookup(self) -> None:
        """Build ensembl_gene_id -> hgnc_id lookup dictionary."""
        if self._ensembl_to_hgnc is not None:
            return

        logger.info("Building Ensembl to HGNC ID lookup")
        row_fields = set(self._ht.row)
        if "ensembl_gene_id" not in row_fields:
            self._ensembl_to_hgnc = {}
            return

        data = (
            self._ht.filter(
                hl.is_defined(self._ht.ensembl_gene_id)
                & (self._ht.ensembl_gene_id != "")
            )
            .select("ensembl_gene_id")
            .collect()
        )
        self._ensembl_to_hgnc = {row.ensembl_gene_id: row.hgnc_id for row in data}

    def _build_entrez_lookup(self) -> None:
        """Build entrez_id -> hgnc_id lookup dictionary."""
        if self._entrez_to_hgnc is not None:
            return

        logger.info("Building Entrez to HGNC ID lookup")
        row_fields = set(self._ht.row)
        if "entrez_id" not in row_fields:
            self._entrez_to_hgnc = {}
            return

        data = (
            self._ht.filter(
                hl.is_defined(self._ht.entrez_id) & (self._ht.entrez_id != "")
            )
            .select("entrez_id")
            .collect()
        )
        self._entrez_to_hgnc = {row.entrez_id: row.hgnc_id for row in data}

    def _build_uniprot_lookup(self) -> None:
        """Build uniprot_id -> hgnc_id lookup dictionary."""
        if self._uniprot_to_hgnc is not None:
            return

        logger.info("Building UniProt to HGNC ID lookup")
        row_fields = set(self._ht.row)
        if "uniprot_ids" not in row_fields:
            self._uniprot_to_hgnc = {}
            return

        uniprot_ht = self._ht.select("uniprot_ids")
        data = uniprot_ht.order_by(uniprot_ht.hgnc_id).collect()
        self._uniprot_to_hgnc = {}
        for row in data:
            if row.uniprot_ids:
                for uid in row.uniprot_ids:
                    existing_hgnc = self._uniprot_to_hgnc.get(uid)
                    if existing_hgnc is None:
                        self._uniprot_to_hgnc[uid] = row.hgnc_id
                    elif existing_hgnc != row.hgnc_id:
                        logger.warning(
                            "UniProt collision for %s: existing=%s new=%s",
                            uid,
                            existing_hgnc,
                            row.hgnc_id,
                        )

    # === Mapping Methods ===

    def map_to_hgnc(
        self,
        ids: List[str],
        source_type: ID_TYPE,
    ) -> Dict[str, Optional[str]]:
        """
        Map external IDs to HGNC IDs.

        Parameters
        ----------
        ids : list of str
            List of IDs to map.
        source_type : str
            Type of input IDs. One of: "gene_symbol", "ensembl_gene_id",
            "entrez_id", "uniprot_id".

        Returns
        -------
        dict
            Mapping from input ID to HGNC ID (or None if not found).
        """
        if source_type == "hgnc_id":
            # Identity mapping
            return {id_: id_ for id_ in ids}

        if source_type == "gene_symbol":
            self._build_symbol_lookup()
            return {id_: self._symbol_to_hgnc.get(id_) for id_ in ids}

        if source_type == "ensembl_gene_id":
            self._build_ensembl_lookup()
            return {id_: self._ensembl_to_hgnc.get(id_) for id_ in ids}

        if source_type == "entrez_id":
            self._build_entrez_lookup()
            return {id_: self._entrez_to_hgnc.get(id_) for id_ in ids}

        if source_type == "uniprot_id":
            self._build_uniprot_lookup()
            return {id_: self._uniprot_to_hgnc.get(id_) for id_ in ids}

        raise ValueError(f"Unknown source_type: {source_type}")

    def map_from_hgnc(
        self,
        hgnc_ids: List[str],
        target_type: ID_TYPE,
    ) -> Dict[str, Optional[str]]:
        """
        Map HGNC IDs to external IDs.

        Parameters
        ----------
        hgnc_ids : list of str
            List of HGNC IDs to map.
        target_type : str
            Type of output IDs. One of: "gene_symbol", "ensembl_gene_id",
            "entrez_id", "uniprot_id".

        Returns
        -------
        dict
            Mapping from HGNC ID to target ID (or None if not found).
            For uniprot_id, returns the first UniProt ID if multiple exist.
        """
        if target_type == "hgnc_id":
            return {id_: id_ for id_ in hgnc_ids}

        # Build lookup from HGNC table
        row_fields = set(self._ht.row)
        target_field = MAPPABLE_FIELDS.get(target_type, target_type)

        if target_field not in row_fields:
            logger.warning(f"Field {target_field} not in HGNC table")
            return {id_: None for id_ in hgnc_ids}

        # Filter to requested HGNC IDs and collect
        hgnc_set = hl.literal(set(hgnc_ids))
        data = (
            self._ht.filter(hgnc_set.contains(self._ht.hgnc_id))
            .select(target_field)
            .collect()
        )

        result = {}
        for row in data:
            value = getattr(row, target_field)
            if isinstance(value, list):
                # For array fields, return first element
                result[row.hgnc_id] = value[0] if value else None
            else:
                result[row.hgnc_id] = value if value else None

        # Fill in missing IDs with None
        for id_ in hgnc_ids:
            if id_ not in result:
                result[id_] = None

        return result

    def map_ids(
        self,
        ids: List[str],
        source_type: ID_TYPE,
        target_type: ID_TYPE,
    ) -> Dict[str, Optional[str]]:
        """
        Map between any two supported ID types (via HGNC as hub).

        Parameters
        ----------
        ids : list of str
            List of IDs to map.
        source_type : str
            Type of input IDs.
        target_type : str
            Type of output IDs.

        Returns
        -------
        dict
            Mapping from input ID to target ID (or None if not found).
        """
        if source_type == target_type:
            return {id_: id_ for id_ in ids}

        # Map source -> HGNC -> target
        to_hgnc = self.map_to_hgnc(ids, source_type)
        hgnc_ids = [v for v in to_hgnc.values() if v is not None]

        if not hgnc_ids:
            return {id_: None for id_ in ids}

        from_hgnc = self.map_from_hgnc(hgnc_ids, target_type)

        # Combine mappings
        result = {}
        for id_ in ids:
            hgnc_id = to_hgnc.get(id_)
            if hgnc_id is not None:
                result[id_] = from_hgnc.get(hgnc_id)
            else:
                result[id_] = None

        return result

    # === Symbol Resolution ===

    def resolve_symbol(self, symbol: str) -> Optional[str]:
        """
        Resolve a gene symbol (including aliases/previous) to the current approved symbol.

        Parameters
        ----------
        symbol : str
            Gene symbol to resolve (may be current, alias, or previous symbol).

        Returns
        -------
        str or None
            Current approved symbol, or None if symbol cannot be resolved.
        """
        self._build_symbol_lookup()

        # Check if it's already a current symbol
        if symbol in self._symbol_to_hgnc:
            return symbol

        # Check aliases
        self._build_alias_lookup()
        hgnc_id = self._alias_to_hgnc.get(symbol)
        if hgnc_id is not None:
            return self._hgnc_to_symbol.get(hgnc_id)

        return None

    def resolve_symbols(self, symbols: List[str]) -> Dict[str, Optional[str]]:
        """
        Batch resolve multiple symbols to current approved symbols.

        Parameters
        ----------
        symbols : list of str
            List of gene symbols to resolve.

        Returns
        -------
        dict
            Mapping from input symbol to current approved symbol (or None).
        """
        return {symbol: self.resolve_symbol(symbol) for symbol in symbols}

    # === Table Annotation ===

    def annotate_table(
        self,
        ht: hl.Table,
        source_field: str,
        source_type: ID_TYPE,
        fields_to_add: Optional[List[str]] = None,
    ) -> hl.Table:
        """
        Add HGNC annotations to a Hail Table by joining on a gene ID field.

        Parameters
        ----------
        ht : hl.Table
            Input Hail Table to annotate.
        source_field : str
            Field name in the input table containing gene IDs.
        source_type : str
            Type of ID in source_field (gene_symbol, ensembl_gene_id, etc.).
        fields_to_add : list of str, optional
            HGNC fields to add. If None, adds common fields:
            hgnc_id, gene_symbol, gene_name, locus_group.

        Returns
        -------
        hl.Table
            Annotated table with HGNC fields prefixed with "hgnc_".
        """
        if fields_to_add is None:
            fields_to_add = ["hgnc_id", "gene_symbol", "gene_name", "locus_group"]

        # Validate fields exist in HGNC table
        hgnc_fields = set(self._ht.row)
        invalid = set(fields_to_add) - hgnc_fields
        if invalid:
            raise ValueError(f"Fields not in HGNC table: {invalid}")

        # Prepare join key based on source_type
        if source_type == "hgnc_id":
            # Direct join on hgnc_id
            hgnc_select = self._ht.select(*fields_to_add)
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_select[ht[source_field]][f] for f in fields_to_add}
            )
        elif source_type == "ensembl_gene_id":
            # Re-key HGNC table by ensembl_gene_id for join
            hgnc_rekey = (
                self._ht.filter(
                    hl.is_defined(self._ht.ensembl_gene_id)
                    & (self._ht.ensembl_gene_id != "")
                )
                .key_by("ensembl_gene_id")
                .select(*fields_to_add)
            )
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_rekey[ht[source_field]][f] for f in fields_to_add}
            )
        elif source_type == "gene_symbol":
            # Re-key HGNC table by gene_symbol for join
            hgnc_rekey = self._ht.key_by("gene_symbol").select(*fields_to_add)
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_rekey[ht[source_field]][f] for f in fields_to_add}
            )
        elif source_type == "entrez_id":
            # Re-key HGNC table by entrez_id for join
            hgnc_rekey = (
                self._ht.filter(
                    hl.is_defined(self._ht.entrez_id) & (self._ht.entrez_id != "")
                )
                .key_by("entrez_id")
                .select(*fields_to_add)
            )
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_rekey[ht[source_field]][f] for f in fields_to_add}
            )
        else:
            raise ValueError(
                f"source_type '{source_type}' not supported for table annotation. "
                f"Use hgnc_id, gene_symbol, ensembl_gene_id, or entrez_id."
            )

        return ht

    # === Lookup Methods ===

    def get_gene_info(self, hgnc_id: str) -> Optional[Dict]:
        """
        Get full HGNC record for a gene.

        Parameters
        ----------
        hgnc_id : str
            HGNC ID to look up.

        Returns
        -------
        dict or None
            Dictionary of all HGNC fields for the gene, or None if not found.
        """
        rows = self._ht.filter(self._ht.hgnc_id == hgnc_id).collect()
        if not rows:
            return None
        row = rows[0]
        return {field: getattr(row, field) for field in self._ht.row}

    def get_genes_by_locus_group(
        self,
        locus_group: Literal[
            "protein-coding gene", "non-coding RNA", "pseudogene", "other"
        ],
    ) -> List[str]:
        """
        Get all HGNC IDs for genes in a locus group.

        Parameters
        ----------
        locus_group : str
            Locus group to filter by. Common values:
            "protein-coding gene", "non-coding RNA", "pseudogene", "other".

        Returns
        -------
        list of str
            List of HGNC IDs in the specified locus group.
        """
        row_fields = set(self._ht.row)
        if "locus_group" not in row_fields:
            logger.warning("locus_group field not in HGNC table")
            return []

        data = self._ht.filter(self._ht.locus_group == locus_group).collect()
        return [row.hgnc_id for row in data]

    def get_genes_by_gene_group(self, gene_group: str) -> List[str]:
        """
        Get all HGNC IDs for genes in a gene family/group.

        Parameters
        ----------
        gene_group : str
            Gene group/family name to search for.

        Returns
        -------
        list of str
            List of HGNC IDs in the specified gene group.
        """
        row_fields = set(self._ht.row)
        if "gene_group" not in row_fields:
            logger.warning("gene_group field not in HGNC table")
            return []

        # gene_group is an array field
        data = self._ht.filter(
            hl.any(lambda g: g == gene_group, self._ht.gene_group)
        ).collect()
        return [row.hgnc_id for row in data]

    # === Validation ===

    def validate_ids(self, ids: List[str], id_type: ID_TYPE) -> Dict[str, bool]:
        """
        Check which IDs are valid (exist in HGNC).

        Parameters
        ----------
        ids : list of str
            List of IDs to validate.
        id_type : str
            Type of IDs.

        Returns
        -------
        dict
            Mapping from ID to boolean (True if valid, False if not found).
        """
        mapping = self.map_to_hgnc(ids, id_type)
        return {id_: (hgnc_id is not None) for id_, hgnc_id in mapping.items()}

    # === Statistics ===

    def get_coverage_stats(self) -> Dict[str, int]:
        """
        Get statistics on ID coverage in the HGNC table.

        Returns
        -------
        dict
            Dictionary with counts:
            - total_genes: Total number of genes
            - with_ensembl: Genes with Ensembl ID
            - with_entrez: Genes with Entrez ID
            - with_uniprot: Genes with at least one UniProt ID
            - with_omim: Genes with OMIM ID
            - protein_coding: Protein-coding genes
        """
        row_fields = set(self._ht.row)

        stats = {"total_genes": self._ht.count()}

        if "ensembl_gene_id" in row_fields:
            stats["with_ensembl"] = self._ht.filter(
                hl.is_defined(self._ht.ensembl_gene_id)
                & (self._ht.ensembl_gene_id != "")
            ).count()

        if "entrez_id" in row_fields:
            stats["with_entrez"] = self._ht.filter(
                hl.is_defined(self._ht.entrez_id) & (self._ht.entrez_id != "")
            ).count()

        if "uniprot_ids" in row_fields:
            stats["with_uniprot"] = self._ht.filter(
                hl.is_defined(self._ht.uniprot_ids) & (hl.len(self._ht.uniprot_ids) > 0)
            ).count()

        if "omim_id" in row_fields:
            stats["with_omim"] = self._ht.filter(
                hl.is_defined(self._ht.omim_id) & (hl.len(self._ht.omim_id) > 0)
            ).count()

        if "locus_group" in row_fields:
            stats["protein_coding"] = self._ht.filter(
                self._ht.locus_group == "protein-coding gene"
            ).count()

        return stats
