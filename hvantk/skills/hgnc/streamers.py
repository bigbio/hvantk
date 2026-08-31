"""HGNC concrete streamers.

Implements ``GeneCatalogStreamer`` for HGNC gene metadata. Absorbs the
logic previously in ``core/utils/gene_mapper.py`` and
``core/utils/gene_aliases.py``.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Literal, Optional

import hail as hl

from hvantk.core.models import AnnotationTable
from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

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


class HGNCGeneCatalogStreamer(GeneCatalogStreamer):
    """HGNC gene catalog streamer.

    Implements GeneCatalogStreamer over an HGNC AnnotationTable. The
    underlying schema includes hgnc_id, symbol, alias_symbols,
    prev_symbols, ensembl_gene_id, entrez_id, uniprot_ids, locus_group,
    gene_group, and other fields documented in the HGNC plugin.

    Absorbs the logic previously in ``core/utils/gene_mapper.py``
    (``GeneMapper``) and ``core/utils/gene_aliases.py``
    (``_load_hgnc_symbol_maps``, ``expand_gene_set_with_aliases``).
    """

    def __init__(self, artifact: AnnotationTable) -> None:
        super().__init__(artifact)
        self._ht: hl.Table = artifact.to_hail()
        self._validate_schema()

        # Lazy-loaded lookup dictionaries (absorbed from GeneMapper)
        self._symbol_to_hgnc: Optional[Dict[str, str]] = None
        self._alias_to_hgnc: Optional[Dict[str, str]] = None
        self._ensembl_to_hgnc: Optional[Dict[str, str]] = None
        self._entrez_to_hgnc: Optional[Dict[str, str]] = None
        self._uniprot_to_hgnc: Optional[Dict[str, str]] = None
        self._hgnc_to_symbol: Optional[Dict[str, str]] = None

        # Build symbol maps eagerly (absorbed from gene_aliases._load_hgnc_symbol_maps)
        self._canonical_symbols, self._alias_to_canonical = self._build_symbol_maps()

    def _validate_schema(self) -> None:
        """Validate that the HGNC table has required fields."""
        row_fields = set(self._ht.row)
        required = {"hgnc_id", "gene_symbol"}
        missing = required - row_fields
        if missing:
            raise ValueError(f"HGNC table missing required fields: {missing}")
        if list(self._ht.key) != ["hgnc_id"]:
            raise ValueError("HGNC table must be keyed by hgnc_id")

    def _build_symbol_maps(self) -> tuple[set[str], dict[str, str]]:
        """Build canonical-symbol set and alias-to-canonical dict.

        Absorbed from ``core/utils/gene_aliases._load_hgnc_symbol_maps``
        (Hail Table path only; the TSV path is no longer needed since the
        artifact layer always wraps a Hail Table at this point).
        """
        row_fields = set(self._ht.row)

        fields_to_select = ["gene_symbol"]
        if "alias_symbols" in row_fields:
            fields_to_select.append("alias_symbols")
        if "prev_symbols" in row_fields:
            fields_to_select.append("prev_symbols")

        data = self._ht.select(*fields_to_select).collect()

        canonical_symbols: set[str] = set()
        alias_to_canonical: dict[str, str] = {}

        # Pass 1: collect all canonical symbols
        for row in data:
            canonical_symbols.add(row.gene_symbol)

        # Pass 2: build alias mappings (canonical set must be complete first)
        for row in data:
            symbol = row.gene_symbol

            if "alias_symbols" in fields_to_select and row.alias_symbols:
                for alias in row.alias_symbols:
                    if alias and alias not in canonical_symbols:
                        if alias not in alias_to_canonical:
                            alias_to_canonical[alias] = symbol

            if "prev_symbols" in fields_to_select and row.prev_symbols:
                for prev in row.prev_symbols:
                    if prev and prev not in canonical_symbols:
                        if prev not in alias_to_canonical:
                            alias_to_canonical[prev] = symbol

        return canonical_symbols, alias_to_canonical

    # ---- Lazy lookup builders (absorbed from GeneMapper) ----

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
        self._alias_to_hgnc = {}

        if "alias_symbols" in row_fields:
            alias_ht = self._ht.select("alias_symbols")
            data = alias_ht.order_by(alias_ht.hgnc_id).collect()
            for row in data:
                if row.alias_symbols:
                    for alias in row.alias_symbols:
                        existing_hgnc = self._alias_to_hgnc.get(alias)
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

    # ---- ABC contract methods ----

    def is_canonical(self, symbol: str) -> bool:
        return symbol in self._canonical_symbols

    def resolve_alias(self, symbol: str) -> str | None:
        return self._alias_to_canonical.get(symbol)

    def expand_with_aliases(self, symbols: set[str]) -> tuple[set[str], dict[str, str]]:
        """Expand a gene set to include all known HGNC aliases.

        Absorbed from ``core/utils/gene_aliases.expand_gene_set_with_aliases``.

        For each provided gene symbol:
        - If it is a canonical symbol, also include its known aliases/prev_symbols
        - If it is an alias/prev_symbol, also include the canonical symbol

        Parameters
        ----------
        symbols : set of str
            Input gene symbols (canonical or aliases).

        Returns
        -------
        (expanded_set, alias_map)
            expanded_set contains the input plus any resolved canonical/alias forms;
            alias_map maps each alias that was resolved to its canonical symbol.
        """
        # Build the canonical_to_aliases reverse map on the fly
        row_fields = set(self._ht.row)
        fields_to_select = ["gene_symbol"]
        if "alias_symbols" in row_fields:
            fields_to_select.append("alias_symbols")
        if "prev_symbols" in row_fields:
            fields_to_select.append("prev_symbols")

        data = self._ht.select(*fields_to_select).collect()
        canonical_to_aliases: dict[str, list[str]] = {}
        for row in data:
            symbol = row.gene_symbol
            all_aliases: list[str] = []
            if "alias_symbols" in fields_to_select and row.alias_symbols:
                for alias in row.alias_symbols:
                    if alias and alias not in self._canonical_symbols:
                        if alias in self._alias_to_canonical:
                            all_aliases.append(alias)
            if "prev_symbols" in fields_to_select and row.prev_symbols:
                for prev in row.prev_symbols:
                    if prev and prev not in self._canonical_symbols:
                        if prev in self._alias_to_canonical:
                            all_aliases.append(prev)
            if all_aliases:
                canonical_to_aliases[symbol] = all_aliases

        expanded: set[str] = set(symbols)
        alias_map: dict[str, str] = {}

        for gene in symbols:
            if gene in self._canonical_symbols:
                # User provided a canonical symbol -> add its aliases
                aliases = canonical_to_aliases.get(gene, [])
                for alias in aliases:
                    if alias not in expanded:
                        expanded.add(alias)
                        alias_map[alias] = gene
            elif gene in self._alias_to_canonical:
                # User provided an alias -> add the canonical symbol
                canonical = self._alias_to_canonical[gene]
                alias_map[gene] = canonical
                if canonical not in expanded:
                    expanded.add(canonical)
                # Also add other aliases of the same canonical symbol
                for alias in canonical_to_aliases.get(canonical, []):
                    if alias not in expanded:
                        expanded.add(alias)
                        alias_map[alias] = canonical

        return expanded, alias_map

    def map_ids(
        self,
        ids: list[str],
        source_type: str,
        target_type: str,
    ) -> dict[str, str | None]:
        """Map between any two supported ID types (via HGNC as hub).

        Absorbed from ``core/utils/gene_mapper.GeneMapper.map_ids``.

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

        result = {}
        for id_ in ids:
            hgnc_id = to_hgnc.get(id_)
            if hgnc_id is not None:
                result[id_] = from_hgnc.get(hgnc_id)
            else:
                result[id_] = None

        return result

    # ---- HGNC-specific extension methods (not on the ABC) ----

    def map_to_hgnc(
        self,
        ids: List[str],
        source_type: ID_TYPE,
    ) -> Dict[str, Optional[str]]:
        """Map external IDs to HGNC IDs.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.map_to_hgnc``.

        Parameters
        ----------
        ids : list of str
            List of IDs to map.
        source_type : str
            Type of input IDs. One of: "hgnc_id", "gene_symbol",
            "ensembl_gene_id", "entrez_id", "uniprot_id".

        Returns
        -------
        dict
            Mapping from input ID to HGNC ID (or None if not found).
        """
        if source_type == "hgnc_id":
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
        """Map HGNC IDs to external IDs.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.map_from_hgnc``.

        Parameters
        ----------
        hgnc_ids : list of str
            List of HGNC IDs to map.
        target_type : str
            Type of output IDs. One of: "hgnc_id", "gene_symbol",
            "ensembl_gene_id", "entrez_id", "uniprot_id".

        Returns
        -------
        dict
            Mapping from HGNC ID to target ID (or None if not found).
            For uniprot_id, returns the first UniProt ID if multiple exist.
        """
        if not hgnc_ids:
            return {}

        if target_type == "hgnc_id":
            return {id_: id_ for id_ in hgnc_ids}

        row_fields = set(self._ht.row)
        target_field = MAPPABLE_FIELDS.get(target_type, target_type)

        if target_field not in row_fields:
            logger.warning(f"Field {target_field} not in HGNC table")
            return {id_: None for id_ in hgnc_ids}

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
                result[row.hgnc_id] = value[0] if value else None
            else:
                result[row.hgnc_id] = value if value else None

        for id_ in hgnc_ids:
            if id_ not in result:
                result[id_] = None

        return result

    def resolve_symbol(self, symbol: str) -> Optional[str]:
        """Resolve a gene symbol (including aliases/previous) to the current
        approved symbol.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.resolve_symbol``.

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

        if symbol in self._symbol_to_hgnc:
            return symbol

        self._build_alias_lookup()
        hgnc_id = self._alias_to_hgnc.get(symbol)
        if hgnc_id is not None:
            return self._hgnc_to_symbol.get(hgnc_id)

        return None

    def resolve_symbols(self, symbols: List[str]) -> Dict[str, Optional[str]]:
        """Batch resolve multiple symbols to current approved symbols.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.resolve_symbols``.
        """
        return {symbol: self.resolve_symbol(symbol) for symbol in symbols}

    def resolve_to_canonical(self, symbol: str) -> str:
        """Convenience: resolve_alias or return input unchanged.

        Preserved for builders that called GeneMapper.resolve_to_canonical.
        """
        return self._alias_to_canonical.get(symbol, symbol)

    def annotate_table(
        self,
        ht: hl.Table,
        source_field: str,
        source_type: ID_TYPE,
        fields_to_add: Optional[List[str]] = None,
    ) -> hl.Table:
        """Add HGNC annotations to a Hail Table by joining on a gene ID field.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.annotate_table``.

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

        hgnc_fields = set(self._ht.row)
        invalid = set(fields_to_add) - hgnc_fields
        if invalid:
            raise ValueError(f"Fields not in HGNC table: {invalid}")

        if source_type == "hgnc_id":
            hgnc_select = self._ht.select(*fields_to_add)
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_select[ht[source_field]][f] for f in fields_to_add}
            )
        elif source_type == "ensembl_gene_id":
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
            hgnc_rekey = self._ht.key_by("gene_symbol").select(*fields_to_add)
            ht = ht.annotate(
                **{f"hgnc_{f}": hgnc_rekey[ht[source_field]][f] for f in fields_to_add}
            )
        elif source_type == "entrez_id":
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

    def get_gene_info(self, hgnc_id: str) -> Optional[Dict]:
        """Get full HGNC record for a gene.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.get_gene_info``.
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
        """Get all HGNC IDs for genes in a locus group.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.get_genes_by_locus_group``.
        """
        row_fields = set(self._ht.row)
        if "locus_group" not in row_fields:
            logger.warning("locus_group field not in HGNC table")
            return []
        data = self._ht.filter(self._ht.locus_group == locus_group).collect()
        return [row.hgnc_id for row in data]

    def get_genes_by_gene_group(self, gene_group: str) -> List[str]:
        """Get all HGNC IDs for genes in a gene family/group.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.get_genes_by_gene_group``.
        """
        row_fields = set(self._ht.row)
        if "gene_group" not in row_fields:
            logger.warning("gene_group field not in HGNC table")
            return []
        data = self._ht.filter(
            hl.any(lambda g: g == gene_group, self._ht.gene_group)
        ).collect()
        return [row.hgnc_id for row in data]

    def validate_ids(self, ids: List[str], id_type: ID_TYPE) -> Dict[str, bool]:
        """Check which IDs are valid (exist in HGNC).

        Absorbed from ``core/utils/gene_mapper.GeneMapper.validate_ids``.
        """
        mapping = self.map_to_hgnc(ids, id_type)
        return {id_: (hgnc_id is not None) for id_, hgnc_id in mapping.items()}

    def get_coverage_stats(self) -> Dict[str, int]:
        """Get statistics on ID coverage in the HGNC table.

        Absorbed from ``core/utils/gene_mapper.GeneMapper.get_coverage_stats``.
        """
        ht = self._ht
        row_fields = set(ht.row)

        # Fuse total + per-field coverage counts into a single Hail action
        # (one Spark job) instead of one ``filter(...).count()`` per field.
        agg_exprs = {"total_genes": hl.agg.count()}

        if "ensembl_gene_id" in row_fields:
            agg_exprs["with_ensembl"] = hl.agg.count_where(
                hl.is_defined(ht.ensembl_gene_id) & (ht.ensembl_gene_id != "")
            )

        if "entrez_id" in row_fields:
            agg_exprs["with_entrez"] = hl.agg.count_where(
                hl.is_defined(ht.entrez_id) & (ht.entrez_id != "")
            )

        if "uniprot_ids" in row_fields:
            agg_exprs["with_uniprot"] = hl.agg.count_where(
                hl.is_defined(ht.uniprot_ids) & (hl.len(ht.uniprot_ids) > 0)
            )

        if "omim_id" in row_fields:
            agg_exprs["with_omim"] = hl.agg.count_where(
                hl.is_defined(ht.omim_id) & (hl.len(ht.omim_id) > 0)
            )

        if "locus_group" in row_fields:
            agg_exprs["protein_coding"] = hl.agg.count_where(
                ht.locus_group == "protein-coding gene"
            )

        agg = ht.aggregate(hl.struct(**agg_exprs))
        return {key: agg[key] for key in agg_exprs}
