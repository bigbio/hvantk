"""COSMIC Cancer Gene Census (CGC) data streamer.

Subclass of :class:`GeneDiseaseValidityStreamer` that adds COSMIC-specific
somatic/germline filtering, tumour-type and role-based gene set extraction.

CGC tables have one row per gene (not per gene-disease assertion), so methods
that depend on MONDO ontology IDs raise ``NotImplementedError``.
"""

from __future__ import annotations

import logging
from typing import Dict, Optional, Set

import hail as hl
import pandas as pd

from hvantk.core.constants import (
    COSMIC_CGC_CLASSIFICATION_LEVELS,
    COSMIC_MUTATION_CONTEXTS,
)
from hvantk.data.gene_disease_streamer import GeneDiseaseValidityStreamer

logger = logging.getLogger(__name__)


class CosmicCGCStreamer(GeneDiseaseValidityStreamer):
    """Streamer for querying and aggregating COSMIC Cancer Gene Census data.

    Inherits query, aggregation, and export machinery from
    :class:`GeneDiseaseValidityStreamer` in gene-keyed mode.  Adds
    COSMIC-specific methods for somatic/germline filtering and gene set
    extraction by tumour type, role in cancer, and tissue type.
    """

    source_name = "COSMIC CGC"
    classification_levels = COSMIC_CGC_CLASSIFICATION_LEVELS
    annotation_prefix = "cosmic_cgc_"

    # ------------------------------------------------------------------
    # Extension points
    # ------------------------------------------------------------------

    @property
    def _gene_disease_required_fields(self) -> Set[str]:
        return {"somatic", "germline", "role_in_cancer", "tissue_type"}

    @property
    def _grouping_field(self) -> str:
        return "tissue_type"

    @property
    def _grouping_field_short(self) -> str:
        return "tissue"

    @property
    def _date_field(self) -> Optional[str]:
        return None

    # ------------------------------------------------------------------
    # Keying and validation
    # ------------------------------------------------------------------

    def _detect_keying_mode(self, ht: hl.Table) -> str:
        key_fields = set(ht.key)
        if "hgnc_id" in key_fields:
            return "gene"
        if "gene_symbol" in key_fields:
            return "gene"
        return "unknown"

    def validate_table(self) -> bool:
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        row_fields = set(ht.row)
        required = {
            "gene_symbol",
            "classification",
            "classification_level",
            "somatic",
            "germline",
        }
        missing = sorted(required - row_fields)
        if missing:
            raise ValueError(
                f"COSMIC CGC table missing required fields: " f"{', '.join(missing)}"
            )
        return True

    # ------------------------------------------------------------------
    # Classification filter override
    # ------------------------------------------------------------------

    def _apply_min_classification_filter(
        self, ht: hl.Table, min_classification: str
    ) -> hl.Table:
        """CGC tables use classification_level directly (one tier per gene)."""
        min_level = self._classification_level(min_classification)
        cache_key = None
        if ht is self._table:
            cache_key = f"min:{self._keying_mode}:{min_classification}"
            cached = self._filtered_cache.get(cache_key)
            if cached is not None:
                return cached
        filtered = ht.filter(ht.classification_level <= min_level)
        if cache_key:
            self._filtered_cache[cache_key] = filtered
        return filtered

    # ------------------------------------------------------------------
    # Somatic / Germline filtering
    # ------------------------------------------------------------------

    def filter_by_mutation_context(
        self,
        context: str = "both",
    ) -> hl.Table:
        """Return a filtered view of the table by mutation context.

        Parameters
        ----------
        context : str
            ``"somatic"`` keeps genes with somatic==True,
            ``"germline"`` keeps genes with germline==True,
            ``"both"`` returns all genes (no filter).

        Returns
        -------
        hl.Table
        """
        if context not in COSMIC_MUTATION_CONTEXTS:
            raise ValueError(
                f"context must be one of {COSMIC_MUTATION_CONTEXTS}, " f"got: {context}"
            )
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        if context == "somatic":
            return ht.filter(ht.somatic)
        elif context == "germline":
            return ht.filter(ht.germline)
        return ht

    # ------------------------------------------------------------------
    # Gene set extraction
    # ------------------------------------------------------------------

    def get_geneset_per_tumour_type(
        self,
        min_classification: Optional[str] = None,
        min_genes: int = 0,
        mutation_context: str = "both",
    ) -> Dict[str, Set[str]]:
        """Gene sets grouped by tumour type.

        Explodes the multi-value tumour_types field so each tumour type
        becomes a gene set key.  Respects mutation_context to select
        somatic, germline, or merged tumour types.
        """
        self._ensure_table_loaded()
        ht = self.filter_by_mutation_context(mutation_context)

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        # Select the tumour types field(s) based on context
        row_fields = set(ht.row)
        if mutation_context == "somatic" and "tumour_types_somatic" in row_fields:
            ht = ht.annotate(_tumour_types=ht.tumour_types_somatic)
        elif mutation_context == "germline" and "tumour_types_germline" in row_fields:
            ht = ht.annotate(_tumour_types=ht.tumour_types_germline)
        else:
            # Merge both somatic and germline tumour types
            somatic = (
                ht.tumour_types_somatic
                if "tumour_types_somatic" in row_fields
                else hl.empty_array(hl.tstr)
            )
            germline = (
                ht.tumour_types_germline
                if "tumour_types_germline" in row_fields
                else hl.empty_array(hl.tstr)
            )
            ht = ht.annotate(_tumour_types=hl.array(hl.set(somatic.extend(germline))))

        # Explode and group
        exploded = ht.explode("_tumour_types")
        exploded = exploded.filter(
            hl.is_defined(exploded._tumour_types) & (exploded._tumour_types != "")
        )
        grouped = exploded.group_by(tumour_type=exploded._tumour_types).aggregate(
            genes=hl.agg.collect_as_set(exploded.gene_symbol)
        )
        rows = grouped.collect()

        result: Dict[str, Set[str]] = {}
        for row in rows:
            if not row.tumour_type:
                continue
            genes = set(row.genes)
            if len(genes) < min_genes:
                continue
            result[row.tumour_type] = genes

        logger.info(
            f"Built {len(result)} tumour-type gene sets"
            + (f" (min_genes={min_genes})" if min_genes > 0 else "")
        )
        return result

    def get_geneset_per_role(
        self,
        min_classification: Optional[str] = None,
        min_genes: int = 0,
        mutation_context: str = "both",
    ) -> Dict[str, Set[str]]:
        """Gene sets grouped by Role in Cancer (oncogene, TSG, fusion)."""
        self._ensure_table_loaded()
        ht = self.filter_by_mutation_context(mutation_context)

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        row_fields = set(ht.row)
        if "role_in_cancer" not in row_fields:
            logger.warning("role_in_cancer field not found in table")
            return {}

        exploded = ht.explode("role_in_cancer")
        exploded = exploded.filter(
            hl.is_defined(exploded.role_in_cancer) & (exploded.role_in_cancer != "")
        )
        grouped = exploded.group_by(role=exploded.role_in_cancer).aggregate(
            genes=hl.agg.collect_as_set(exploded.gene_symbol)
        )
        rows = grouped.collect()

        result: Dict[str, Set[str]] = {}
        for row in rows:
            if not row.role:
                continue
            genes = set(row.genes)
            if len(genes) < min_genes:
                continue
            result[row.role] = genes

        logger.info(
            f"Built {len(result)} role-based gene sets"
            + (f" (min_genes={min_genes})" if min_genes > 0 else "")
        )
        return result

    def get_geneset_per_tissue(
        self,
        min_classification: Optional[str] = None,
        min_genes: int = 0,
        mutation_context: str = "both",
    ) -> Dict[str, Set[str]]:
        """Gene sets grouped by tissue type (E, L, M, O)."""
        self._ensure_table_loaded()
        ht = self.filter_by_mutation_context(mutation_context)

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        row_fields = set(ht.row)
        if "tissue_type" not in row_fields:
            logger.warning("tissue_type field not found in table")
            return {}

        grouped = ht.group_by(tissue=ht.tissue_type).aggregate(
            genes=hl.agg.collect_as_set(ht.gene_symbol)
        )
        rows = grouped.collect()

        result: Dict[str, Set[str]] = {}
        for row in rows:
            tissue = row.tissue
            if not tissue or tissue.strip() == "":
                continue
            genes = set(row.genes)
            if len(genes) < min_genes:
                continue
            result[tissue] = genes

        logger.info(
            f"Built {len(result)} tissue-type gene sets"
            + (f" (min_genes={min_genes})" if min_genes > 0 else "")
        )
        return result

    # ------------------------------------------------------------------
    # Summary methods
    # ------------------------------------------------------------------

    def mutation_context_summary(self) -> pd.DataFrame:
        """Breakdown of genes by somatic-only, germline-only, and both."""
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        rows = ht.select("gene_symbol", "somatic", "germline").collect()
        somatic_only = set()
        germline_only = set()
        both = set()
        for row in rows:
            s = row.somatic
            g = row.germline
            if s and g:
                both.add(row.gene_symbol)
            elif s:
                somatic_only.add(row.gene_symbol)
            elif g:
                germline_only.add(row.gene_symbol)

        return pd.DataFrame(
            [
                {"context": "somatic_only", "n_genes": len(somatic_only)},
                {"context": "germline_only", "n_genes": len(germline_only)},
                {"context": "both", "n_genes": len(both)},
            ]
        )

    def role_summary(self) -> pd.DataFrame:
        """Breakdown of genes by oncogene / TSG / fusion role."""
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        row_fields = set(ht.row)
        if "role_in_cancer" not in row_fields:
            return pd.DataFrame(columns=["role", "n_genes"])

        exploded = ht.explode("role_in_cancer")
        exploded = exploded.filter(
            hl.is_defined(exploded.role_in_cancer) & (exploded.role_in_cancer != "")
        )
        grouped = exploded.group_by(role=exploded.role_in_cancer).aggregate(
            n_genes=hl.agg.count()
        )
        return (
            grouped.to_pandas()
            .sort_values("n_genes", ascending=False)
            .reset_index(drop=True)
        )

    # ------------------------------------------------------------------
    # Inherited methods that work as-is
    # ------------------------------------------------------------------
    # get_genes_by_classification(), to_gene_set(), annotate_gene_table(),
    # export_for_enrichex(), _normalize_classification(), setup(), stream()

    # ------------------------------------------------------------------
    # Methods that don't apply (MONDO-dependent)
    # ------------------------------------------------------------------

    def get_genes_by_mondo_id(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not use MONDO ontology IDs. "
            "Use filter_by_mutation_context() or get_geneset_per_tumour_type() instead."
        )

    def categorize_by_ontology(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not use MONDO ontology IDs. "
            "Use get_geneset_per_tumour_type() or get_geneset_per_role() instead."
        )

    def get_genes_by_ontology_category(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not use MONDO ontology IDs. "
            "Use get_geneset_per_tumour_type() or get_geneset_per_role() instead."
        )

    def categorize_by_ontology_summary(self, *args, **kwargs):
        raise NotImplementedError("COSMIC CGC does not use MONDO ontology IDs.")

    # Override base methods that assume disease_label / mode_of_inheritance
    def get_genes_by_disease(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not have per-disease assertions. "
            "Use get_geneset_per_tumour_type() instead."
        )

    def get_genes_by_moi(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not include mode of inheritance data."
        )

    def get_geneset_per_disease(self, *args, **kwargs):
        raise NotImplementedError(
            "COSMIC CGC does not have per-disease assertions. "
            "Use get_geneset_per_tumour_type() instead."
        )

    def aggregate_by_disease_category(
        self,
        categories,
        min_classification=None,
    ):
        """Aggregate genes into keyword categories by matching tumour types.

        Since COSMIC CGC lacks MONDO disease labels, this searches the
        merged tumour types (somatic + germline) for keyword matches.
        """
        from typing import Dict, List, Set

        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        row_fields = set(ht.row)
        somatic = (
            ht.tumour_types_somatic
            if "tumour_types_somatic" in row_fields
            else hl.empty_array(hl.tstr)
        )
        germline = (
            ht.tumour_types_germline
            if "tumour_types_germline" in row_fields
            else hl.empty_array(hl.tstr)
        )
        ht = ht.annotate(_all_tumour_types=hl.array(hl.set(somatic.extend(germline))))
        rows = ht.select("gene_symbol", "_all_tumour_types").collect()

        gene_sets: Dict[str, Set[str]] = {}
        for cat_name, terms in categories.items():
            matched_genes: Set[str] = set()
            lower_terms = [t.lower() for t in terms]
            for row in rows:
                for tt in row._all_tumour_types or []:
                    if any(term in tt.lower() for term in lower_terms):
                        matched_genes.add(row.gene_symbol)
                        break
            gene_sets[cat_name] = matched_genes

        return gene_sets

    def compute_stats(self) -> Dict:
        """Compute COSMIC CGC-specific statistics."""
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        total = ht.count()
        classification_counts = ht.aggregate(hl.agg.counter(ht.classification))

        return {
            "total_genes": total,
            "classification_counts": classification_counts,
        }

    def classification_summary(self) -> pd.DataFrame:
        """Summary of genes by tier classification."""
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("COSMIC CGC table not loaded.")

        grouped = ht.group_by(
            classification=ht.classification,
            classification_level=ht.classification_level,
        ).aggregate(
            n_genes=hl.agg.count(),
            genes=hl.agg.collect_as_set(ht.gene_symbol),
        )

        df = grouped.to_pandas()
        total = ht.count()
        df["pct_genes"] = (df["n_genes"] / total) * 100
        df = df.sort_values("classification_level").reset_index(drop=True)
        return df.drop(columns=["classification_level", "genes"])

    def describe(self) -> str:
        """Return a human-readable description of the dataset."""
        stats = self.compute_stats()
        total = stats["total_genes"]
        lines = [
            "COSMIC Cancer Gene Census Dataset",
            "=" * 34,
            f"Total genes: {total:,}",
            "",
            "Tier distribution:",
        ]
        for level in self.classification_levels:
            count = stats["classification_counts"].get(level, 0)
            pct = (count / total * 100) if total else 0.0
            lines.append(f"  {level}: {count:,} ({pct:.1f}%)")
        return "\n".join(lines)
