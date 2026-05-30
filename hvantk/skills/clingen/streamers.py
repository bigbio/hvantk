"""ClinGen Gene-Disease Validity data streamer.

Thin subclass of :class:`GeneDiseaseTableStreamer` that adds
ClinGen-specific GCEP (Gene Curation Expert Panel) methods.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Set

if TYPE_CHECKING:
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

import hail as hl
import pandas as pd

from hvantk.skills.clingen.shared.constants import CLINGEN_CLASSIFICATION_LEVELS
from hvantk.core.streamers.gene_disease_table import GeneDiseaseTableStreamer

logger = logging.getLogger(__name__)


class ClinGenGeneDiseaseTableStreamer(GeneDiseaseTableStreamer):
    """Streamer for querying and aggregating ClinGen Gene-Disease Validity data.

    Inherits all generic query, aggregation, ontology, integration, and
    summary methods from :class:`GeneDiseaseTableStreamer`.  Adds
    ClinGen-specific methods for working with Gene Curation Expert Panels.
    """

    source_name = "ClinGen"
    classification_levels = CLINGEN_CLASSIFICATION_LEVELS
    annotation_prefix = "clingen_"

    # ------------------------------------------------------------------
    # ClinGen-specific extension points
    # ------------------------------------------------------------------

    @property
    def _gene_disease_required_fields(self) -> Set[str]:
        return {"classification_date", "gene_curation_expert_panel"}

    @property
    def _grouping_field(self) -> str:
        return "gene_curation_expert_panel"

    @property
    def _grouping_field_short(self) -> str:
        return "gcep"

    @property
    def _date_field(self) -> Optional[str]:
        return "classification_date"

    # ------------------------------------------------------------------
    # GCEP-specific methods
    # ------------------------------------------------------------------

    def gcep_summary(self) -> pd.DataFrame:
        """Generate a summary table by Gene Curation Expert Panel.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: gcep, n_associations, n_genes,
            n_diseases, top_classification, last_curation.
        """
        self._ensure_gene_disease_mode("gcep_summary")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        summary_ht = ht.group_by(gcep=ht.gene_curation_expert_panel).aggregate(
            n_associations=hl.agg.count(),
            genes=hl.agg.collect_as_set(ht.gene_symbol),
            diseases=hl.agg.collect_as_set(ht.disease_label),
            classification_counts=hl.agg.counter(ht.classification),
            last_curation=hl.agg.max(ht.classification_date),
        )

        df = summary_ht.to_pandas()
        df["n_genes"] = df["genes"].apply(len)
        df["n_diseases"] = df["diseases"].apply(len)
        df["top_classification"] = df["classification_counts"].apply(
            lambda counts: max(counts, key=counts.get) if counts else None
        )
        return df.drop(
            columns=["classification_counts", "genes", "diseases"]
        ).sort_values("n_associations", ascending=False)

    def get_geneset_per_gcep(
        self,
        min_classification: Optional[str] = None,
        min_genes: int = 0,
        shorten_names: bool = True,
    ) -> Dict[str, Set[str]]:
        """Get genesets grouped by Gene Curation Expert Panel (GCEP).

        Parameters
        ----------
        min_classification : str, optional
            Filter to minimum classification level before grouping.
        min_genes : int
            Exclude GCEPs with fewer than this many genes (default: 0).
        shorten_names : bool
            If True (default), remove the " Gene Curation Expert Panel"
            suffix from GCEP names.

        Returns
        -------
        dict
            Mapping of GCEP name -> set of gene symbols.
        """
        self._ensure_gene_disease_mode("get_geneset_per_gcep")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        grouped = ht.group_by(gcep=ht.gene_curation_expert_panel).aggregate(
            genes=hl.agg.collect_as_set(ht.gene_symbol)
        )
        rows = grouped.collect()

        result: Dict[str, Set[str]] = {}
        for row in rows:
            gcep_name = row.gcep
            if not gcep_name:
                continue
            genes = set(row.genes)
            if len(genes) < min_genes:
                continue
            if shorten_names:
                gcep_name = gcep_name.replace(" Gene Curation Expert Panel", "")
            result[gcep_name] = genes

        logger.info(
            f"Built {len(result)} GCEP-based gene sets"
            + (f" (min_genes={min_genes})" if min_genes > 0 else "")
        )
        return result

    def group_by_gcep(
        self,
        min_classification: Optional[str] = None,
    ) -> Dict[str, hl.Table]:
        """Group data by Gene Curation Expert Panel.

        Parameters
        ----------
        min_classification : str, optional
            Filter to minimum classification level.

        Returns
        -------
        dict
            Mapping of GCEP name to Hail Table subset.
        """
        self._ensure_gene_disease_mode("group_by_gcep")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        gceps = sorted(
            ht.aggregate(hl.agg.collect_as_set(ht.gene_curation_expert_panel))
        )
        return {
            gcep: ht.filter(ht.gene_curation_expert_panel == gcep) for gcep in gceps
        }
