"""GenCC (Gene Curation Coalition) submissions data streamer.

Thin subclass of :class:`GeneDiseaseValidityStreamer` that adds
GenCC-specific submitter methods and the ``gene_disease_submitter``
keying mode.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Union

if TYPE_CHECKING:
    from hvantk.data.gene_mapper import GeneMapper

import hail as hl
import pandas as pd

from hvantk.core.constants import GENCC_CLASSIFICATION_LEVELS
from hvantk.data.gene_disease_streamer import GeneDiseaseValidityStreamer

logger = logging.getLogger(__name__)


class GenCCStreamer(GeneDiseaseValidityStreamer):
    """Streamer for querying and aggregating GenCC submissions data.

    Inherits all generic query, aggregation, ontology, integration, and
    summary methods from :class:`GeneDiseaseValidityStreamer`.  Adds
    GenCC-specific methods for working with submitting organizations.
    """

    source_name = "GenCC"
    classification_levels = GENCC_CLASSIFICATION_LEVELS
    annotation_prefix = "gencc_"

    # ------------------------------------------------------------------
    # GenCC-specific extension points
    # ------------------------------------------------------------------

    @property
    def _gene_disease_required_fields(self) -> Set[str]:
        return {"submission_date", "submitter"}

    @property
    def _grouping_field(self) -> str:
        return "submitter"

    @property
    def _grouping_field_short(self) -> str:
        return "submitter"

    @property
    def _date_field(self) -> Optional[str]:
        return "submission_date"

    def _detect_keying_mode(self, ht: hl.Table) -> str:
        """Extend base detection with gene_disease_submitter mode."""
        key_fields = set(ht.key)
        if {"hgnc_id", "mondo_id", "submitter"} <= key_fields:
            return "gene_disease_submitter"
        return super()._detect_keying_mode(ht)

    # ------------------------------------------------------------------
    # Submitter-specific methods
    # ------------------------------------------------------------------

    def get_genes_by_submitter(
        self,
        submitter: str,
        min_classification: Optional[str] = None,
        as_set: bool = True,
    ) -> Union[Set[str], hl.Table]:
        """Get genes asserted by a specific submitting organization.

        Parameters
        ----------
        submitter : str
            Name of the submitting organization (case-insensitive contains).
        min_classification : str, optional
            Filter to minimum classification level.
        as_set : bool
            If True, return a Python set of gene symbols.

        Returns
        -------
        set[str] or hl.Table
            Gene symbols asserted by the submitter.
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        submitter_lower = submitter.lower()
        if self._keying_mode == "gene_disease_submitter":
            ht = ht.filter(hl.str(ht.submitter).lower().contains(submitter_lower))
        elif self._keying_mode == "gene_disease":
            ht = ht.filter(
                hl.array(ht.submitters).any(
                    lambda s: hl.str(s).lower().contains(submitter_lower)
                )
            )
        else:
            row_fields = set(ht.row)
            if "submitters" in row_fields:
                ht = ht.filter(
                    hl.array(ht.submitters).any(
                        lambda s: hl.str(s).lower().contains(submitter_lower)
                    )
                )

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        return self._return_gene_symbols(ht, as_set=as_set)

    def submitter_summary(self) -> pd.DataFrame:
        """Generate a summary table by submitting organization.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: submitter, n_assertions, n_genes,
            n_diseases, top_classification, last_submission.
        """
        self._ensure_submitter_mode("submitter_summary")
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        summary_ht = ht.group_by(sub=ht.submitter).aggregate(
            n_assertions=hl.agg.count(),
            genes=hl.agg.collect_as_set(ht.gene_symbol),
            diseases=hl.agg.collect_as_set(ht.disease_label),
            classification_counts=hl.agg.counter(ht.classification),
            last_submission=hl.agg.max(ht.submission_date),
        )

        df = summary_ht.to_pandas()
        df = df.rename(columns={"sub": "submitter"})
        df["n_genes"] = df["genes"].apply(len)
        df["n_diseases"] = df["diseases"].apply(len)
        df["top_classification"] = df["classification_counts"].apply(
            lambda counts: max(counts, key=counts.get) if counts else None
        )
        return df.drop(
            columns=["classification_counts", "genes", "diseases"]
        ).sort_values("n_assertions", ascending=False)

    def get_geneset_per_submitter(
        self,
        min_classification: Optional[str] = None,
        min_genes: int = 0,
    ) -> Dict[str, Set[str]]:
        """Get genesets grouped by submitting organization.

        Parameters
        ----------
        min_classification : str, optional
            Filter to minimum classification level before grouping.
        min_genes : int
            Exclude submitters with fewer than this many genes (default: 0).

        Returns
        -------
        dict
            Mapping of submitter name -> set of gene symbols.
        """
        self._ensure_submitter_mode("get_geneset_per_submitter")
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        grouped = ht.group_by(sub=ht.submitter).aggregate(
            genes=hl.agg.collect_as_set(ht.gene_symbol)
        )
        rows = grouped.collect()

        result: Dict[str, Set[str]] = {}
        for row in rows:
            submitter_name = row.sub
            if not submitter_name:
                continue
            genes = set(row.genes)
            if len(genes) < min_genes:
                continue
            result[submitter_name] = genes

        logger.info(
            f"Built {len(result)} submitter-based gene sets"
            + (f" (min_genes={min_genes})" if min_genes > 0 else "")
        )
        return result

    def compare_submitters(self, gene_symbol: str) -> pd.DataFrame:
        """Compare assertions from different submitters for a specific gene.

        Parameters
        ----------
        gene_symbol : str
            Gene symbol to look up.

        Returns
        -------
        pd.DataFrame
            DataFrame with one row per submitter assertion, showing
            submitter, disease, classification, and submission_date.
        """
        self._ensure_submitter_mode("compare_submitters")
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        ht = ht.filter(ht.gene_symbol == gene_symbol)
        rows = ht.collect()

        data = []
        for row in rows:
            data.append(
                {
                    "submitter": row.submitter,
                    "disease": row.disease_label,
                    "classification": row.classification,
                    "submission_date": row.submission_date,
                }
            )

        df = pd.DataFrame(data)
        if not df.empty:
            df = df.sort_values(["disease", "submitter"]).reset_index(drop=True)
        return df

    def consensus_genes(
        self,
        min_submitters: int = 2,
        min_classification: Optional[str] = None,
    ) -> Set[str]:
        """Get genes with assertions from multiple submitters.

        Parameters
        ----------
        min_submitters : int
            Minimum number of distinct submitters required (default: 2).
        min_classification : str, optional
            Filter to minimum classification level before counting.

        Returns
        -------
        set[str]
            Gene symbols meeting the consensus threshold.
        """
        self._ensure_submitter_mode("consensus_genes")
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        grouped = ht.group_by(ht.gene_symbol).aggregate(
            submitters=hl.agg.collect_as_set(ht.submitter),
        )
        # Filter by number of distinct submitters
        filtered = grouped.filter(hl.len(grouped.submitters) >= min_submitters)
        genes = set(filtered.aggregate(hl.agg.collect_as_set(filtered.gene_symbol)))

        logger.info(
            f"Found {len(genes)} consensus genes "
            f"(min_submitters={min_submitters})"
        )
        return genes

    # ------------------------------------------------------------------
    # Overrides for gene_disease_submitter keying mode
    # ------------------------------------------------------------------

    def _apply_min_classification_filter(
        self, ht: hl.Table, min_classification: str
    ) -> hl.Table:
        """Extend base filter to handle gene_disease_submitter mode."""
        if self._keying_mode == "gene_disease_submitter":
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
        return super()._apply_min_classification_filter(ht, min_classification)

    def validate_table(self) -> bool:
        """Extend base validation for gene_disease_submitter keying."""
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("GenCC table not loaded.")

        if self._keying_mode == "gene_disease_submitter":
            row_fields = set(ht.row)
            required = {
                "hgnc_id",
                "mondo_id",
                "submitter",
                "gene_symbol",
                "disease_label",
                "mode_of_inheritance",
                "classification",
                "classification_level",
                "submission_date",
            }
            missing = sorted(required - row_fields)
            if missing:
                raise ValueError(
                    f"GenCC table missing required fields: {', '.join(missing)}"
                )
            return True

        return super().validate_table()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _ensure_submitter_mode(self, method: str) -> None:
        """Ensure table is in a mode with per-row scalar submitter access."""
        self._ensure_table_loaded()
        if self._keying_mode != "gene_disease_submitter":
            raise ValueError(
                f"{method} requires a GenCC table keyed by gene_disease_submitter."
            )
