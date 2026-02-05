"""ClinGen Gene-Disease Validity data streamer."""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import hail as hl
import pandas as pd

from hvantk.core.constants import CLINGEN_CLASSIFICATION_LEVELS
from hvantk.data.data_streamer import HailDataStreamer
from hvantk.data.dataset import get_gene_ann_ht
from hvantk.utils.gene_sets import load_gene_sets_from_dict
from hvantk.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


class ClinGenStreamer(HailDataStreamer):
    """
    Streamer for querying and aggregating ClinGen Gene-Disease Validity data.

    Provides methods for:
    - Filtering genes by classification levels and disease categories
    - Computing summary statistics
    - Aggregating by various dimensions (MOI, GCEP, disease)
    - Integration with other hvantk gene tables
    """

    _MOI_ALIASES = {
        "AD": "Autosomal dominant",
        "AR": "Autosomal recessive",
        "XL": "X-linked",
        "XLD": "X-linked dominant",
        "XLR": "X-linked recessive",
    }

    def __init__(
        self,
        table_path: str,
        chunk_size: int = 10000,
        init_hail: bool = True,
    ):
        """
        Initialize ClinGen streamer.

        Parameters
        ----------
        table_path : str
            Path to the ClinGen Hail Table (.ht directory)
        chunk_size : int
            Chunk size for streaming operations (default: 10000)
        init_hail : bool
            Whether to initialize Hail if not already done (default: True)
        """
        super().__init__("ClinGenStreamer", chunk_size=chunk_size, init_hail=init_hail)
        self.table_path = table_path
        self._table: Optional[hl.Table] = None
        self._keying_mode: Optional[str] = None
        self._filtered_cache: Dict[str, hl.Table] = {}

    def setup(self) -> None:
        super().setup()
        if self._table is not None:
            return
        if not self.table_path:
            raise ValueError("table_path must be a non-empty string")
        self._ensure_table_path_exists()
        self._table = hl.read_table(self.table_path)
        self._keying_mode = self._detect_keying_mode(self._table)
        self.validate_table()

    def stream(self):
        self._ensure_table_loaded()
        if self._table is None:
            raise ValueError("ClinGen table not loaded. Call setup() first.")
        total = self._table.count()
        if total == 0:
            return
        ht = self._table.add_index("row_idx")
        n_chunks = (total + self.chunk_size - 1) // self.chunk_size
        for i in range(n_chunks):
            start = i * self.chunk_size
            end = min(start + self.chunk_size, total)
            chunk = ht.filter((ht.row_idx >= start) & (ht.row_idx < end)).drop(
                "row_idx"
            )
            yield self.process_chunk(chunk)

    def process_chunk(self, chunk: hl.Table) -> hl.Table:
        return chunk

    def get_genes_by_classification(
        self,
        min_classification: str = "Moderate",
        classifications: Optional[List[str]] = None,
        as_set: bool = True,
    ) -> Union[Set[str], hl.Table]:
        """
        Get genes filtered by classification level(s).

        Parameters
        ----------
        min_classification : str
            Minimum classification level (inclusive). Options:
            "Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted"
        classifications : list of str, optional
            Specific classification levels to include. If provided,
            overrides min_classification.
        as_set : bool
            If True, return a Python set of gene symbols.
            If False, return filtered Hail Table.

        Returns
        -------
        set[str] or hl.Table
            Gene symbols matching criteria, or filtered table.
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        if classifications:
            normalized = [self._normalize_classification(c) for c in classifications]
            desired = hl.literal(set(normalized))
            if self._keying_mode == "gene_disease":
                ht = ht.filter(desired.contains(ht.classification))
            else:
                ht = ht.filter(hl.len(ht.classifications.intersection(desired)) > 0)
        else:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        return self._return_gene_symbols(ht, as_set=as_set)

    def get_genes_by_disease(
        self,
        disease_terms: Union[str, List[str]],
        match_mode: str = "contains",
        min_classification: Optional[str] = None,
        as_set: bool = True,
    ) -> Union[Set[str], hl.Table]:
        """
        Get genes associated with disease(s) matching given terms.

        Parameters
        ----------
        disease_terms : str or list of str
            Disease label(s) or search term(s) to match.
        match_mode : str
            How to match disease labels:
            - "exact": Exact match (case-insensitive)
            - "contains": Label contains term (case-insensitive)
            - "regex": Regular expression match
        min_classification : str, optional
            Filter to minimum classification level.
        as_set : bool
            If True, return Python set of gene symbols.

        Returns
        -------
        set[str] or hl.Table
            Genes associated with matching diseases.
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        terms = self._normalize_terms(disease_terms)
        match_mode = match_mode.lower()
        if match_mode not in ("exact", "contains", "regex"):
            raise ValueError("match_mode must be one of: exact, contains, regex")

        if self._keying_mode == "gene_disease":
            match_expr = self._match_disease_expr(ht.disease_label, terms, match_mode)
        else:
            labels = hl.or_else(
                hl.array(ht.disease_labels), hl.empty_array(hl.tstr)
            )
            match_expr = labels.any(
                lambda label: self._match_disease_expr(label, terms, match_mode)
            )

        ht = ht.filter(match_expr)
        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        return self._return_gene_symbols(ht, as_set=as_set)

    def get_genes_by_mondo_id(
        self,
        mondo_ids: Union[str, List[str]],
        min_classification: Optional[str] = None,
        as_set: bool = True,
    ) -> Union[Set[str], hl.Table]:
        """
        Get genes associated with specific MONDO disease ID(s).

        Parameters
        ----------
        mondo_ids : str or list of str
            MONDO ID(s) (with or without "MONDO:" prefix).
        min_classification : str, optional
            Filter to minimum classification level.
        as_set : bool
            If True, return Python set of gene symbols.

        Returns
        -------
        set[str] or hl.Table
            Genes associated with specified MONDO IDs.
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        normalized_ids = {self._normalize_mondo_id(m) for m in self._ensure_list(mondo_ids)}
        mondo_set = hl.literal(normalized_ids)
        if self._keying_mode == "gene_disease":
            ht = ht.filter(mondo_set.contains(ht.mondo_id))
        else:
            ht = ht.filter(hl.len(ht.mondo_ids.intersection(mondo_set)) > 0)

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        return self._return_gene_symbols(ht, as_set=as_set)

    def get_genes_by_moi(
        self,
        modes: Union[str, List[str]],
        min_classification: Optional[str] = None,
        as_set: bool = True,
    ) -> Union[Set[str], hl.Table]:
        """
        Get genes by mode of inheritance.

        Parameters
        ----------
        modes : str or list of str
            Mode(s) of inheritance: "AD" (autosomal dominant),
            "AR" (autosomal recessive), "XL" (X-linked), etc.
        min_classification : str, optional
            Filter to minimum classification level.
        as_set : bool
            If True, return Python set of gene symbols.

        Returns
        -------
        set[str] or hl.Table
            Genes with specified inheritance mode(s).
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        normalized = self._normalize_moi_terms(modes)
        normalized_set = hl.literal(normalized)

        if self._keying_mode == "gene_disease":
            ht = ht.filter(
                normalized_set.contains(hl.str(ht.mode_of_inheritance).lower())
            )
        else:
            modes_array = hl.or_else(
                hl.array(ht.modes_of_inheritance), hl.empty_array(hl.tstr)
            )
            ht = ht.filter(
                modes_array.any(
                    lambda moi: normalized_set.contains(hl.str(moi).lower())
                )
            )

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        return self._return_gene_symbols(ht, as_set=as_set)

    def get_geneset_per_disease(
        self,
        min_classification: Optional[str] = None,
        include_mondo_id: bool = False,
    ) -> Dict[str, Set[str]]:
        """
        Get genesets grouped by individual disease labels.

        Returns a dictionary where each key is a disease label and the value
        is a set of gene symbols associated with that disease.

        Parameters
        ----------
        min_classification : str, optional
            Filter to minimum classification level before grouping.
        include_mondo_id : bool
            If True, returns dict with tuples (disease_label, mondo_id) as keys.

        Returns
        -------
        dict
            Mapping of disease_label -> set of gene symbols.
            If include_mondo_id is True, returns {(disease_label, mondo_id): genes}.

        Examples
        --------
        >>> streamer = ClinGenStreamer("clingen.ht")
        >>> streamer.setup()
        >>> genesets = streamer.get_geneset_per_disease(min_classification="Moderate")
        >>> print(genesets["breast-ovarian cancer, familial 1"])
        {'BRCA1'}
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        # Apply classification filter if specified
        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        if self._keying_mode == "gene_disease":
            # Group by disease and collect genes
            if include_mondo_id:
                grouped = ht.group_by("disease_label", "mondo_id").aggregate(
                    genes=hl.agg.collect_as_set(ht.gene_symbol)
                )
                rows = grouped.collect()
                return {
                    (row.disease_label, row.mondo_id): set(row.genes)
                    for row in rows
                }
            else:
                grouped = ht.group_by("disease_label").aggregate(
                    genes=hl.agg.collect_as_set(ht.gene_symbol)
                )
                rows = grouped.collect()
                return {row.disease_label: set(row.genes) for row in rows}
        else:
            # Gene-keyed table: disease_labels is already a set per gene
            # Need to invert the mapping
            if include_mondo_id:
                row_fields = get_row_fields(ht)
                if "disease_mondo_pairs" not in row_fields:
                    raise ValueError(
                        "include_mondo_id requires disease_mondo_pairs on gene-keyed "
                        "ClinGen tables. Rebuild with create_clingen_gene_disease_tb "
                        "(key_by='gene') from a recent hvantk version."
                    )
                rows = ht.select("gene_symbol", "disease_mondo_pairs").collect()
                disease_to_genes: Dict[Tuple[str, Optional[str]], Set[str]] = {}
                for row in rows:
                    for pair in row.disease_mondo_pairs or []:
                        key = (pair.disease_label, pair.mondo_id)
                        if key not in disease_to_genes:
                            disease_to_genes[key] = set()
                        disease_to_genes[key].add(row.gene_symbol)
                return disease_to_genes

            rows = ht.select("gene_symbol", "disease_labels").collect()
            disease_to_genes = {}
            for row in rows:
                for disease in row.disease_labels or []:
                    if disease not in disease_to_genes:
                        disease_to_genes[disease] = set()
                    disease_to_genes[disease].add(row.gene_symbol)
            return disease_to_genes

    def compute_stats(self) -> Dict[str, Any]:
        """
        Compute comprehensive statistics about the ClinGen dataset.

        Returns
        -------
        dict
            Statistics dictionary containing:
            - total_associations: Total gene-disease pairs
            - unique_genes: Number of unique genes
            - unique_diseases: Number of unique diseases
            - classification_counts: Dict of classification -> count
            - moi_counts: Dict of MOI -> count
            - gcep_counts: Dict of GCEP -> count
            - genes_per_classification: Dict of classification -> gene count
            - diseases_per_classification: Dict of classification -> disease count
            - top_genes_by_diseases: List of (gene, n_diseases) tuples
            - last_update: Most recent classification date
        """
        self._ensure_gene_disease_mode("compute_stats")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        total_associations = ht.count()
        unique_genes = len(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol)))
        unique_diseases = len(ht.aggregate(hl.agg.collect_as_set(ht.disease_label)))
        classification_counts = ht.aggregate(hl.agg.counter(ht.classification))
        moi_counts = ht.aggregate(hl.agg.counter(ht.mode_of_inheritance))
        gcep_counts = ht.aggregate(hl.agg.counter(ht.gene_curation_expert_panel))

        genes_per_classification = ht.aggregate(
            hl.agg.group_by(
                ht.classification, hl.agg.collect_as_set(ht.gene_symbol)
            )
        )
        genes_per_classification = {
            k: len(v) for k, v in genes_per_classification.items()
        }
        diseases_per_classification = ht.aggregate(
            hl.agg.group_by(
                ht.classification, hl.agg.collect_as_set(ht.disease_label)
            )
        )
        diseases_per_classification = {
            k: len(v) for k, v in diseases_per_classification.items()
        }

        top_genes_rows = (
            ht.group_by(ht.gene_symbol)
            .aggregate(n_diseases=hl.agg.count())
            .order_by(hl.desc("n_diseases"))
            .take(10)
        )
        top_genes_by_diseases = [
            (row.gene_symbol, row.n_diseases) for row in top_genes_rows
        ]

        # Get most recent classification date using Hail aggregator (ignores missing)
        last_update = ht.aggregate(hl.agg.max(ht.classification_date))

        return {
            "total_associations": total_associations,
            "unique_genes": unique_genes,
            "unique_diseases": unique_diseases,
            "classification_counts": classification_counts,
            "moi_counts": moi_counts,
            "gcep_counts": gcep_counts,
            "genes_per_classification": genes_per_classification,
            "diseases_per_classification": diseases_per_classification,
            "top_genes_by_diseases": top_genes_by_diseases,
            "last_update": last_update,
        }

    def classification_summary(self) -> pd.DataFrame:
        """
        Generate a summary table of classifications.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - classification: Classification level
            - n_associations: Number of gene-disease pairs
            - n_genes: Unique genes at this level
            - n_diseases: Unique diseases at this level
            - pct_associations: Percentage of total associations
        """
        self._ensure_gene_disease_mode("classification_summary")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        summary_ht = ht.group_by(
            classification=ht.classification,
            classification_level=ht.classification_level,
        ).aggregate(
            n_associations=hl.agg.count(),
            genes=hl.agg.collect_as_set(ht.gene_symbol),
            diseases=hl.agg.collect_as_set(ht.disease_label),
        )

        total = ht.count()
        df = summary_ht.to_pandas()
        df["n_genes"] = df["genes"].apply(len)
        df["n_diseases"] = df["diseases"].apply(len)
        df["pct_associations"] = (df["n_associations"] / total) * 100
        df = df.sort_values("classification_level").reset_index(drop=True)
        return df.drop(columns=["classification_level", "genes", "diseases"])

    def gcep_summary(self) -> pd.DataFrame:
        """
        Generate a summary table by Gene Curation Expert Panel.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - gcep: Expert panel name
            - n_associations: Number of curations
            - n_genes: Unique genes curated
            - n_diseases: Unique diseases curated
            - top_classification: Most common classification
            - last_curation: Most recent curation date
        """
        self._ensure_gene_disease_mode("gcep_summary")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        summary_ht = ht.group_by(
            gcep=ht.gene_curation_expert_panel
        ).aggregate(
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
        return df.drop(columns=["classification_counts", "genes", "diseases"]).sort_values(
            "n_associations", ascending=False
        )

    def aggregate_by_disease_category(
        self,
        categories: Dict[str, List[str]],
        min_classification: Optional[str] = None,
    ) -> Dict[str, Set[str]]:
        """
        Aggregate genes into user-defined disease categories.

        Parameters
        ----------
        categories : dict
            Mapping of category name to list of search terms.
        min_classification : str, optional
            Filter to minimum classification level.

        Returns
        -------
        dict
            Mapping of category name to set of gene symbols.
        """
        gene_sets: Dict[str, Set[str]] = {}
        for name, terms in categories.items():
            genes = self.get_genes_by_disease(
                terms,
                match_mode="contains",
                min_classification=min_classification,
                as_set=True,
            )
            gene_sets[name] = genes
        return gene_sets

    def categorize_by_ontology(
        self,
        mondo_obo_path: str,
        min_classification: Optional[str] = None,
        categories: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Dict[str, Set[str]]]:
        """
        Categorize diseases using MONDO ontology hierarchy.

        This method uses the MONDO disease ontology to properly categorize
        diseases based on their ontological relationships (is_a hierarchy),
        rather than keyword matching.

        Parameters
        ----------
        mondo_obo_path : str
            Path to the MONDO OBO file. Can be downloaded using:
            `hvantk.utils.mondo_parser.download_mondo_obo()`
        min_classification : str, optional
            Filter to minimum classification level.
        categories : dict, optional
            Custom category mapping {MONDO_ID: category_name}.
            If None, uses default MONDO_DISEASE_CATEGORIES.

        Returns
        -------
        dict
            Mapping of category_name -> {
                "genes": set of gene symbols,
                "diseases": set of disease labels,
                "mondo_ids": set of MONDO IDs
            }

        Examples
        --------
        >>> streamer = ClinGenStreamer("clingen.ht")
        >>> streamer.setup()
        >>> results = streamer.categorize_by_ontology("mondo.obo")
        >>> print(results["cardiovascular disease"]["genes"])
        {'TTN', 'MYH7', 'MYBPC3', ...}
        """
        from hvantk.utils.mondo_parser import MondoOntology, MONDO_DISEASE_CATEGORIES

        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        # Apply classification filter if specified
        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        # Load the MONDO ontology
        logger.info(f"Loading MONDO ontology from {mondo_obo_path}")
        mondo = MondoOntology(mondo_obo_path)

        # Use default categories if not provided
        if categories is None:
            categories = MONDO_DISEASE_CATEGORIES

        # Collect all gene-disease-mondo associations
        if self._keying_mode == "gene_disease":
            # Key fields can't be selected directly, so collect the full rows
            # and extract the needed fields
            rows = ht.collect()
        else:
            # For gene-keyed tables, we need to expand
            row_fields = get_row_fields(ht)
            if "disease_mondo_pairs" not in row_fields:
                raise ValueError(
                    "Gene-keyed ClinGen tables must include disease_mondo_pairs to "
                    "categorize by ontology. Rebuild with "
                    "create_clingen_gene_disease_tb(key_by='gene') from a recent "
                    "hvantk version."
                )
            rows = []
            gene_rows = ht.select("gene_symbol", "disease_mondo_pairs").collect()
            for row in gene_rows:
                for pair in row.disease_mondo_pairs or []:
                    rows.append(
                        {
                            "gene_symbol": row.gene_symbol,
                            "disease_label": pair.disease_label,
                            "mondo_id": pair.mondo_id,
                        }
                    )

        # Categorize each disease and aggregate
        results: Dict[str, Dict[str, Set[str]]] = {}
        for cat_name in categories.values():
            results[cat_name] = {"genes": set(), "diseases": set(), "mondo_ids": set()}

        uncategorized = {"genes": set(), "diseases": set(), "mondo_ids": set()}

        for row in rows:
            gene = row.gene_symbol if hasattr(row, 'gene_symbol') else row["gene_symbol"]
            disease = row.disease_label if hasattr(row, 'disease_label') else row["disease_label"]
            mondo_id = row.mondo_id if hasattr(row, 'mondo_id') else row["mondo_id"]

            if not mondo_id:
                uncategorized["genes"].add(gene)
                uncategorized["diseases"].add(disease)
                continue

            # Normalize MONDO ID
            if not mondo_id.startswith("MONDO:"):
                mondo_id = f"MONDO:{mondo_id}"

            # Get categories for this disease
            matched_cats = mondo.categorize(mondo_id, categories)

            if matched_cats:
                for cat_id, cat_name in matched_cats:
                    results[cat_name]["genes"].add(gene)
                    results[cat_name]["diseases"].add(disease)
                    results[cat_name]["mondo_ids"].add(mondo_id)
            else:
                uncategorized["genes"].add(gene)
                uncategorized["diseases"].add(disease)
                uncategorized["mondo_ids"].add(mondo_id)

        # Add uncategorized if non-empty
        if uncategorized["genes"]:
            results["uncategorized"] = uncategorized

        logger.info(f"Categorized diseases into {len(results)} categories")
        return results

    def categorize_by_ontology_summary(
        self,
        mondo_obo_path: str,
        min_classification: Optional[str] = None,
        categories: Optional[Dict[str, str]] = None,
    ) -> pd.DataFrame:
        """
        Get a summary DataFrame of ontology-based disease categorization.

        Parameters
        ----------
        mondo_obo_path : str
            Path to the MONDO OBO file.
        min_classification : str, optional
            Filter to minimum classification level.
        categories : dict, optional
            Custom category mapping {MONDO_ID: category_name}.

        Returns
        -------
        pd.DataFrame
            Summary with columns: category, n_genes, n_diseases, sample_genes
        """
        results = self.categorize_by_ontology(
            mondo_obo_path, min_classification, categories
        )

        summary_data = []
        for category, data in sorted(results.items(), key=lambda x: -len(x[1]["genes"])):
            genes = data["genes"]
            summary_data.append({
                "category": category,
                "n_genes": len(genes),
                "n_diseases": len(data["diseases"]),
                "sample_genes": ", ".join(sorted(genes)[:10]) + ("..." if len(genes) > 10 else ""),
            })

        return pd.DataFrame(summary_data)

    def get_genes_by_ontology_category(
        self,
        mondo_obo_path: str,
        category_id: str,
        min_classification: Optional[str] = None,
        as_set: bool = True,
    ) -> Union[Set[str], List[Tuple[str, str, str]]]:
        """
        Get genes belonging to a specific ontology category.

        Parameters
        ----------
        mondo_obo_path : str
            Path to the MONDO OBO file.
        category_id : str
            MONDO ID of the category (e.g., "MONDO:0004995" for cardiovascular).
        min_classification : str, optional
            Filter to minimum classification level.
        as_set : bool
            If True, return set of gene symbols. If False, return list of
            (gene, disease, mondo_id) tuples.

        Returns
        -------
        set or list
            Genes in the specified category.
        """
        from hvantk.utils.mondo_parser import MondoOntology

        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        if min_classification:
            min_classification = self._normalize_classification(min_classification)
            ht = self._apply_min_classification_filter(ht, min_classification)

        # Load ontology
        mondo = MondoOntology(mondo_obo_path)

        # Get all descendants of the category (all diseases under this category)
        category_diseases = mondo.get_descendants(category_id, include_self=True)

        # Collect genes from ClinGen data (can't select key fields directly)
        genes = set()
        full_results = []

        if self._keying_mode == "gene_disease":
            rows = ht.collect()
            for row in rows:
                mondo_id = row.mondo_id
                if not mondo_id:
                    continue
                if not mondo_id.startswith("MONDO:"):
                    mondo_id = f"MONDO:{mondo_id}"

                if mondo_id in category_diseases:
                    genes.add(row.gene_symbol)
                    full_results.append((row.gene_symbol, row.disease_label, mondo_id))
        else:
            row_fields = get_row_fields(ht)
            if "disease_mondo_pairs" not in row_fields:
                raise ValueError(
                    "Gene-keyed ClinGen tables must include disease_mondo_pairs to "
                    "query ontology categories. Rebuild with "
                    "create_clingen_gene_disease_tb(key_by='gene') from a recent "
                    "hvantk version."
                )
            rows = ht.select("gene_symbol", "disease_mondo_pairs").collect()
            for row in rows:
                gene_symbol = row.gene_symbol
                for pair in row.disease_mondo_pairs or []:
                    mondo_id = pair.mondo_id
                    if not mondo_id:
                        continue
                    if not mondo_id.startswith("MONDO:"):
                        mondo_id = f"MONDO:{mondo_id}"

                    if mondo_id in category_diseases:
                        genes.add(gene_symbol)
                        full_results.append(
                            (gene_symbol, pair.disease_label, mondo_id)
                        )

        if as_set:
            return genes
        return full_results

    def group_by_gcep(
        self,
        min_classification: Optional[str] = None,
    ) -> Dict[str, hl.Table]:
        """
        Group data by Gene Curation Expert Panel.

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
        return {gcep: ht.filter(ht.gene_curation_expert_panel == gcep) for gcep in gceps}

    def pivot_genes_by_classification(self) -> pd.DataFrame:
        """
        Create a pivot table of genes x classifications.

        Returns
        -------
        pd.DataFrame
            Pivot table where:
            - Rows: gene symbols
            - Columns: classification levels
            - Values: count of diseases at each level
        """
        self._ensure_gene_disease_mode("pivot_genes_by_classification")
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        grouped = ht.group_by(ht.gene_symbol).aggregate(
            classification_counts=hl.agg.counter(ht.classification)
        )
        df = grouped.to_pandas()
        for level in CLINGEN_CLASSIFICATION_LEVELS:
            df[level] = df["classification_counts"].apply(
                lambda counts, lvl=level: counts.get(lvl, 0) if counts else 0
            )
        df = df.drop(columns=["classification_counts"]).set_index("gene_symbol")
        return df

    def to_gene_set(
        self,
        min_classification: str = "Moderate",
        id_type: str = "symbol",
    ) -> Set[str]:
        """
        Export as gene set for integration with other hvantk components.

        Parameters
        ----------
        min_classification : str
            Minimum classification level to include.
        id_type : str
            Type of gene identifier: "symbol", "hgnc_id", or "both".

        Returns
        -------
        set[str]
            Gene identifiers meeting criteria.
        """
        id_type = id_type.lower()
        if id_type not in ("symbol", "hgnc_id", "both"):
            raise ValueError("id_type must be one of: symbol, hgnc_id, both")

        min_classification = self._normalize_classification(min_classification)
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        ht = self._apply_min_classification_filter(ht, min_classification)

        if id_type == "symbol":
            return self._collect_set(ht, ht.gene_symbol)
        if id_type == "hgnc_id":
            return self._collect_set(ht, ht.hgnc_id)

        symbols = self._collect_set(ht, ht.gene_symbol)
        hgnc_ids = self._collect_set(ht, ht.hgnc_id)
        return symbols | hgnc_ids

    def annotate_gene_table(
        self,
        gene_table: hl.Table,
        gene_id_field: str = "gene_id",
        annotation_fields: Optional[List[str]] = None,
    ) -> hl.Table:
        """
        Annotate an existing gene table with ClinGen data.

        Parameters
        ----------
        gene_table : hl.Table
            Input gene table to annotate.
        gene_id_field : str
            Name of the gene identifier field in input table.
            Supports: "gene_id" (Ensembl), "gene_symbol", "hgnc_id".
        annotation_fields : list of str, optional
            Fields to add from ClinGen. If None, adds:
            - clingen_classifications
            - clingen_diseases
            - clingen_max_classification
            - clingen_n_diseases

        Returns
        -------
        hl.Table
            Annotated gene table with ClinGen fields.
        """
        self._ensure_table_loaded()
        if gene_id_field not in gene_table.row:
            raise ValueError(
                f"Field {gene_id_field} not found in gene_table row fields"
            )

        clingen_gene_ht = self._get_gene_level_table()
        gene_table_keyed, join_expr, clingen_key = self._prepare_gene_table_key(
            gene_table, gene_id_field
        )
        clingen_gene_ht = clingen_gene_ht.key_by(clingen_key)

        if annotation_fields is None:
            mapping = {
                "clingen_classifications": "classifications",
                "clingen_diseases": "disease_labels",
                "clingen_max_classification": "max_classification_label",
                "clingen_n_diseases": "n_diseases",
            }
        else:
            mapping = {f"clingen_{field}": field for field in annotation_fields}

        annotations = {}
        for out_field, in_field in mapping.items():
            annotations[out_field] = clingen_gene_ht[join_expr][in_field]
        return gene_table_keyed.annotate(**annotations)

    def export_for_enrichex(
        self,
        output_path: str,
        categories: Optional[Dict[str, List[str]]] = None,
        min_classification: str = "Moderate",
    ) -> None:
        """
        Export gene sets in format compatible with EnrichEx.

        Parameters
        ----------
        output_path : str
            Path to write gene set file (GMT or JSON format).
        categories : dict, optional
            Disease categories to export. If None, exports
            one set per classification level.
        min_classification : str
            Minimum classification level to include.
        """
        min_classification = self._normalize_classification(min_classification)
        if categories:
            gene_sets = self.aggregate_by_disease_category(
                categories, min_classification=min_classification
            )
        else:
            min_level = self._classification_level(min_classification)
            gene_sets = {}
            for level in CLINGEN_CLASSIFICATION_LEVELS:
                if self._classification_level(level) <= min_level:
                    genes = self.get_genes_by_classification(
                        classifications=[level], as_set=True
                    )
                    gene_sets[level] = genes

        if output_path.lower().endswith(".gmt"):
            self._write_gmt(output_path, gene_sets)
            return

        collection = load_gene_sets_from_dict(
            {k: sorted(v) for k, v in gene_sets.items()},
            source="clingen_gene_disease_validity",
        )
        collection.save(output_path)

    def describe(self) -> str:
        """
        Return a human-readable description of the dataset.

        Returns
        -------
        str
            Multi-line string with dataset summary.
        """
        stats = self.compute_stats()
        total = stats["total_associations"]
        lines = [
            "ClinGen Gene-Disease Validity Dataset",
            "=" * 36,
            f"Total associations: {total:,}",
            f"Unique genes: {stats['unique_genes']:,}",
            f"Unique diseases: {stats['unique_diseases']:,}",
            "",
            "Classification distribution:",
        ]
        for level in CLINGEN_CLASSIFICATION_LEVELS:
            count = stats["classification_counts"].get(level, 0)
            pct = (count / total * 100) if total else 0.0
            lines.append(f"  {level}: {count:,} ({pct:.1f}%)")
        return "\n".join(lines)

    def validate_table(self) -> bool:
        """
        Validate that the loaded table has expected ClinGen schema.

        Returns
        -------
        bool
            True if schema is valid.

        Raises
        ------
        ValueError
            If required fields are missing or have wrong types.
        """
        self._ensure_table_loaded()
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")

        row_fields = set(ht.row)
        if self._keying_mode == "gene_disease":
            required = {
                "hgnc_id",
                "mondo_id",
                "gene_symbol",
                "disease_label",
                "mode_of_inheritance",
                "classification",
                "classification_level",
                "classification_date",
                "gene_curation_expert_panel",
            }
        elif self._keying_mode == "gene":
            required = {
                "hgnc_id",
                "gene_symbol",
                "disease_labels",
                "mondo_ids",
                "classifications",
                "modes_of_inheritance",
                "max_classification_level",
                "max_classification_label",
                "n_diseases",
            }
        else:
            raise ValueError(
                "ClinGen table has unsupported keying. Expected key_by gene_disease "
                "or gene."
            )

        missing = sorted(required - row_fields)
        if missing:
            raise ValueError(
                f"ClinGen table missing required fields: {', '.join(missing)}"
            )

        return True

    def _ensure_table_loaded(self) -> None:
        if self._table is None:
            self.setup()

    def _ensure_table_path_exists(self) -> None:
        if hl.hadoop_exists(self.table_path):
            return
        if os.path.exists(self.table_path):
            return
        raise FileNotFoundError(f"ClinGen table path not found: {self.table_path}")

    def _detect_keying_mode(self, ht: hl.Table) -> str:
        key_fields = set(ht.key)
        if {"hgnc_id", "mondo_id"} <= key_fields:
            return "gene_disease"
        if "hgnc_id" in key_fields and "mondo_id" not in key_fields:
            return "gene"
        return "unknown"

    def _ensure_gene_disease_mode(self, method: str) -> None:
        self._ensure_table_loaded()
        if self._keying_mode != "gene_disease":
            raise ValueError(
                f"{method} requires a ClinGen table keyed by gene-disease."
            )

    def _classification_level(self, classification: str) -> int:
        return CLINGEN_CLASSIFICATION_LEVELS.index(classification)

    def _normalize_classification(self, classification: str) -> str:
        if not classification:
            raise ValueError("classification must be a non-empty string")
        normalized = classification.strip()
        lookup = {c.lower(): c for c in CLINGEN_CLASSIFICATION_LEVELS}
        key = normalized.lower()
        if key not in lookup:
            raise ValueError(
                f"Invalid classification '{classification}'. "
                f"Expected one of {CLINGEN_CLASSIFICATION_LEVELS}."
            )
        return lookup[key]

    def _apply_min_classification_filter(
        self, ht: hl.Table, min_classification: str
    ) -> hl.Table:
        min_level = self._classification_level(min_classification)
        cache_key = None
        if ht is self._table:
            cache_key = f"min:{self._keying_mode}:{min_classification}"
            cached = self._filtered_cache.get(cache_key)
            if cached is not None:
                return cached
        if self._keying_mode == "gene_disease":
            filtered = ht.filter(ht.classification_level <= min_level)
        else:
            filtered = ht.filter(ht.max_classification_level <= min_level)
        if cache_key:
            self._filtered_cache[cache_key] = filtered
        return filtered

    def _return_gene_symbols(self, ht: hl.Table, as_set: bool) -> Union[Set[str], hl.Table]:
        if not as_set:
            return ht
        return self._collect_set(ht, ht.gene_symbol)

    def _collect_set(self, ht: hl.Table, expr: hl.expr.Expression) -> Set[str]:
        return set(ht.aggregate(hl.agg.collect_as_set(expr)))

    @staticmethod
    def _ensure_list(values: Union[str, List[str]]) -> List[str]:
        if isinstance(values, str):
            return [values]
        return list(values)

    @staticmethod
    def _normalize_mondo_id(mondo_id: str) -> str:
        mondo_id = mondo_id.strip()
        if mondo_id.startswith("MONDO:"):
            return mondo_id.replace("MONDO:", "")
        return mondo_id

    @staticmethod
    def _normalize_terms(terms: Union[str, List[str]]) -> List[str]:
        if isinstance(terms, str):
            return [terms.strip()] if terms.strip() else []
        return [t.strip() for t in terms if t and t.strip()]

    def _match_disease_expr(
        self, label_expr: hl.expr.Expression, terms: List[str], match_mode: str
    ) -> hl.expr.BooleanExpression:
        label = hl.str(label_expr)
        if match_mode == "exact":
            terms_set = hl.literal({t.lower() for t in terms})
            return terms_set.contains(label.lower())
        if match_mode == "contains":
            terms_array = hl.literal([t.lower() for t in terms])
            return terms_array.any(lambda term: label.lower().contains(term))
        terms_array = hl.literal(list(terms))
        return terms_array.any(lambda term: hl.regex_match(term, label))

    def _normalize_moi_terms(self, modes: Union[str, List[str]]) -> Set[str]:
        modes_list = self._ensure_list(modes)
        normalized: Set[str] = set()
        for mode in modes_list:
            raw = mode.strip()
            if not raw:
                continue
            normalized.add(raw.lower())
            alias = self._MOI_ALIASES.get(raw.upper())
            if alias:
                normalized.add(alias.lower())
        return normalized

    def _get_gene_level_table(self) -> hl.Table:
        ht = self._table
        if ht is None:
            raise ValueError("ClinGen table not loaded.")
        if self._keying_mode == "gene":
            return ht

        gene_ht = (
            ht.group_by("hgnc_id", "gene_symbol")
            .aggregate(
                disease_labels=hl.agg.collect_as_set(ht.disease_label),
                disease_mondo_pairs=hl.agg.collect_as_set(
                    hl.struct(disease_label=ht.disease_label, mondo_id=ht.mondo_id)
                ),
                mondo_ids=hl.agg.collect_as_set(ht.mondo_id),
                classifications=hl.agg.collect_as_set(ht.classification),
                modes_of_inheritance=hl.agg.collect_as_set(ht.mode_of_inheritance),
                max_classification_level=hl.agg.min(ht.classification_level),
                n_diseases=hl.agg.count(),
            )
            .key_by("hgnc_id")
        )
        classification_labels = hl.literal(CLINGEN_CLASSIFICATION_LEVELS)
        safe_index = hl.min(
            gene_ht.max_classification_level, hl.len(classification_labels) - 1
        )
        gene_ht = gene_ht.annotate(
            max_classification_label=classification_labels[safe_index]
        )
        return gene_ht

    def _prepare_gene_table_key(
        self, gene_table: hl.Table, gene_id_field: str
    ) -> tuple[hl.Table, hl.expr.Expression, str]:
        if gene_id_field == "gene_symbol":
            return gene_table, gene_table[gene_id_field], "gene_symbol"
        if gene_id_field == "hgnc_id":
            return gene_table, gene_table[gene_id_field], "hgnc_id"
        if gene_id_field != "gene_id":
            raise ValueError(
                "gene_id_field must be one of: gene_id, gene_symbol, hgnc_id"
            )

        from hvantk.data import dataset as dataset_module

        if not dataset_module.source_dir:
            raise ValueError(
                "hvantk.data.dataset.source_dir is not set; cannot map gene_id to "
                "gene_symbol. Set source_dir or use gene_symbol/hgnc_id."
            )
        gene_ann_ht = get_gene_ann_ht()
        if "gene_name" not in gene_ann_ht.row:
            raise ValueError(
                "Ensembl gene annotation table missing gene_name field for mapping."
            )
        mapped = gene_table.annotate(
            clingen_gene_symbol=gene_ann_ht[gene_table[gene_id_field]].gene_name
        )
        return mapped, mapped.clingen_gene_symbol, "gene_symbol"

    @staticmethod
    def _write_gmt(output_path: str, gene_sets: Dict[str, Set[str]]) -> None:
        with open(output_path, "w", encoding="utf-8") as handle:
            for name, genes in gene_sets.items():
                if not genes:
                    continue
                line = [name, "ClinGen"] + sorted(genes)
                handle.write("\t".join(line) + "\n")
