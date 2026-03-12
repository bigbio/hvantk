"""
BurdenPipeline: orchestrated multi-tissue × multi-variant-class burden testing.

This module provides a ``BurdenConfig`` dataclass and ``BurdenPipeline`` class
for running burden analysis across multiple gene set collections and variant
classes in a single coordinated run.  The pipeline reuses a single cohort
MatrixTable across all runs and generates per-run TSVs plus a combined summary.

Example::

    from hvantk.enrichex.pipeline import BurdenConfig, BurdenPipeline
    from hvantk.enrichex.burden import build_variant_classes_from_presets

    config = BurdenConfig(
        cohort_mt_path="cohort.mt",
        phenotype_ht_path="phenotypes.ht",
        gene_set_collections={
            "heart": "gene_sets/heart_cell_types.json",
            "brain": "gene_sets/brain_cell_types.json",
        },
        variant_classes=build_variant_classes_from_presets(
            ["lof", "missense_constrained", "synonymous"]
        ),
        covariate_fields=["PC1", "PC2", "PC3", "sex"],
        competitive=True,
        output_dir="results/chd_celltype_burden",
    )
    pipeline = BurdenPipeline(config)
    pipeline.show_plan()
    combined_df = pipeline.run()
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

try:
    import hail as hl
except ModuleNotFoundError:  # pragma: no cover
    hl = None  # type: ignore

from hvantk.enrichex.constants import (
    CORRECTION_METHODS,
    GENOTYPE_AGGREGATION_METHODS,
    PHENOTYPE_TYPES,
    _DEPRECATED_AGGREGATION_ALIASES,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class BurdenConfig:
    """Configuration for a multi-tissue, multi-variant-class burden pipeline.

    Attributes
    ----------
    cohort_mt_path : str
        Path to the cohort MatrixTable.
    phenotype_ht_path : str
        Path to phenotype Hail Table (``.ht``), TSV file, or empty string.
        When empty, phenotype and covariates are extracted from the cohort
        MT column fields (see ``phenotype_field`` and ``covariate_fields``).
    phenotype_field : str
        Column field containing the phenotype.  Supports dot notation for
        nested struct fields (e.g., ``"phe.is_case"``).  When using an
        external phenotype HT, this is the field name in that table.
        When using MT column fields (``phenotype_ht_path=""``), this is the
        path in the MT column schema.
    phenotype_type : str
        ``"binary"`` or ``"continuous"``.
    covariate_fields : List[str]
        Column fields to use as covariates (e.g., PCs, sex).  Supports
        dot notation for nested struct fields.
    gene_set_collections : Dict[str, str]
        Mapping of collection name to JSON path
        (as produced by ``GeneSetCollection.save()``).
    variant_classes : Dict[str, VariantFilter]
        Mapping of variant class name to filter criteria.  When empty,
        a single unstratified run is performed.
    gene_field : str
        Row field for gene symbol in the cohort MatrixTable.
    genotype_aggregation : str
        Genotype aggregation method.
    normalize_by_length : bool
        Normalize per-gene burden by gene length.
    gene_lengths_path : Optional[str]
        TSV with columns ``gene`` and ``length_bp``.
    min_carriers : int
        Minimum carriers per gene set for regression.
    correction_method : str
        Multiple testing correction method.
    alpha : float
        Significance threshold.
    competitive : bool
        Run permutation-based competitive test.
    n_permutations : int
        Number of permutations for competitive test.
    permutation_seed : Optional[int]
        Random seed for permutation reproducibility.
    sample_id_field : str
        Column name for sample IDs in TSV phenotype files.
    output_dir : str
        Output directory.
    generate_report : bool
        Generate HTML report after pipeline completes.
    """

    # -- Required inputs --
    cohort_mt_path: str = ""
    phenotype_ht_path: str = ""  # empty → extract from MT column fields

    # -- Phenotype --
    phenotype_field: str = "is_case"
    phenotype_type: str = "binary"
    covariate_fields: List[str] = field(default_factory=list)
    sample_id_field: str = "sample_id"

    # -- Gene sets --
    gene_set_collections: Dict[str, str] = field(default_factory=dict)

    # -- Variant classes (empty → single unstratified run) --
    variant_classes: Dict[str, Any] = field(default_factory=dict)

    # -- Analysis settings --
    gene_field: str = "SYMBOL"
    genotype_aggregation: str = "hets"
    normalize_by_length: bool = False
    gene_lengths_path: Optional[str] = None
    min_carriers: int = 5

    # -- Correction --
    correction_method: str = "benjamini-hochberg"
    alpha: float = 0.05

    # -- Competitive testing --
    competitive: bool = False
    n_permutations: int = 10000
    permutation_seed: Optional[int] = None

    # -- Output --
    output_dir: str = "."
    generate_report: bool = True

    def validate(self) -> List[str]:
        """Validate configuration and return a list of error messages.

        Returns
        -------
        List[str]
            Error messages.  Empty if the configuration is valid.
        """
        errors: List[str] = []

        if not self.cohort_mt_path:
            errors.append("cohort_mt_path is required.")
        elif not Path(self.cohort_mt_path).exists():
            errors.append(f"Cohort MT not found: {self.cohort_mt_path}")

        # phenotype_ht_path is optional — when empty, phenotype and covariates
        # are extracted from MT column fields.
        if self.phenotype_ht_path:
            if not Path(self.phenotype_ht_path).exists():
                errors.append(f"Phenotype file not found: {self.phenotype_ht_path}")

        if not self.gene_set_collections:
            errors.append("At least one gene_set_collection is required.")
        else:
            for name, path in self.gene_set_collections.items():
                if not Path(path).exists():
                    errors.append(f"Gene set collection '{name}' not found: {path}")

        if self.phenotype_type not in PHENOTYPE_TYPES:
            errors.append(
                f"Invalid phenotype_type '{self.phenotype_type}'. "
                f"Must be one of: {PHENOTYPE_TYPES}"
            )

        agg = self.genotype_aggregation
        valid_agg = GENOTYPE_AGGREGATION_METHODS + list(_DEPRECATED_AGGREGATION_ALIASES)
        if agg not in valid_agg:
            errors.append(
                f"Invalid genotype_aggregation '{agg}'. "
                f"Must be one of: {GENOTYPE_AGGREGATION_METHODS}"
            )

        if self.correction_method not in CORRECTION_METHODS:
            errors.append(
                f"Invalid correction_method '{self.correction_method}'. "
                f"Must be one of: {CORRECTION_METHODS}"
            )

        if self.gene_lengths_path and not Path(self.gene_lengths_path).exists():
            errors.append(f"Gene lengths file not found: {self.gene_lengths_path}")

        if self.min_carriers < 0:
            errors.append("min_carriers must be >= 0.")

        if not 0 < self.alpha <= 1:
            errors.append("alpha must be in (0, 1].")

        if self.n_permutations < 1:
            errors.append("n_permutations must be >= 1.")

        return errors


# ---------------------------------------------------------------------------
# Per-run result container
# ---------------------------------------------------------------------------


@dataclass
class BurdenRunResult:
    """Result from a single (variant_class, gene_set_collection) run.

    Attributes
    ----------
    variant_class : str
        Variant class name (or ``"all"`` if unstratified).
    collection_name : str
        Gene set collection name.
    results_df : pd.DataFrame
        Corrected burden results.
    permutation_df : Optional[pd.DataFrame]
        Permutation test results (when ``competitive=True``).
    n_gene_sets_tested : int
        Number of gene sets tested.
    n_significant : int
        Number of significant gene sets.
    """

    variant_class: str
    collection_name: str
    results_df: pd.DataFrame
    permutation_df: Optional[pd.DataFrame] = None
    n_gene_sets_tested: int = 0
    n_significant: int = 0


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


class BurdenPipeline:
    """Orchestrated multi-tissue, multi-variant-class burden testing.

    The pipeline loads the cohort MatrixTable once and iterates over the
    cross-product of ``(variant_class × gene_set_collection)``.  For each
    combination it calls the existing ``run_burden_analysis()`` function,
    applies multiple testing correction, and optionally runs a competitive
    permutation test.

    Parameters
    ----------
    config : BurdenConfig
        Pipeline configuration.  Validated on construction.

    Raises
    ------
    ValueError
        If configuration validation fails.
    """

    def __init__(self, config: BurdenConfig) -> None:
        errors = config.validate()
        if errors:
            raise ValueError(
                "BurdenConfig validation failed:\n"
                + "\n".join(f"  - {e}" for e in errors)
            )
        self.config = config

        # Set up output paths (no Hail needed)
        self._output_dir = Path(config.output_dir)
        self._per_run_dir = self._output_dir / "per_run"
        self._plots_dir = self._output_dir / "plots"
        self._output_dir.mkdir(parents=True, exist_ok=True)
        self._per_run_dir.mkdir(exist_ok=True)
        self._plots_dir.mkdir(exist_ok=True)

        # Internal state (populated during run)
        self._mt: Any = None
        self._phenotype_ht: Any = None
        self._gene_lengths: Optional[Dict[str, float]] = None
        self._collections: Dict[str, Any] = {}  # GeneSetCollections
        self._run_results: List[BurdenRunResult] = []

    # ------------------------------------------------------------------
    # Dry-run
    # ------------------------------------------------------------------

    def show_plan(self) -> str:
        """Display the execution plan without initialising Hail.

        Returns
        -------
        str
            The formatted plan text (also printed to stdout).
        """
        cfg = self.config

        # Determine runs
        vc_names = list(cfg.variant_classes.keys()) if cfg.variant_classes else ["all"]
        coll_names = list(cfg.gene_set_collections.keys())
        n_runs = len(vc_names) * len(coll_names)

        lines = [
            "",
            "=" * 70,
            "BURDEN PIPELINE EXECUTION PLAN",
            "=" * 70,
            "",
            "Inputs:",
            f"  Cohort MT:        {cfg.cohort_mt_path}",
            f"  Phenotypes:       {cfg.phenotype_ht_path or '(from MT column fields)'}",
            f"    Field:          {cfg.phenotype_field} ({cfg.phenotype_type})",
            "",
            "Gene Set Collections:",
        ]
        for name, path in cfg.gene_set_collections.items():
            lines.append(f"  {name}: {path}")

        lines.append("")
        lines.append("Variant Classes:")
        if cfg.variant_classes:
            for vc_name, vf in cfg.variant_classes.items():
                csq = getattr(vf, "consequences", None) or "any"
                lines.append(f"  {vc_name}: consequences={csq}")
        else:
            lines.append("  (none — single unstratified run)")

        lines.append("")
        lines.append("Analysis Settings:")
        lines.append(f"  Gene field:              {cfg.gene_field}")
        lines.append(f"  Genotype aggregation:    {cfg.genotype_aggregation}")
        lines.append(f"  Min carriers:            {cfg.min_carriers}")
        if cfg.covariate_fields:
            lines.append(
                f"  Covariates:              {', '.join(cfg.covariate_fields)}"
            )
        if cfg.normalize_by_length:
            src = cfg.gene_lengths_path or "variant site count proxy"
            lines.append(f"  Length normalisation:     {src}")
        lines.append(f"  Correction:              {cfg.correction_method}")
        lines.append(f"  Alpha:                   {cfg.alpha}")
        if cfg.competitive:
            lines.append(
                f"  Competitive testing:     {cfg.n_permutations} permutations"
            )

        lines.append("")
        lines.append(
            f"Total runs: {n_runs}  ({len(vc_names)} classes x {len(coll_names)} collections)"
        )
        lines.append(f"Output dir: {cfg.output_dir}")
        lines.append("=" * 70)
        lines.append("")

        plan = "\n".join(lines)
        print(plan)
        return plan

    # ------------------------------------------------------------------
    # Full execution
    # ------------------------------------------------------------------

    def run(self) -> pd.DataFrame:
        """Execute the pipeline.

        Returns
        -------
        pd.DataFrame
            Combined results across all (variant_class, collection) runs.
        """
        t0 = time.time()
        logger.info("=" * 70)
        logger.info("BURDEN PIPELINE — START")
        logger.info("=" * 70)

        self._initialize_hail()

        t_load = time.time()
        self._load_cohort()
        self._load_phenotypes()
        self._load_gene_sets()
        self._load_gene_lengths()
        logger.info("Data loading completed in %.1fs", time.time() - t_load)

        # Log summary of loaded data
        n_rows, n_cols = self._mt.count()
        n_pheno = self._phenotype_ht.count()
        total_gene_sets = sum(len(c) for c in self._collections.values())
        all_genes = set()
        for coll in self._collections.values():
            for gs in coll:
                all_genes.update(gs.genes)
        logger.info(
            "Loaded: %d variants x %d samples, %d phenotypes, "
            "%d gene sets (%d unique genes) across %d collections",
            n_rows,
            n_cols,
            n_pheno,
            total_gene_sets,
            len(all_genes),
            len(self._collections),
        )

        # Determine iteration
        vc_items: Dict[str, Any]
        if self.config.variant_classes:
            vc_items = dict(self.config.variant_classes)
        else:
            vc_items = {"all": None}

        total_runs = len(vc_items) * len(self._collections)
        run_idx = 0

        t_runs = time.time()
        for vc_name, vf in vc_items.items():
            for coll_name, coll in self._collections.items():
                run_idx += 1
                logger.info(
                    "\n[Run %d/%d] variant_class=%s, collection=%s",
                    run_idx,
                    total_runs,
                    vc_name,
                    coll_name,
                )
                gene_sets_dict = {gs.name: list(gs.genes) for gs in coll}
                self._run_single(vc_name, vf, coll_name, gene_sets_dict)
        logger.info("All %d runs completed in %.1fs", total_runs, time.time() - t_runs)

        combined_df = self._generate_summary()

        if self.config.generate_report:
            self._generate_outputs(combined_df)

        elapsed = time.time() - t0
        logger.info("=" * 70)
        logger.info("BURDEN PIPELINE — COMPLETE (%.1fs)", elapsed)
        logger.info("=" * 70)

        return combined_df

    # ------------------------------------------------------------------
    # Internal stages
    # ------------------------------------------------------------------

    def _initialize_hail(self) -> None:
        if hl is None:  # pragma: no cover
            raise ImportError(
                "Hail is required for BurdenPipeline. "
                "Install hvantk with the 'hail' extra."
            )
        try:
            hl.current_backend()
            logger.info("Hail already initialised")
        except Exception:
            from hvantk.core.hail_context import init_hail

            init_hail()

    def _load_cohort(self) -> None:
        logger.info("Loading cohort MT: %s", self.config.cohort_mt_path)
        self._mt = hl.read_matrix_table(self.config.cohort_mt_path)
        n_rows, n_cols = self._mt.count()
        logger.info("  %d variants, %d samples", n_rows, n_cols)

    @staticmethod
    def _resolve_nested_field(obj: Any, field_path: str) -> Any:
        """Resolve a field path on a Hail table/struct.

        When ``field_path`` contains dots (e.g. ``"phe.is_case"``), this
        navigates the nested struct: ``obj["phe"]["is_case"]``.

        This is necessary because Hail's ``obj["phe.is_case"]`` looks up
        a single field literally named ``"phe.is_case"`` — it does NOT
        traverse into a struct named ``phe``.  Our convention is that
        dots always mean struct navigation.
        """
        expr = obj
        for part in field_path.split("."):
            expr = expr[part]
        return expr

    @staticmethod
    def _leaf_name(field_path: str) -> str:
        """Return the last component of a dot-delimited path."""
        return field_path.rsplit(".", 1)[-1]

    def _load_phenotypes(self) -> None:
        path = self.config.phenotype_ht_path

        if not path:
            # Extract phenotype and covariates from MT column fields
            logger.info("Extracting phenotype from MT column fields")
            pheno_field = self.config.phenotype_field
            cov_fields = self.config.covariate_fields or []

            # Get cols table first, then resolve fields against it
            cols_ht = self._mt.cols()

            annotations = {}
            pheno_leaf = self._leaf_name(pheno_field)
            annotations[pheno_leaf] = self._resolve_nested_field(cols_ht, pheno_field)

            for cov in cov_fields:
                cov_leaf = self._leaf_name(cov)
                annotations[cov_leaf] = self._resolve_nested_field(cols_ht, cov)

            self._phenotype_ht = cols_ht.select(**annotations)

            # Update config fields to use flattened leaf names so downstream
            # regression uses the correct field names in the phenotype HT.
            self.config.phenotype_field = pheno_leaf
            self.config.covariate_fields = [self._leaf_name(c) for c in cov_fields]

            logger.info(
                "  Extracted phenotype '%s' and %d covariates from MT cols",
                pheno_leaf,
                len(cov_fields),
            )
        elif path.endswith(".ht"):
            logger.info("Loading phenotypes: %s", path)
            self._phenotype_ht = hl.read_table(path)
        else:
            logger.info("Loading phenotypes: %s", path)
            self._phenotype_ht = hl.import_table(path, impute=True).key_by(
                self.config.sample_id_field
            )

        logger.info("  %d samples", self._phenotype_ht.count())

    def _load_gene_sets(self) -> None:
        from hvantk.utils.gene_sets import GeneSetCollection

        all_genes: set = set()
        for name, path in self.config.gene_set_collections.items():
            logger.info("Loading gene set collection '%s': %s", name, path)
            coll = GeneSetCollection.load(path)
            logger.info("  %d gene sets", len(coll))
            for gs in coll:
                all_genes.update(gs.genes)
            self._collections[name] = coll
        logger.info(
            "Total: %d unique genes across all collections",
            len(all_genes),
        )

    def _load_gene_lengths(self) -> None:
        if not self.config.gene_lengths_path:
            return
        logger.info("Loading gene lengths: %s", self.config.gene_lengths_path)
        df = pd.read_csv(self.config.gene_lengths_path, sep="\t")
        self._gene_lengths = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        logger.info("  %d gene lengths loaded", len(self._gene_lengths))

    def _run_single(
        self,
        vc_name: str,
        variant_filter: Any,
        coll_name: str,
        gene_sets_dict: Dict[str, List[str]],
    ) -> None:
        from hvantk.enrichex.burden import (
            permutation_burden_test,
            run_burden_analysis,
        )
        from hvantk.enrichex.correction import apply_correction

        t_run = time.time()

        result_ht = run_burden_analysis(
            cohort_mt=self._mt,
            gene_sets=gene_sets_dict,
            phenotype_ht=self._phenotype_ht,
            phenotype_field=self.config.phenotype_field,
            covariate_fields=self.config.covariate_fields or None,
            phenotype_type=self.config.phenotype_type,
            gene_field=self.config.gene_field,
            genotype_aggregation=self.config.genotype_aggregation,
            variant_filter=variant_filter,
            normalize_by_length=self.config.normalize_by_length,
            gene_lengths=self._gene_lengths,
            min_carriers=self.config.min_carriers,
        )

        if result_ht is None:
            logger.warning(
                "No results for variant_class=%s, collection=%s (%.1fs)",
                vc_name,
                coll_name,
                time.time() - t_run,
            )
            self._run_results.append(
                BurdenRunResult(
                    variant_class=vc_name,
                    collection_name=coll_name,
                    results_df=pd.DataFrame(),
                )
            )
            return

        df = result_ht.to_pandas()

        # Coerce p_value to float (Hail may produce pandas NA)
        df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

        # Apply correction only to non-NaN p-values
        valid_mask = df["p_value"].notna()
        df["p_adjusted"] = float("nan")
        if valid_mask.any():
            p_adjusted = apply_correction(
                df.loc[valid_mask, "p_value"].tolist(),
                method=self.config.correction_method,
            )
            df.loc[valid_mask, "p_adjusted"] = p_adjusted
        df["significant"] = df["p_adjusted"] < self.config.alpha
        df["variant_class"] = vc_name
        df["collection"] = coll_name
        df = df.sort_values("p_value")

        # Write per-run TSV
        run_path = self._per_run_dir / f"burden_{coll_name}_{vc_name}.tsv"
        df.to_csv(run_path, sep="\t", index=False)
        logger.info("  Wrote %s (%d gene sets)", run_path.name, len(df))

        # Optional permutation test
        perm_df = None
        if self.config.competitive:
            logger.info("  Running competitive permutation test...")
            perm_df = permutation_burden_test(
                cohort_mt=self._mt,
                gene_sets=gene_sets_dict,
                phenotype_ht=self._phenotype_ht,
                phenotype_field=self.config.phenotype_field,
                covariate_fields=self.config.covariate_fields or None,
                phenotype_type=self.config.phenotype_type,
                gene_field=self.config.gene_field,
                genotype_aggregation=self.config.genotype_aggregation,
                variant_filter=variant_filter,
                n_permutations=self.config.n_permutations,
                seed=self.config.permutation_seed,
                normalize_by_length=self.config.normalize_by_length,
                gene_lengths=self._gene_lengths,
            )
            if not perm_df.empty:
                perm_df["variant_class"] = vc_name
                perm_df["collection"] = coll_name
                perm_path = self._per_run_dir / f"permutation_{coll_name}_{vc_name}.tsv"
                perm_df.to_csv(perm_path, sep="\t", index=False)
                logger.info("  Wrote %s", perm_path.name)

        n_sig = int(df["significant"].sum()) if not df.empty else 0
        run_elapsed = time.time() - t_run
        logger.info(
            "  Run complete: %d gene sets tested, %d significant (%.1fs)",
            len(df),
            n_sig,
            run_elapsed,
        )
        self._run_results.append(
            BurdenRunResult(
                variant_class=vc_name,
                collection_name=coll_name,
                results_df=df,
                permutation_df=perm_df,
                n_gene_sets_tested=len(df),
                n_significant=n_sig,
            )
        )

    # ------------------------------------------------------------------
    # Summary and outputs
    # ------------------------------------------------------------------

    def _generate_summary(self) -> pd.DataFrame:
        """Combine all run results and write summary."""
        dfs = [r.results_df for r in self._run_results if not r.results_df.empty]
        if not dfs:
            logger.warning("No runs produced results.")
            return pd.DataFrame()

        combined = pd.concat(dfs, ignore_index=True)
        summary_path = self._output_dir / "burden_combined.tsv"
        combined.to_csv(summary_path, sep="\t", index=False)
        logger.info("Combined results: %s (%d rows)", summary_path, len(combined))

        # Print summary table
        self._print_summary()

        return combined

    def _print_summary(self) -> None:
        """Print a summary table to the log and stdout."""
        lines = [
            "",
            "=" * 70,
            "BURDEN PIPELINE SUMMARY",
            "=" * 70,
            f"  {'Variant Class':<25} {'Collection':<20} {'Tested':>7} {'Sig':>5}",
            "  " + "-" * 60,
        ]
        for r in self._run_results:
            lines.append(
                f"  {r.variant_class:<25} {r.collection_name:<20} "
                f"{r.n_gene_sets_tested:>7} {r.n_significant:>5}"
            )
        lines.append("=" * 70)
        summary = "\n".join(lines)
        logger.info(summary)
        print(summary)

    def _generate_outputs(self, combined_df: pd.DataFrame) -> None:
        """Generate HTML report from combined results."""
        if combined_df.empty:
            logger.warning("No results to report.")
            return

        # Write combined TSV to output dir (already done in _generate_summary)
        results_path = self._output_dir / "burden_combined.tsv"

        try:
            from hvantk.enrichex.report import generate_report

            report_path = self._output_dir / "enrichex_burden_report.html"
            generate_report(
                output_path=str(report_path),
                burden_results=str(results_path),
                top_n=25,
                embed_static_plots=True,
            )
            logger.info("Report: %s", report_path)
        except Exception:
            logger.warning("Report generation failed.", exc_info=True)
