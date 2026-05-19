"""
PSROC Pipeline Module

This module provides orchestration for PSROC (Prediction Score ROC Analysis)
workflows, from loading variant data through score annotation, ROC computation,
and report generation.

The pipeline evaluates the discriminative power of variant pathogenicity
prediction scores (e.g., CADD, REVEL, MetaLR) against ClinVar truth labels.

Example:
    >>> from hvantk.algorithms.psroc.pipeline import PSROCConfig, PSROCPipeline
    >>> config = PSROCConfig(
    ...     clinvar_ht="/data/clinvar.ht",
    ...     dbnsfp_ht="/data/dbnsfp.ht",
    ...     scores=["CADD_phred", "REVEL_score"],
    ...     output_dir="/results/psroc",
    ... )
    >>> pipeline = PSROCPipeline(config)
    >>> result = pipeline.run()
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any, List, Set
from enum import Enum
import hashlib
import json
import logging
import re
from datetime import datetime

import numpy as np
import hail as hl

from hvantk.algorithms.psroc.roc import (
    ROCResult,
    ScoreMissingness,
    compute_roc_metrics,
    compute_all_missingness,
    filter_scores_by_missingness,
)
from hvantk.algorithms.psroc.plots import (
    plot_roc_curves,
    plot_auc_comparison,
    plot_missingness_summary,
    plot_psroc_summary_dashboard,
    plot_collection_heatmap,
)
from hvantk.core.utils.gene_sets import load_gene_set

logger = logging.getLogger(__name__)


# Re-export from core for backward compatibility
from hvantk.core.constants import (  # noqa: E402
    CLINVAR_PATHOGENIC_LABELS as PATHOGENIC_LABELS,
    CLINVAR_BENIGN_LABELS as BENIGN_LABELS,
)

# Score directionality: True means higher values indicate pathogenicity.
# Scores not in this map default to higher_is_pathogenic=True.
SCORE_DIRECTIONALITY: Dict[str, bool] = {
    # Higher = pathogenic
    "CADD_phred": True,
    "CADD_raw": True,
    "REVEL_score": True,
    "MetaLR_score": True,
    "MetaSVM_score": True,
    "MetaRNN_score": True,
    "M-CAP_score": True,
    "MPC_score": True,
    "PrimateAI_score": True,
    "DEOGEN2_score": True,
    "BayesDel_addAF_score": True,
    "BayesDel_noAF_score": True,
    "ClinPred_score": True,
    "VEST4_score": True,
    "Eigen-raw_coding": True,
    "Eigen-PC-raw_coding": True,
    "GenoCanyon_score": True,
    "integrated_fitCons_score": True,
    "GM12878_fitCons_score": True,
    "H1-hESC_fitCons_score": True,
    "HUVEC_fitCons_score": True,
    "LINSIGHT": True,
    "GERP++_RS": True,
    "phyloP100way_vertebrate": True,
    "phyloP30way_mammalian": True,
    "phyloP17way_primate": True,
    "phastCons100way_vertebrate": True,
    "phastCons30way_mammalian": True,
    "phastCons17way_primate": True,
    "fathmm-MKL_coding_score": True,
    "fathmm-XF_coding_score": True,
    "MutationAssessor_score": True,
    "MutPred_score": True,
    "MVP_score": True,
    "gMVP_score": True,
    # Lower = pathogenic
    "SIFT_score": False,
    "SIFT4G_score": False,
    "PROVEAN_score": False,
    "FATHMM_score": False,
    "LRT_score": False,
}

# ClinVar CLNREVSTAT string → review-star mapping
# See: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
CLNREVSTAT_STAR_MAP = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_for_the_individual_variant": 0,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_assertion_provided": 0,
    "no_classifications_from_unflagged_records": 0,
    "flagged_submission": 0,
}


class PSROCStage(Enum):
    """Enumeration of PSROC pipeline stages."""

    LOAD_TABLES = "load_tables"
    FILTER_CLINVAR = "filter_clinvar"
    ASSIGN_LABELS = "assign_labels"
    ANNOTATE_SCORES = "annotate_scores"
    COMPUTE_MISSINGNESS = "compute_missingness"
    COMPUTE_ROC = "compute_roc"
    GENERATE_OUTPUTS = "generate_outputs"


@dataclass
class PSROCConfig:
    """Configuration for PSROC pipeline execution.

    Attributes:
        genes: List of gene symbols to filter ClinVar variants.
        genes_file: Path to a file containing gene symbols (one per line).
        variants_path: Path to a variant list file (chr:pos:ref:alt format).
        clinvar_ht: Path to pre-built ClinVar Hail Table.
        dbnsfp_ht: Path to pre-built dbNSFP Hail Table.
        scores: List of dbNSFP score field names to evaluate.
        output_dir: Directory for output files.
        reference_genome: Reference genome (default: GRCh38).
        min_stars: Minimum ClinVar review status stars (default: 1).
        max_missingness: Maximum allowed missingness rate per score (default: 0.3).
        threshold_method: Method for finding optimal threshold (default: "youden").
        hgnc_path: Path to HGNC data (TSV or .ht) for gene alias resolution.
        export_tsv: Whether to export annotated variants as TSV.
        overwrite: Whether to overwrite existing output files.
        generate_plots: Whether to generate visualization plots.
        output_prefix: Prefix for output file names.
    """

    # Input sources (one of genes/genes_file/variants_path/gene_set_collection required)
    genes: Optional[List[str]] = None
    genes_file: Optional[str] = None
    variants_path: Optional[str] = None
    gene_set_collection: Optional[Dict[str, Set[str]]] = None

    # Required table paths
    clinvar_ht: str = ""
    dbnsfp_ht: str = ""

    # Score configuration
    scores: List[str] = field(default_factory=list)
    score_directions: Optional[
        Dict[str, bool]
    ] = None  # per-score override: True=higher_is_pathogenic

    # Output configuration
    output_dir: str = ""
    output_prefix: str = "psroc"

    # Processing options
    reference_genome: str = "GRCh38"
    min_stars: int = 1
    max_missingness: float = 0.3
    threshold_method: str = "youden"
    hgnc_path: Optional[str] = None

    # Output options
    export_tsv: bool = False
    overwrite: bool = False
    generate_plots: bool = True
    group_name: Optional[str] = None

    # Minimum number of labeled (P + B) variants required for ROC analysis
    min_variants: int = 10

    # Bootstrap confidence interval for AUC
    n_bootstrap: int = 2000

    def validate(self) -> List[str]:
        """Validate configuration and return list of errors.

        Returns:
            List of validation error messages (empty if valid).
        """
        errors = []

        # Check variant source
        sources = [
            self.genes is not None and len(self.genes) > 0,
            self.genes_file is not None,
            self.variants_path is not None,
            self.gene_set_collection is not None and len(self.gene_set_collection) > 0,
        ]
        if sum(sources) == 0:
            errors.append(
                "Must provide at least one of: --genes, --genes-file, "
                "--variants, or --gene-sets"
            )
        if sum(sources) > 1:
            errors.append(
                "Cannot provide multiple variant sources. Use only one of: "
                "--genes, --genes-file, --variants, or --gene-sets"
            )

        # Check required table paths
        if not self.clinvar_ht:
            errors.append("--clinvar-ht is required")
        elif not Path(self.clinvar_ht).exists():
            errors.append(f"ClinVar table not found: {self.clinvar_ht}")

        if not self.dbnsfp_ht:
            errors.append("--dbnsfp-ht is required")
        elif not Path(self.dbnsfp_ht).exists():
            errors.append(f"dbNSFP table not found: {self.dbnsfp_ht}")

        # Check scores
        if not self.scores:
            errors.append("Must provide at least one score via --scores")

        # Check output directory
        if not self.output_dir:
            errors.append("--output-dir is required")

        # Check file paths
        if self.genes_file and not Path(self.genes_file).exists():
            errors.append(f"Genes file not found: {self.genes_file}")

        if self.variants_path and not Path(self.variants_path).exists():
            errors.append(f"Variants file not found: {self.variants_path}")

        if self.gene_set_collection is not None:
            empty_groups = [
                name for name, genes in self.gene_set_collection.items() if not genes
            ]
            if empty_groups:
                errors.append(
                    f"Gene set collection contains empty groups: "
                    f"{', '.join(empty_groups)}"
                )

        if self.hgnc_path and not Path(self.hgnc_path).exists():
            errors.append(f"HGNC data not found: {self.hgnc_path}")

        # Validate thresholds
        if not 0.0 <= self.max_missingness <= 1.0:
            errors.append("max_missingness must be between 0.0 and 1.0")

        if self.min_stars < 0:
            errors.append("min_stars must be non-negative")

        # Validate threshold method
        valid_methods = {"youden", "closest_to_corner", "f1"}
        if self.threshold_method not in valid_methods:
            errors.append(
                f"Invalid threshold_method: {self.threshold_method}. "
                f"Must be one of: {valid_methods}"
            )

        # Validate min_variants and n_bootstrap
        if self.min_variants < 1:
            errors.append("min_variants must be >= 1")
        if self.n_bootstrap < 0:
            errors.append("n_bootstrap must be >= 0")

        return errors

    def get_gene_set(self) -> Optional[Set[str]]:
        """Load and return the gene set from configured sources.

        Returns:
            Set of gene symbols, or None if using variants_path.
        """
        if self.variants_path:
            return None

        gene_set: Set[str] = set()

        if self.genes:
            gene_set.update(self.genes)

        if self.genes_file:
            gene_set.update(load_gene_set(path=self.genes_file))

        return gene_set if gene_set else None


@dataclass
class PSROCState:
    """Tracks the state of PSROC pipeline execution."""

    config: PSROCConfig
    current_stage: Optional[PSROCStage] = None
    completed_stages: List[str] = field(default_factory=list)
    outputs: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    start_time: Optional[str] = None
    end_time: Optional[str] = None

    @staticmethod
    def _serialize_outputs(outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Convert pipeline outputs to JSON-serializable form.

        Handles ScoreMissingness, ROCResult, and dicts thereof by calling
        their ``.to_dict()`` method.  Scalars and strings pass through
        unchanged.
        """
        serialized: Dict[str, Any] = {}
        for key, value in outputs.items():
            if hasattr(value, "to_dict"):
                serialized[key] = value.to_dict()
            elif isinstance(value, dict):
                serialized[key] = {
                    k: v.to_dict() if hasattr(v, "to_dict") else v
                    for k, v in value.items()
                }
            else:
                serialized[key] = value
        return serialized

    def save(self, path: Path) -> None:
        """Save pipeline state to JSON file.

        Args:
            path: Path to save state file.
        """
        state_dict = {
            "config": asdict(self.config),
            "current_stage": self.current_stage.value if self.current_stage else None,
            "completed_stages": self.completed_stages,
            "outputs": self._serialize_outputs(self.outputs),
            "errors": self.errors,
            "start_time": self.start_time,
            "end_time": self.end_time,
        }

        with open(path, "w") as f:
            json.dump(state_dict, f, indent=2)

        logger.info(f"Pipeline state saved to {path}")

    @classmethod
    def load(cls, path: Path) -> "PSROCState":
        """Load pipeline state from JSON file.

        Args:
            path: Path to state file.

        Returns:
            PSROCState object.
        """
        with open(path, "r") as f:
            state_dict = json.load(f)

        config = PSROCConfig(**state_dict["config"])

        state = cls(
            config=config,
            current_stage=(
                PSROCStage(state_dict["current_stage"])
                if state_dict["current_stage"]
                else None
            ),
            completed_stages=state_dict["completed_stages"],
            outputs=state_dict["outputs"],
            errors=state_dict["errors"],
            start_time=state_dict.get("start_time"),
            end_time=state_dict.get("end_time"),
        )

        logger.info(f"Pipeline state loaded from {path}")
        return state

    def mark_stage_complete(
        self, stage: PSROCStage, output_path: Optional[str] = None
    ) -> None:
        """Mark a stage as completed."""
        stage_name = stage.value
        if stage_name not in self.completed_stages:
            self.completed_stages.append(stage_name)

        if output_path:
            self.outputs[stage_name] = output_path

        logger.info(f"Stage {stage_name} marked as complete")

    def is_stage_complete(self, stage: PSROCStage) -> bool:
        """Check if a stage has been completed."""
        return stage.value in self.completed_stages


@dataclass
class PSROCResult:
    """Results from PSROC pipeline execution.

    Attributes:
        annotated_ht_path: Path to the annotated Hail Table.
        metrics: ROC results for scores that passed missingness threshold.
        missingness: Missingness statistics for all requested scores.
        n_pathogenic: Number of pathogenic variants.
        n_benign: Number of benign variants.
        n_excluded: Number of excluded variants (uncertain/conflicting).
        n_total: Total number of variants processed.
        n_out_of_scope: Variants dropped with no dbNSFP scores (non-missense).
        n_genes: Number of genes in the panel.
        scores_included: Scores that passed the missingness threshold.
        scores_excluded: Scores excluded due to high missingness.
        max_missingness_threshold: The threshold used for this run.
        output_dir: Path to the output directory.
    """

    annotated_ht_path: str
    metrics: Dict[str, ROCResult]
    missingness: Dict[str, ScoreMissingness]
    n_pathogenic: int
    n_benign: int
    n_excluded: int
    n_total: int
    scores_included: List[str]
    scores_excluded: List[str]
    max_missingness_threshold: float
    output_dir: str
    n_out_of_scope: int = 0
    n_genes: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert results to a dictionary for JSON serialization."""
        return {
            "annotated_ht_path": self.annotated_ht_path,
            "metrics": {name: r.to_dict() for name, r in self.metrics.items()},
            "missingness": {name: m.to_dict() for name, m in self.missingness.items()},
            "n_pathogenic": self.n_pathogenic,
            "n_benign": self.n_benign,
            "n_excluded": self.n_excluded,
            "n_total": self.n_total,
            "n_out_of_scope": self.n_out_of_scope,
            "n_genes": self.n_genes,
            "scores_included": self.scores_included,
            "scores_excluded": self.scores_excluded,
            "max_missingness_threshold": self.max_missingness_threshold,
            "output_dir": self.output_dir,
        }

    def summary(self) -> str:
        """Generate a human-readable summary of results."""
        lines = [
            "=" * 60,
            "PSROC Analysis Summary",
            "=" * 60,
            "",
            f"Gene panel size: {self.n_genes}",
            f"Variants analyzed: {self.n_total}",
            f"  Pathogenic (P): {self.n_pathogenic}",
            f"  Benign (B): {self.n_benign}",
            f"  B/P ratio: {self.n_benign / self.n_pathogenic:.1f}:1"
            if self.n_pathogenic > 0
            else "  B/P ratio: N/A",
            f"  Excluded: {self.n_excluded}",
            f"  Out of scope (no dbNSFP scores): {self.n_out_of_scope}",
            "",
            f"Scores requested: {len(self.missingness)}",
            f"  Included: {len(self.scores_included)}",
            f"  Excluded: {len(self.scores_excluded)}",
            f"  Max missingness threshold: {self.max_missingness_threshold:.0%}",
            "",
        ]

        if self.metrics:
            lines.append("ROC Results (sorted by AUC):")
            lines.append("-" * 40)
            sorted_metrics = sorted(
                self.metrics.items(), key=lambda x: x[1].auc, reverse=True
            )
            for name, roc in sorted_metrics:
                lines.append(
                    f"  {name}: AUC={roc.auc:.3f}, "
                    f"threshold={roc.optimal_threshold:.3f}, "
                    f"sens={roc.sensitivity_at_optimal:.3f}, "
                    f"spec={roc.specificity_at_optimal:.3f}"
                )
            lines.append("")

        if self.scores_excluded:
            lines.append("Excluded scores (high missingness):")
            for name in self.scores_excluded:
                miss = self.missingness[name]
                lines.append(f"  {name}: {miss.missingness_rate:.1%} missing")
            lines.append("")

        lines.append(f"Output directory: {self.output_dir}")
        lines.append("=" * 60)

        return "\n".join(lines)


def parse_variant_list(variants_path: str):
    """Parse a variant list file (chr:pos:ref:alt per line) and normalize chromosome and position."""
    variant_keys = []
    with open(variants_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(":")
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
                # Normalize chromosome (add chr prefix if needed)
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                variant_keys.append((chrom, int(pos), ref, alt))
    return variant_keys


class PSROCPipeline:
    """PSROC: Prediction Score ROC Analysis Pipeline.

    End-to-end workflow for evaluating variant pathogenicity prediction scores:
        1. Load ClinVar and dbNSFP Hail Tables
        2. Filter ClinVar to genes/variants of interest
        3. Assign binary labels (pathogenic/benign)
        4. Annotate with dbNSFP scores
        5. Compute missingness statistics per score
        6. Filter scores by missingness threshold
        7. Compute ROC metrics for included scores
        8. Generate plots, metrics, and missingness reports

    Example:
        >>> config = PSROCConfig(
        ...     genes=["BRCA1", "BRCA2"],
        ...     clinvar_ht="/data/clinvar.ht",
        ...     dbnsfp_ht="/data/dbnsfp.ht",
        ...     scores=["CADD_phred", "REVEL_score"],
        ...     output_dir="/results/psroc",
        ... )
        >>> pipeline = PSROCPipeline(config)
        >>> result = pipeline.run()
    """

    def __init__(self, config: PSROCConfig):
        """Initialize PSROC pipeline.

        Args:
            config: Pipeline configuration.

        Raises:
            ValueError: If configuration validation fails.
        """
        self.config = config
        errors = config.validate()
        if errors:
            raise ValueError(
                "Configuration validation failed:\n"
                + "\n".join(f"  - {e}" for e in errors)
            )

        self.state = PSROCState(config=config)
        self._setup_output_paths()
        self._initialize_hail()

        # Internal state
        self._clinvar_ht: Optional[hl.Table] = None
        self._dbnsfp_ht: Optional[hl.Table] = None
        self._labeled_ht: Optional[hl.Table] = None
        self._annotated_ht: Optional[hl.Table] = None
        self._n_out_of_scope: int = 0
        self._n_genes: int = 0

    def _setup_output_paths(self) -> None:
        """Initialize output directory structure."""
        output_dir = Path(self.config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix = self.config.output_prefix
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        logs_dir = output_dir / "logs"
        logs_dir.mkdir(exist_ok=True)

        self.paths = {
            "annotated_ht": str(output_dir / f"{prefix}_annotated.ht"),
            "annotated_tsv": str(output_dir / f"{prefix}_annotated.tsv"),
            "metrics_json": str(output_dir / f"{prefix}_metrics.json"),
            "missingness_json": str(output_dir / f"{prefix}_missingness.json"),
            "roc_curves_png": str(plots_dir / f"{prefix}_roc_curves.png"),
            "auc_comparison_png": str(plots_dir / f"{prefix}_auc_comparison.png"),
            "missingness_png": str(plots_dir / f"{prefix}_missingness.png"),
            "dashboard_png": str(plots_dir / f"{prefix}_dashboard.png"),
            "state": str(output_dir / ".pipeline_state.json"),
            "log": str(
                logs_dir / f'psroc_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            ),
        }

        logger.info(f"Output paths configured in: {output_dir}")

    def _initialize_hail(self) -> None:
        """Initialize Hail if not already initialized."""
        try:
            hl.current_backend()
            logger.info("Hail already initialized")
        except Exception:
            from hvantk.core.hail_context import init_hail

            init_hail()
            logger.info("Hail initialized via hvantk")

    def show_plan(self) -> None:
        """Display the execution plan without running."""
        print("\n" + "=" * 70)
        print("PSROC PIPELINE EXECUTION PLAN")
        print("=" * 70)

        print("\n📁 Input Sources:")
        if self.config.genes:
            print(f"  Genes: {', '.join(self.config.genes[:5])}", end="")
            if len(self.config.genes) > 5:
                print(f" ... (+{len(self.config.genes) - 5} more)")
            else:
                print()
        elif self.config.genes_file:
            print(f"  Genes file: {self.config.genes_file}")
        elif self.config.variants_path:
            print(f"  Variants file: {self.config.variants_path}")

        print(f"  ClinVar table: {self.config.clinvar_ht}")
        print(f"  dbNSFP table: {self.config.dbnsfp_ht}")

        print("\n🔧 Configuration:")
        print(f"  Scores: {', '.join(self.config.scores)}")
        print(f"  Reference genome: {self.config.reference_genome}")
        print(f"  Min review stars: {self.config.min_stars}")
        print(f"  Max missingness: {self.config.max_missingness:.0%}")
        print(f"  Threshold method: {self.config.threshold_method}")
        if self.config.hgnc_path:
            print(f"  HGNC alias resolution: {self.config.hgnc_path}")

        print("\n📊 Stages to execute:")
        stages = [
            ("1", "Load Hail Tables", "Read ClinVar and dbNSFP tables"),
            ("2", "Filter ClinVar", "Filter to target genes/variants"),
            ("3", "Assign Labels", "Convert CLNSIG to binary P/B labels"),
            ("4", "Annotate Scores", "Join with dbNSFP prediction scores"),
            ("5", "Compute Missingness", "Calculate per-score missingness"),
            ("6", "Compute ROC", "Calculate ROC metrics for included scores"),
            ("7", "Generate Outputs", "Create plots, metrics, and reports"),
        ]

        for num, name, desc in stages:
            print(f"  [{num}] {name}")
            print(f"      {desc}")

        print("\n📂 Output Directory:")
        print(f"  {self.config.output_dir}")
        if self.config.generate_plots:
            print("  Plots: enabled")
        if self.config.export_tsv:
            print("  TSV export: enabled")

        print("\n" + "=" * 70 + "\n")

    def run(self) -> PSROCResult:
        """Execute the PSROC pipeline for a single gene set.

        For multi-group analysis with ``gene_set_collection``, use
        ``run_collection()`` instead.

        Returns:
            PSROCResult with computed metrics and output paths.

        Raises:
            ValueError: If gene_set_collection is the active source.
            Exception: If pipeline execution fails.
        """
        if self.config.gene_set_collection:
            raise ValueError(
                "Config uses gene_set_collection. "
                "Use run_collection() for multi-group execution."
            )

        self.state.start_time = datetime.now().isoformat()
        logger.info("Starting PSROC pipeline execution")

        try:
            # Stage 1: Load tables
            self._run_stage(PSROCStage.LOAD_TABLES)

            # Stage 2: Filter ClinVar
            self._run_stage(PSROCStage.FILTER_CLINVAR)

            # Stage 3: Assign labels
            self._run_stage(PSROCStage.ASSIGN_LABELS)

            # Stage 4: Annotate scores
            self._run_stage(PSROCStage.ANNOTATE_SCORES)

            # Stage 5: Compute missingness
            self._run_stage(PSROCStage.COMPUTE_MISSINGNESS)

            # Stage 6: Compute ROC metrics
            self._run_stage(PSROCStage.COMPUTE_ROC)

            # Stage 7: Generate outputs
            result = self._run_stage(PSROCStage.GENERATE_OUTPUTS)

            self.state.end_time = datetime.now().isoformat()
            logger.info("PSROC pipeline execution completed successfully")

            return result

        except Exception as e:
            self.state.errors.append(str(e))
            self.state.end_time = datetime.now().isoformat()
            logger.exception(f"PSROC pipeline execution failed: {e}")
            raise

        finally:
            self.state.save(Path(self.paths["state"]))

    def run_collection(self) -> Dict[str, PSROCResult]:
        """Execute PSROC pipeline independently for each named gene set.

        Requires ``gene_set_collection`` to be set in the config.  For each
        group, a per-group ``PSROCConfig`` is created with the group's genes
        and a dedicated output subdirectory, then the standard ``run()``
        pipeline is executed.

        Groups that yield zero ClinVar variants (or otherwise fail) are
        skipped with a warning.

        Returns:
            Dictionary mapping group names to their PSROCResult.

        Raises:
            ValueError: If ``gene_set_collection`` is not configured.
            RuntimeError: If all groups fail.
        """
        collection = self.config.gene_set_collection
        if not collection:
            raise ValueError(
                "gene_set_collection is not configured. "
                "Use run() for single gene set execution."
            )

        results: Dict[str, PSROCResult] = {}
        failed_groups: List[str] = []

        for group_name in sorted(collection):
            gene_set = collection[group_name]
            logger.info(
                f"Running PSROC for group '{group_name}' " f"({len(gene_set)} genes)"
            )

            # Sanitize group name for filesystem paths
            safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", group_name)
            safe_name = re.sub(r"_+", "_", safe_name).strip("_")
            if not safe_name or safe_name in (".", ".."):
                slug = hashlib.sha256(group_name.encode()).hexdigest()[:8]
                safe_name = f"group_{slug}"

            group_config = PSROCConfig(
                genes=sorted(gene_set),
                clinvar_ht=self.config.clinvar_ht,
                dbnsfp_ht=self.config.dbnsfp_ht,
                scores=list(self.config.scores),
                output_dir=str(Path(self.config.output_dir) / safe_name),
                output_prefix=f"{self.config.output_prefix}_{safe_name}",
                reference_genome=self.config.reference_genome,
                min_stars=self.config.min_stars,
                max_missingness=self.config.max_missingness,
                threshold_method=self.config.threshold_method,
                hgnc_path=self.config.hgnc_path,
                export_tsv=self.config.export_tsv,
                overwrite=self.config.overwrite,
                generate_plots=self.config.generate_plots,
                group_name=group_name,
                min_variants=self.config.min_variants,
                n_bootstrap=self.config.n_bootstrap,
            )

            try:
                pipeline = PSROCPipeline(group_config)
                result = pipeline.run()
                results[group_name] = result
                logger.info(
                    f"Group '{group_name}' completed: "
                    f"{result.n_total} variants, "
                    f"{len(result.scores_included)} scores"
                )
            except Exception as e:
                logger.warning(f"Group '{group_name}' failed: {e}")
                failed_groups.append(group_name)
                continue

        if not results:
            raise RuntimeError(
                "All gene set groups failed. "
                f"Failed groups: {', '.join(failed_groups)}"
            )

        if failed_groups:
            logger.warning(
                f"{len(failed_groups)} group(s) failed: " f"{', '.join(failed_groups)}"
            )

        logger.info(
            f"Gene set collection complete: "
            f"{len(results)}/{len(collection)} groups succeeded"
        )

        # Generate cross-panel heatmap if plots are enabled and we have results
        if self.config.generate_plots and len(results) >= 2:
            self._generate_collection_heatmap(results)

        return results

    def _generate_collection_heatmap(self, results: Dict[str, "PSROCResult"]) -> None:
        """Generate a cross-panel AUC heatmap from run_collection results."""
        import matplotlib

        matplotlib.use("Agg")

        collection_metrics = {
            name: result.metrics for name, result in results.items() if result.metrics
        }

        if len(collection_metrics) < 2:
            logger.info("Skipping collection heatmap: fewer than 2 groups with metrics")
            return

        output_dir = Path(self.config.output_dir)
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        heatmap_path = str(
            plots_dir / f"{self.config.output_prefix}_collection_heatmap"
        )

        try:
            plot_collection_heatmap(
                collection_metrics,
                output_path=heatmap_path,
                title="AUC Across Gene Set Panels",
            )
            logger.info(f"Collection heatmap saved: {heatmap_path}.png")
        except Exception as e:
            logger.warning(f"Collection heatmap generation failed: {e}")

    def _run_stage(self, stage: PSROCStage) -> Any:
        """Execute a single pipeline stage.

        Args:
            stage: Pipeline stage to execute.

        Returns:
            Stage-specific output.

        Raises:
            Exception: If stage execution fails.
        """
        self.state.current_stage = stage
        logger.info(f"Running stage: {stage.value}")

        try:
            if stage == PSROCStage.LOAD_TABLES:
                output = self._load_tables()
            elif stage == PSROCStage.FILTER_CLINVAR:
                output = self._filter_clinvar()
            elif stage == PSROCStage.ASSIGN_LABELS:
                output = self._assign_labels()
            elif stage == PSROCStage.ANNOTATE_SCORES:
                output = self._annotate_scores()
            elif stage == PSROCStage.COMPUTE_MISSINGNESS:
                output = self._compute_missingness()
            elif stage == PSROCStage.COMPUTE_ROC:
                output = self._compute_roc_metrics()
            elif stage == PSROCStage.GENERATE_OUTPUTS:
                output = self._generate_outputs()
            else:
                raise ValueError(f"Unknown stage: {stage}")

            self.state.mark_stage_complete(stage)
            return output

        except Exception as e:
            error_msg = f"Stage {stage.value} failed: {e}"
            self.state.errors.append(error_msg)
            logger.exception(error_msg)
            raise

    def _load_tables(self) -> None:
        """Stage 1: Load ClinVar and dbNSFP Hail Tables."""
        logger.info("🔄 [1/7] Loading Hail Tables...")

        self._clinvar_ht = hl.read_table(self.config.clinvar_ht)
        logger.info(f"   ✓ ClinVar table loaded: {self._clinvar_ht.count()} variants")

        self._dbnsfp_ht = hl.read_table(self.config.dbnsfp_ht)
        logger.info(f"   ✓ dbNSFP table loaded: {self._dbnsfp_ht.count()} variants")

    def _filter_clinvar(self) -> hl.Table:
        """Stage 2: Filter ClinVar to target genes or variants."""
        logger.info("🔄 [2/7] Filtering ClinVar to target variants...")

        if self._clinvar_ht is None:
            raise RuntimeError("ClinVar table not loaded. Run _load_tables first.")

        ht = self._clinvar_ht

        # Extract gene from GENEINFO field if available
        if "info" in ht.row and "GENEINFO" in ht.info:
            geneinfo_dtype = ht.info.GENEINFO.dtype
            if isinstance(geneinfo_dtype, hl.tarray):
                # VCF import with Number=. → array<str>: take first element
                geneinfo_str = hl.or_missing(
                    hl.is_defined(ht.info.GENEINFO) & (hl.len(ht.info.GENEINFO) > 0),
                    ht.info.GENEINFO[0],
                )
            else:
                # Plain string field
                geneinfo_str = ht.info.GENEINFO

            ht = ht.annotate(
                gene=hl.or_missing(
                    hl.is_defined(geneinfo_str) & (geneinfo_str != ""),
                    geneinfo_str.split(":")[0],
                )
            )

        # Filter based on input source
        if self.config.variants_path:
            # Load variant list and filter
            ht = self._filter_by_variant_list(ht)
        else:
            # Filter by gene set
            gene_set = self.config.get_gene_set()
            if gene_set:
                # Expand with HGNC aliases if configured
                if self.config.hgnc_path:
                    from hvantk.core.utils.gene_aliases import (
                        expand_gene_set_with_aliases,
                    )

                    pre_expand_count = len(gene_set)
                    gene_set, alias_map = expand_gene_set_with_aliases(
                        list(gene_set), self.config.hgnc_path
                    )
                    if alias_map:
                        logger.info(
                            f"   Expanded gene set with {len(alias_map)} "
                            f"aliases from HGNC "
                            f"({pre_expand_count} → {len(gene_set)} symbols)"
                        )
                        for alias, canonical in sorted(alias_map.items()):
                            logger.info(f"     {alias} → {canonical}")

                self._n_genes = (
                    pre_expand_count if self.config.hgnc_path else len(gene_set)
                )
                logger.info(f"   Filtering to {self._n_genes} genes")
                gene_literal = hl.literal(gene_set)
                ht = ht.filter(gene_literal.contains(ht.gene))

        # Filter by review stars if the field exists
        if "info" in ht.row and "CLNREVSTAT" in ht.info:
            if self.config.min_stars > 0:
                ht = self._filter_by_review_stars(ht)

        count = ht.count()
        logger.info(f"   ✓ Filtered ClinVar: {count} variants")

        self._labeled_ht = ht
        return ht

    def _filter_by_variant_list(self, ht: hl.Table) -> hl.Table:
        """Filter table by a list of variants from file.

        Variant file format: chr:pos:ref:alt (one per line)
        """
        variants_path = self.config.variants_path
        logger.info(f"   Loading variant list from {variants_path}")

        # Use the new helper
        variant_keys = parse_variant_list(variants_path)

        logger.info(f"   Loaded {len(variant_keys)} variants from file")

        if not variant_keys:
            raise ValueError(f"No valid variants found in {variants_path}")

        # Create a set of variant strings for filtering
        variant_set = {f"{v[0]}:{v[1]}:{v[2]}:{v[3]}" for v in variant_keys}
        variant_literal = hl.literal(variant_set)

        # Create variant key expression
        ht = ht.annotate(
            _variant_key=hl.delimit(
                [
                    hl.str(ht.locus.contig),
                    hl.str(ht.locus.position),
                    ht.alleles[0],
                    ht.alleles[1],
                ],
                ":",
            )
        )

        # Filter and drop
        ht = ht.filter(variant_literal.contains(ht._variant_key)).drop("_variant_key")

        logger.info(f"   ✓ Filtered to {ht.count()} variants")
        return ht

    def _filter_by_review_stars(self, ht: hl.Table) -> hl.Table:
        """Filter ClinVar variants by review status stars (CLNREVSTAT).

        CLNREVSTAT is imported from VCF as a string (or array<str>)
        describing the review status. This method maps it to an integer
        star count using CLNREVSTAT_STAR_MAP and filters accordingly.
        """
        min_stars = self.config.min_stars

        # Convert min_stars to integer if it's a string (e.g., from CLI)
        if isinstance(min_stars, str):
            try:
                min_stars = int(min_stars)
            except ValueError:
                raise ValueError(f"Invalid min_stars value: {min_stars}")

        logger.info(f"   Filtering ClinVar variants with review stars >= {min_stars}")

        star_map = hl.literal(CLNREVSTAT_STAR_MAP)

        # Normalize CLNREVSTAT: extract string from array if needed
        clnrevstat_dtype = ht.info.CLNREVSTAT.dtype
        if isinstance(clnrevstat_dtype, hl.tarray):
            clnrevstat_str = hl.or_missing(
                hl.is_defined(ht.info.CLNREVSTAT) & (hl.len(ht.info.CLNREVSTAT) > 0),
                ht.info.CLNREVSTAT[0],
            )
        else:
            clnrevstat_str = ht.info.CLNREVSTAT

        # Map string to star count, default 0 for unknown values
        stars_expr = star_map.get(clnrevstat_str, 0)

        return ht.filter(hl.is_defined(clnrevstat_str) & (stars_expr >= min_stars))

    def _assign_labels(self) -> hl.Table:
        """Stage 3: Assign binary labels to ClinVar variants."""
        logger.info("🔄 [3/7] Assigning binary labels (Pathogenic/Benign)...")

        if self._labeled_ht is None:
            raise RuntimeError(
                "ClinVar variants not filtered. Run _filter_clinvar first."
            )

        ht = self._labeled_ht

        # Normalize CLNSIG: extract string from array if needed
        clnsig_dtype = ht.info.CLNSIG.dtype
        if isinstance(clnsig_dtype, hl.tarray):
            clnsig = hl.or_missing(
                hl.is_defined(ht.info.CLNSIG) & (hl.len(ht.info.CLNSIG) > 0),
                ht.info.CLNSIG[0],
            )
        else:
            clnsig = ht.info.CLNSIG

        # Assign labels based on CLNSIG values
        ht = ht.annotate(
            label=hl.case()
            .when(
                hl.is_defined(clnsig)
                & (clnsig != "")
                & hl.literal(PATHOGENIC_LABELS).contains(clnsig),
                "Pathogenic",
            )
            .when(
                hl.is_defined(clnsig)
                & (clnsig != "")
                & hl.literal(BENIGN_LABELS).contains(clnsig),
                "Benign",
            )
            .default("Uncertain/Conflicting"),
        )

        # Filter out uncertain/conflicting variants by default
        if self.config.min_stars <= 0:
            ht = ht.filter(ht.label != "Uncertain/Conflicting")

        logger.info(
            f"   ✓ Assigned labels: {len(ht.aggregate(hl.agg.collect_as_set(ht.label)))} classes"
        )

        # Update instance variable so downstream stages see the labeled table
        self._labeled_ht = ht
        return ht

    def _annotate_scores(self) -> hl.Table:
        """Stage 4: Annotate ClinVar variants with dbNSFP scores."""
        logger.info("🔄 [4/7] Annotating with dbNSFP scores...")

        if self._labeled_ht is None:
            raise RuntimeError(
                "ClinVar variants not labeled. Run _assign_labels first."
            )

        if self._dbnsfp_ht is None:
            raise RuntimeError("dbNSFP table not loaded. Run _load_tables first.")

        ht = self._labeled_ht

        # Left-join annotation on shared (locus, alleles) key
        ht = ht.annotate(**self._dbnsfp_ht[ht.key])

        # Verify all requested score fields exist in the annotated table
        available_fields = set(ht.row)
        missing = [s for s in self.config.scores if s not in available_fields]

        if missing:
            raise RuntimeError(
                f"Requested scores not found in dbNSFP table: {missing}. "
                f"Available score-like fields: "
                f"{sorted(available_fields - {'locus', 'alleles'})}"
            )

        score_fields = list(self.config.scores)

        # Resolve dict-typed transcript scores to scalar, respecting directionality.
        # Higher-is-pathogenic scores use max(); lower-is-pathogenic use min().
        directions = dict(SCORE_DIRECTIONALITY)
        if self.config.score_directions:
            directions.update(self.config.score_directions)

        resolve_ann = {}
        for sf in score_fields:
            dtype = ht[sf].dtype
            if isinstance(dtype, hl.tdict):
                higher_is_pathogenic = directions.get(sf, True)
                agg_func = hl.max if higher_is_pathogenic else hl.min
                resolve_ann[sf] = hl.or_missing(
                    hl.is_defined(ht[sf]),
                    agg_func(ht[sf].values()),
                )
        if resolve_ann:
            agg_summary = {
                sf: ("max" if directions.get(sf, True) else "min") for sf in resolve_ann
            }
            logger.info(
                f"   Resolving {len(resolve_ann)} dict-typed scores to scalar: "
                f"{agg_summary}"
            )
            ht = ht.annotate(**resolve_ann)

        # Keep only label and score fields
        ht = ht.select("label", *score_fields)

        # Drop variants outside dbNSFP domain (all scores null)
        n_before = ht.count()
        score_defined = [hl.is_defined(ht[sf]) for sf in score_fields]
        ht = ht.filter(hl.any(lambda x: x, score_defined))
        n_after = ht.count()
        n_dropped = n_before - n_after

        logger.info(
            f"   Dropped {n_dropped} variants with no dbNSFP scores "
            f"(out of scope for score evaluation)"
        )
        logger.info(
            f"   ✓ Annotated with {len(score_fields)} dbNSFP scores: "
            f"{n_after} variants in score domain"
        )

        self._n_out_of_scope = n_dropped
        self._annotated_ht = ht
        return ht

    def _compute_missingness(self) -> Dict[str, ScoreMissingness]:
        """Stage 5: Compute missingness statistics for each score."""
        logger.info("🔄 [5/7] Computing missingness statistics...")

        if self._annotated_ht is None:
            raise RuntimeError("Variants not annotated. Run _annotate_scores first.")

        ht = self._annotated_ht

        # Identify score fields from config that are present in the table
        score_fields = [s for s in self.config.scores if s in ht.row]

        # Single materialization: select scores, convert to pandas
        df = ht.select(*score_fields).to_pandas()

        # Build numpy arrays for each score
        scores_np = {
            sf: df[sf].to_numpy(dtype=float, na_value=np.nan) for sf in score_fields
        }

        # Compute missingness using existing roc.py utility (single pass)
        missingness_results = compute_all_missingness(
            scores_np, max_missingness=self.config.max_missingness
        )

        scores_included, scores_excluded = filter_scores_by_missingness(
            missingness_results, max_missingness=self.config.max_missingness
        )

        logger.info(f"   ✓ Computed missingness for {len(missingness_results)} scores")
        logger.info(
            f"   ✓ Included: {len(scores_included)}, "
            f"Excluded: {len(scores_excluded)}"
        )

        self.state.outputs["missingness"] = missingness_results

        return missingness_results

    def _compute_roc_metrics(self) -> Dict[str, ROCResult]:
        """Stage 6: Compute ROC metrics for included scores."""
        logger.info("🔄 [6/7] Computing ROC metrics for included scores...")

        if self._annotated_ht is None:
            raise RuntimeError("Variants not annotated. Run _annotate_scores first.")

        ht = self._annotated_ht

        # Filter to pathogenic/benign only (exclude uncertain/conflicting)
        ht = ht.filter((ht.label == "Pathogenic") | (ht.label == "Benign"))

        n_labeled = ht.count()
        if n_labeled < self.config.min_variants:
            logger.warning(
                "   ⚠ Too few labeled variants (%d) for ROC analysis "
                "(minimum: %d). Skipping ROC computation.",
                n_labeled,
                self.config.min_variants,
            )
            self.state.outputs["roc_metrics"] = {}
            return {}

        # Get score fields that passed missingness threshold
        missingness = self.state.outputs.get("missingness", {})
        score_fields = [
            s
            for s in self.config.scores
            if s in ht.row and s in missingness and missingness[s].included_in_analysis
        ]

        if not score_fields:
            logger.warning("   ⚠ No scores passed the missingness threshold")
            return {}

        # Select only needed fields
        ht = ht.select("label", *score_fields)

        # Convert to pandas for easier NumPy extraction
        df = ht.to_pandas()

        # Convert labels to binary (1=Pathogenic, 0=Benign)
        labels = np.array([1 if label == "Pathogenic" else 0 for label in df["label"]])

        # Extract scores into dictionary of NumPy arrays
        scores_dict = {}
        for score_field in score_fields:
            scores_dict[score_field] = df[score_field].to_numpy(
                dtype=float, na_value=np.nan
            )

        # Compute ROC metrics for all scores at once
        roc_results = compute_roc_metrics(
            labels=labels,
            scores=scores_dict,
            max_missingness=self.config.max_missingness,
            pos_label=1,
            threshold_method=self.config.threshold_method,
            n_bootstrap=self.config.n_bootstrap,
        )

        if roc_results:
            logger.info(f"   ✓ Computed ROC metrics for {len(roc_results)} scores")
        else:
            logger.warning(
                "   ⚠ No scores produced ROC metrics after filtering to P/B variants"
            )

        self.state.outputs["roc_metrics"] = roc_results

        return roc_results

    def _generate_outputs(self) -> PSROCResult:
        """Stage 7: Generate output files and reports."""
        logger.info("🔄 [7/7] Generating output files and reports...")

        if self._annotated_ht is None:
            raise RuntimeError("Variants not annotated. Run _annotate_scores first.")

        ht = self._annotated_ht

        # Export annotated variants to TSV if requested
        if self.config.export_tsv:
            logger.info(
                f"   Exporting annotated variants to TSV: {self.paths['annotated_tsv']}"
            )
            ht.export(self.paths["annotated_tsv"])

        # Reduce excessive partitions inherited from source tables.
        # naive_coalesce merges adjacent partitions without a data shuffle,
        # unlike repartition() which would do a full re-balance.
        # n_partitions() is free (metadata only, no Hail action triggered).
        n_parts = ht.n_partitions()
        if n_parts > 10:
            ht = ht.naive_coalesce(10)

        # Write the annotated Hail Table to disk
        logger.info(f"   Writing annotated Hail Table: {self.paths['annotated_ht']}")
        ht.write(self.paths["annotated_ht"], overwrite=True)

        # Collect metrics and missingness results
        metrics = self.state.outputs.get("roc_metrics", {})
        missingness = self.state.outputs.get("missingness", {})

        result = PSROCResult(
            annotated_ht_path=self.paths["annotated_ht"],
            metrics=metrics,
            missingness=missingness,
            n_pathogenic=ht.filter(ht.label == "Pathogenic").count(),
            n_benign=ht.filter(ht.label == "Benign").count(),
            n_excluded=ht.filter(ht.label == "Uncertain/Conflicting").count(),
            n_total=ht.count(),
            n_out_of_scope=self._n_out_of_scope,
            n_genes=self._n_genes,
            scores_included=[
                s
                for s in self.config.scores
                if s in missingness and missingness[s].included_in_analysis
            ],
            scores_excluded=[
                s
                for s in self.config.scores
                if s in missingness and not missingness[s].included_in_analysis
            ],
            max_missingness_threshold=self.config.max_missingness,
            output_dir=self.config.output_dir,
        )

        # Save metrics and missingness to JSON
        logger.info(f"   Saving metrics to JSON: {self.paths['metrics_json']}")
        with open(self.paths["metrics_json"], "w") as f:
            json.dump({k: v.to_dict() for k, v in metrics.items()}, f, indent=2)

        logger.info(f"   Saving missingness to JSON: {self.paths['missingness_json']}")
        with open(self.paths["missingness_json"], "w") as f:
            json.dump({k: v.to_dict() for k, v in missingness.items()}, f, indent=2)

        # Generate and save plots
        if self.config.generate_plots:
            self._generate_plots(result)

        logger.info(
            f"   ✓ Output files and reports generated in: {self.config.output_dir}"
        )

        return result

    def _generate_plots(self, result: PSROCResult) -> None:
        """Generate visualization plots."""
        import matplotlib

        matplotlib.use("Agg")  # Non-interactive backend

        # Build title suffix from group name (e.g., " — Hereditary Cancer")
        name = self.config.group_name
        suffix = f" — {name}" if name else ""

        try:
            # ROC curves
            if result.metrics:
                base_roc_path = self.paths["roc_curves_png"].removesuffix(".png")
                plot_roc_curves(
                    result.metrics,
                    output_path=base_roc_path,
                    title=f"PSROC: ROC Curves ({len(result.metrics)} scores){suffix}",
                )
                logger.info(f"   ✓ ROC curves plot: {self.paths['roc_curves_png']}")

                # AUC comparison
                base_auc_path = self.paths["auc_comparison_png"].removesuffix(".png")
                plot_auc_comparison(
                    result.metrics,
                    output_path=base_auc_path,
                    title=f"AUC Comparison{suffix}",
                )
                logger.info(
                    f"   ✓ AUC comparison plot: {self.paths['auc_comparison_png']}"
                )

            # Missingness summary
            if result.missingness:
                base_missingness_path = self.paths["missingness_png"].removesuffix(
                    ".png"
                )
                plot_missingness_summary(
                    result.missingness,
                    output_path=base_missingness_path,
                    max_missingness_threshold=result.max_missingness_threshold,
                    title=f"Score Missingness{suffix}",
                )
                logger.info(f"   ✓ Missingness plot: {self.paths['missingness_png']}")

            # Summary dashboard
            if result.metrics and result.missingness:
                base_dashboard_path = self.paths["dashboard_png"].removesuffix(".png")
                plot_psroc_summary_dashboard(
                    result.metrics,
                    result.missingness,
                    output_path=base_dashboard_path,
                    max_missingness_threshold=result.max_missingness_threshold,
                    n_genes=result.n_genes,
                    n_pathogenic=result.n_pathogenic,
                    n_benign=result.n_benign,
                    title=f"PSROC Analysis Summary{suffix}",
                )
                logger.info(f"   ✓ Dashboard plot: {self.paths['dashboard_png']}")

        except Exception as e:
            logger.warning(f"   ⚠ Plot generation failed: {e}")
