"""
PSROC Pipeline Module

This module provides orchestration for PSROC (Prediction Score ROC Analysis)
workflows, from loading variant data through score annotation, ROC computation,
and report generation.

The pipeline evaluates the discriminative power of variant pathogenicity
prediction scores (e.g., CADD, REVEL, MetaLR) against ClinVar truth labels.

Example:
    >>> from hvantk.psroc.pipeline import PSROCConfig, PSROCPipeline
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
from typing import Optional, Dict, Any, List, Set, Union
from enum import Enum
import json
import logging
from datetime import datetime

import numpy as np
import hail as hl

from hvantk.psroc.roc import (
    ROCResult,
    ScoreMissingness,
    compute_roc_metrics,
    compute_all_missingness,
    filter_scores_by_missingness,
)
from hvantk.psroc.plots import (
    plot_roc_curves,
    plot_auc_comparison,
    plot_missingness_summary,
    plot_psroc_summary_dashboard,
)
from hvantk.utils.gene_sets import load_gene_set

logger = logging.getLogger(__name__)


# Label mapping constants (matching ClinvarDataStreamer in clinvar_streamer.py)
PATHOGENIC_LABELS = [
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic",
    "Pathogenic",
]
BENIGN_LABELS = [
    "Benign/Likely_benign",
    "Likely_benign",
    "Benign",
]


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
        export_tsv: Whether to export annotated variants as TSV.
        overwrite: Whether to overwrite existing output files.
        generate_plots: Whether to generate visualization plots.
        output_prefix: Prefix for output file names.
    """

    # Input sources (one of genes/genes_file/variants_path required)
    genes: Optional[List[str]] = None
    genes_file: Optional[str] = None
    variants_path: Optional[str] = None

    # Required table paths
    clinvar_ht: str = ""
    dbnsfp_ht: str = ""

    # Score configuration
    scores: List[str] = field(default_factory=list)

    # Output configuration
    output_dir: str = ""
    output_prefix: str = "psroc"

    # Processing options
    reference_genome: str = "GRCh38"
    min_stars: int = 1
    max_missingness: float = 0.3
    threshold_method: str = "youden"

    # Output options
    export_tsv: bool = False
    overwrite: bool = False
    generate_plots: bool = True

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
        ]
        if sum(sources) == 0:
            errors.append(
                "Must provide at least one of: --genes, --genes-file, or --variants"
            )
        if sum(sources) > 1:
            errors.append(
                "Cannot provide multiple variant sources. Use only one of: "
                "--genes, --genes-file, or --variants"
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
    outputs: Dict[str, str] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    start_time: Optional[str] = None
    end_time: Optional[str] = None

    def save(self, path: Path) -> None:
        """Save pipeline state to JSON file.

        Args:
            path: Path to save state file.
        """
        state_dict = {
            "config": asdict(self.config),
            "current_stage": self.current_stage.value if self.current_stage else None,
            "completed_stages": self.completed_stages,
            "outputs": self.outputs,
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
            f"Variants analyzed: {self.n_total}",
            f"  Pathogenic: {self.n_pathogenic}",
            f"  Benign: {self.n_benign}",
            f"  Excluded: {self.n_excluded}",
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
                "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
            )

        self.state = PSROCState(config=config)
        self._setup_output_paths()
        self._initialize_hail()

        # Internal state
        self._clinvar_ht: Optional[hl.Table] = None
        self._dbnsfp_ht: Optional[hl.Table] = None
        self._labeled_ht: Optional[hl.Table] = None
        self._annotated_ht: Optional[hl.Table] = None

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
        """Execute the PSROC pipeline.

        Returns:
            PSROCResult with computed metrics and output paths.

        Raises:
            Exception: If pipeline execution fails.
        """
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
            missingness = self._run_stage(PSROCStage.COMPUTE_MISSINGNESS)

            # Stage 6: Compute ROC metrics
            metrics = self._run_stage(PSROCStage.COMPUTE_ROC)

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
            ht = ht.annotate(
                gene=hl.if_else(
                    hl.is_defined(ht.info.GENEINFO)
                    & (hl.len(ht.info.GENEINFO) > 0),
                    ht.info.GENEINFO[0].split(":")[0],
                    hl.missing(hl.tstr),
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
                logger.info(f"   Filtering to {len(gene_set)} genes")
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

        # Read variant list
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

        logger.info(f"   Loaded {len(variant_keys)} variants from file")

        if not variant_keys:
            raise ValueError(f"No valid variants found in {variants_path}")

        # Create a set of variant strings for filtering
        variant_set = {
            f"{v[0]}:{v[1]}:{v[2]}:{v[3]}" for v in variant_keys
        }
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

        # Filter and drop temporary field
        ht = ht.filter(variant_literal.contains(ht._variant_key))
        ht = ht.drop("_variant_key")

        return ht

    def _filter_by_review_stars(self, ht: hl.Table) -> hl.Table:
        """Filter by ClinVar review status stars."""
        # Map CLNREVSTAT to star count
        # Based on ClinVar review status definitions
        star_mapping = {
            "practice_guideline": 4,
            "reviewed_by_expert_panel": 3,
            "criteria_provided,_multiple_submitters,_no_conflicts": 2,
            "criteria_provided,_conflicting_interpretations": 1,
            "criteria_provided,_single_submitter": 1,
            "no_assertion_criteria_provided": 0,
            "no_assertion_provided": 0,
        }

        star_literal = hl.literal(star_mapping)
        min_stars = self.config.min_stars

        # Get first review status and map to stars
        review_status = ht.info.CLNREVSTAT[0]
        normalized_status = review_status.lower().replace(" ", "_")
        stars = star_literal.get(normalized_status, 0)

        ht = ht.filter(stars >= min_stars)
        logger.info(f"   Filtered by review stars >= {min_stars}")

        return ht

    def _assign_labels(self) -> hl.Table:
        """Stage 3: Assign binary pathogenic/benign labels."""
        logger.info("🔄 [3/7] Assigning binary labels...")

        if self._labeled_ht is None:
            raise RuntimeError("Filtered table not available. Run _filter_clinvar first.")

        ht = self._labeled_ht

        # Check if CLNSIG field exists
        if "info" not in ht.row or "CLNSIG" not in ht.info:
            raise ValueError(
                "ClinVar table missing info.CLNSIG field. "
                "Ensure the table was built from ClinVar VCF."
            )

        # Create pathogenic/benign label sets
        pathogenic_set = hl.set(PATHOGENIC_LABELS)
        benign_set = hl.set(BENIGN_LABELS)

        # Check if any element in CLNSIG matches pathogenic or benign labels
        is_pathogenic = ht.info.CLNSIG.any(lambda x: pathogenic_set.contains(x))
        is_benign = ht.info.CLNSIG.any(lambda x: benign_set.contains(x))

        # Assign binary label: 1 = pathogenic, 0 = benign, missing = other
        ht = ht.annotate(
            psroc_label=hl.case()
            .when(is_pathogenic, 1)
            .when(is_benign, 0)
            .or_missing()
        )

        # Filter to variants with valid labels (exclude VUS, conflicting, etc.)
        labeled_ht = ht.filter(hl.is_defined(ht.psroc_label))

        # Count labels
        n_total = ht.count()
        label_counts = labeled_ht.aggregate(
            hl.struct(
                n_pathogenic=hl.agg.count_where(labeled_ht.psroc_label == 1),
                n_benign=hl.agg.count_where(labeled_ht.psroc_label == 0),
            )
        )
        n_labeled = labeled_ht.count()
        n_excluded = n_total - n_labeled

        self._n_pathogenic = label_counts.n_pathogenic
        self._n_benign = label_counts.n_benign
        self._n_excluded = n_excluded
        self._n_total = n_total

        logger.info(f"   ✓ Labels assigned:")
        logger.info(f"     Pathogenic: {self._n_pathogenic}")
        logger.info(f"     Benign: {self._n_benign}")
        logger.info(f"     Excluded (VUS/conflicting): {n_excluded}")

        # Warn if class imbalance is severe
        if self._n_pathogenic < 30:
            logger.warning(f"   ⚠ Low pathogenic count ({self._n_pathogenic}) may affect ROC reliability")
        if self._n_benign < 30:
            logger.warning(f"   ⚠ Low benign count ({self._n_benign}) may affect ROC reliability")

        self._labeled_ht = labeled_ht
        return labeled_ht

    def _annotate_scores(self) -> hl.Table:
        """Stage 4: Join with dbNSFP scores."""
        logger.info("🔄 [4/7] Annotating with dbNSFP scores...")

        if self._labeled_ht is None:
            raise RuntimeError("Labeled table not available. Run _assign_labels first.")
        if self._dbnsfp_ht is None:
            raise RuntimeError("dbNSFP table not loaded. Run _load_tables first.")

        ht = self._labeled_ht
        dbnsfp = self._dbnsfp_ht

        # Check which requested scores exist in dbNSFP
        dbnsfp_fields = set(dbnsfp.row)
        available_scores = []
        missing_scores = []

        for score in self.config.scores:
            if score in dbnsfp_fields:
                available_scores.append(score)
            else:
                missing_scores.append(score)
                logger.warning(f"   ⚠ Score '{score}' not found in dbNSFP table")

        if not available_scores:
            raise ValueError(
                f"None of the requested scores found in dbNSFP table: {self.config.scores}"
            )

        if missing_scores:
            logger.warning(
                f"   Missing scores will be skipped: {missing_scores}"
            )

        # Select only the scores we need from dbNSFP
        dbnsfp_subset = dbnsfp.select(*available_scores)

        # Join tables on variant key
        ht = ht.annotate(**dbnsfp_subset[ht.key])

        # Handle transcript-specific scores (stored as dicts in dbNSFP)
        # For dict scores, we take the max value across transcripts
        for score in available_scores:
            score_type = ht[score].dtype
            if isinstance(score_type, hl.tdict):
                # Take max value across transcripts
                ht = ht.annotate(
                    **{
                        score: hl.if_else(
                            hl.is_defined(ht[score]) & (hl.len(ht[score]) > 0),
                            hl.max(ht[score].values()),
                            hl.missing(hl.tfloat64),
                        )
                    }
                )

        # Checkpoint the annotated table
        ht = ht.checkpoint(
            self.paths["annotated_ht"],
            overwrite=self.config.overwrite,
        )

        logger.info(f"   ✓ Annotated table saved: {self.paths['annotated_ht']}")

        # Export TSV if requested
        if self.config.export_tsv:
            export_fields = ["psroc_label"] + available_scores
            if "gene" in ht.row:
                export_fields = ["gene"] + export_fields

            ht.select(*[f for f in export_fields if f in ht.row]).export(
                self.paths["annotated_tsv"]
            )
            logger.info(f"   ✓ TSV exported: {self.paths['annotated_tsv']}")

        self._annotated_ht = ht
        self._available_scores = available_scores
        return ht

    def _compute_missingness(self) -> Dict[str, ScoreMissingness]:
        """Stage 5: Compute missingness statistics for each score."""
        logger.info("🔄 [5/7] Computing missingness statistics...")

        if self._annotated_ht is None:
            raise RuntimeError("Annotated table not available. Run _annotate_scores first.")

        ht = self._annotated_ht

        # Collect data for missingness computation
        # We need labels and all score values
        data = ht.select(
            "psroc_label",
            *self._available_scores,
        ).collect()

        # Convert to numpy arrays
        labels = np.array([row.psroc_label for row in data])
        scores_dict = {
            score: np.array(
                [
                    float(row[score]) if row[score] is not None else np.nan
                    for row in data
                ],
                dtype=np.float64,
            )
            for score in self._available_scores
        }

        # Compute missingness for all scores
        missingness = compute_all_missingness(
            scores_dict, max_missingness=self.config.max_missingness
        )

        # Filter scores by missingness
        scores_included, scores_excluded = filter_scores_by_missingness(
            missingness, max_missingness=self.config.max_missingness
        )

        logger.info(f"   ✓ Missingness computed:")
        for score_name, stats in missingness.items():
            status = "included" if stats.included_in_analysis else "EXCLUDED"
            logger.info(f"     {score_name}: {stats.missingness_rate:.1%} missing [{status}]")

        if not scores_included:
            raise ValueError(
                f"All scores exceeded missingness threshold ({self.config.max_missingness:.0%}). "
                "Consider increasing --max-missingness."
            )

        # Store for later stages
        self._labels = labels
        self._scores_dict = scores_dict
        self._missingness = missingness
        self._scores_included = scores_included
        self._scores_excluded = scores_excluded

        return missingness

    def _compute_roc_metrics(self) -> Dict[str, ROCResult]:
        """Stage 6: Compute ROC metrics for included scores."""
        logger.info("🔄 [6/7] Computing ROC metrics...")

        if self._labels is None or self._scores_dict is None:
            raise RuntimeError("Data not prepared. Run _compute_missingness first.")

        # Filter to only included scores
        included_scores = {
            name: values
            for name, values in self._scores_dict.items()
            if name in self._scores_included
        }

        # Compute ROC metrics
        metrics = compute_roc_metrics(
            labels=self._labels,
            scores=included_scores,
            max_missingness=self.config.max_missingness,
            threshold_method=self.config.threshold_method,
        )

        logger.info(f"   ✓ ROC metrics computed for {len(metrics)} scores:")
        for name, roc in sorted(metrics.items(), key=lambda x: x[1].auc, reverse=True):
            logger.info(f"     {name}: AUC={roc.auc:.3f}")

        self._metrics = metrics
        return metrics

    def _generate_outputs(self) -> PSROCResult:
        """Stage 7: Generate plots, metrics, and reports."""
        logger.info("🔄 [7/7] Generating outputs...")

        # Create result object
        result = PSROCResult(
            annotated_ht_path=self.paths["annotated_ht"],
            metrics=self._metrics,
            missingness=self._missingness,
            n_pathogenic=self._n_pathogenic,
            n_benign=self._n_benign,
            n_excluded=self._n_excluded,
            n_total=self._n_total,
            scores_included=self._scores_included,
            scores_excluded=self._scores_excluded,
            max_missingness_threshold=self.config.max_missingness,
            output_dir=self.config.output_dir,
        )

        # Save metrics JSON
        metrics_dict = {
            "total_variants": result.n_total,
            "n_pathogenic": result.n_pathogenic,
            "n_benign": result.n_benign,
            "n_excluded": result.n_excluded,
            "scores_included": result.scores_included,
            "scores_excluded": result.scores_excluded,
            "max_missingness_threshold": result.max_missingness_threshold,
            "metrics": {name: r.to_dict() for name, r in result.metrics.items()},
        }

        with open(self.paths["metrics_json"], "w") as f:
            json.dump(metrics_dict, f, indent=2)
        logger.info(f"   ✓ Metrics saved: {self.paths['metrics_json']}")

        # Save missingness report
        missingness_dict = {
            "total_variants": result.n_total,
            "max_missingness_threshold": result.max_missingness_threshold,
            "scores_included": result.scores_included,
            "scores_excluded": result.scores_excluded,
            "scores": {name: m.to_dict() for name, m in result.missingness.items()},
        }

        with open(self.paths["missingness_json"], "w") as f:
            json.dump(missingness_dict, f, indent=2)
        logger.info(f"   ✓ Missingness report saved: {self.paths['missingness_json']}")

        # Generate plots if enabled
        if self.config.generate_plots and self._metrics:
            self._generate_plots(result)

        # Print summary
        print(result.summary())

        return result

    def _generate_plots(self, result: PSROCResult) -> None:
        """Generate visualization plots."""
        import matplotlib
        matplotlib.use("Agg")  # Non-interactive backend

        try:
            # ROC curves
            if result.metrics:
                fig = plot_roc_curves(
                    result.metrics,
                    output_path=self.paths["roc_curves_png"],
                    title=f"PSROC: ROC Curves ({len(result.metrics)} scores)",
                )
                fig.savefig(self.paths["roc_curves_png"], dpi=150, bbox_inches="tight")
                matplotlib.pyplot.close(fig)
                logger.info(f"   ✓ ROC curves plot: {self.paths['roc_curves_png']}")

                # AUC comparison
                fig = plot_auc_comparison(
                    result.metrics,
                    output_path=self.paths["auc_comparison_png"],
                )
                fig.savefig(self.paths["auc_comparison_png"], dpi=150, bbox_inches="tight")
                matplotlib.pyplot.close(fig)
                logger.info(f"   ✓ AUC comparison plot: {self.paths['auc_comparison_png']}")

            # Missingness summary
            if result.missingness:
                fig = plot_missingness_summary(
                    result.missingness,
                    output_path=self.paths["missingness_png"],
                    max_missingness_threshold=result.max_missingness_threshold,
                )
                fig.savefig(self.paths["missingness_png"], dpi=150, bbox_inches="tight")
                matplotlib.pyplot.close(fig)
                logger.info(f"   ✓ Missingness plot: {self.paths['missingness_png']}")

            # Summary dashboard
            if result.metrics and result.missingness:
                fig = plot_psroc_summary_dashboard(
                    result.metrics,
                    result.missingness,
                    output_path=self.paths["dashboard_png"],
                    max_missingness_threshold=result.max_missingness_threshold,
                )
                fig.savefig(self.paths["dashboard_png"], dpi=150, bbox_inches="tight")
                matplotlib.pyplot.close(fig)
                logger.info(f"   ✓ Dashboard plot: {self.paths['dashboard_png']}")

        except Exception as e:
            logger.warning(f"   ⚠ Plot generation failed: {e}")
