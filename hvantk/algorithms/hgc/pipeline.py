"""
HGC Pipeline Module

This module provides orchestration for end-to-end gVCF processing workflows,
from raw gVCF files through combination, quality control, and VCF export.

The pipeline supports:
- Sequential execution of multiple stages
- Stage skipping for flexible workflows
- State persistence for recovery
- Comprehensive validation and error handling
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any, List
from enum import Enum
import json
import logging
import tempfile
from datetime import datetime

# hail_context applies the NumPy ``np.bool`` compatibility shim before it
# imports Hail, so it must be imported before ``import hail`` below.
from hvantk.core.utils.hail_context import init_hail
import hail as hl

from hvantk.algorithms.hgc import (
    combine_gvcfs,
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
    compute_sample_qc,
    compute_variant_qc,
    filter_samples_by_qc,
    filter_variants_by_qc,
    QCMetrics,
)

logger = logging.getLogger(__name__)


def _shown_partitions(n: Optional[int], *, skipped: bool = False) -> str:
    """Render `n_partitions` for the run plan.

    Keyed on `is None`, not truthiness: 0 is a rejected value, and `or` would print it
    as "auto", telling the user a dry run's plan is fine for an invocation that aborts.

    `skipped` covers --skip-vds-to-mt: that stage is the only consumer, so with it off
    the existing MatrixTable keeps whatever layout it already has. Printing the
    requested number there would repeat #208's mistake in a new place -- affirming a
    setting that no stage will read.
    """
    if skipped:
        return "n/a (VDS -> MT stage skipped)"
    return "auto (VDS layout)" if n is None else str(n)


class PipelineStage(Enum):
    """Enumeration of pipeline stages."""

    COMBINE_GVCFS = "combine_gvcfs"
    VDS_TO_MT = "vds_to_mt"
    COMPUTE_SAMPLE_QC = "compute_sample_qc"
    COMPUTE_VARIANT_QC = "compute_variant_qc"
    EXPORT_PVCF = "export_pvcf"


@dataclass
class PipelineConfig:
    """Configuration for HGC pipeline execution."""

    # Required parameters
    input_dir: str
    output_dir: str

    # Optional processing configuration
    tmp_dir: Optional[str] = None
    reference_genome: str = "GRCh38"
    # Target partition count for the VDS -> MatrixTable stage, forwarded to
    # convert_vds_to_mt, which COALESCES the dense MatrixTable before writing it
    # (see #207/#208). Deliberately not applied at the read: that approach was measured
    # and fails -- see the step 3b comment in algorithms/hgc/converters.py before
    # re-attempting it. NOT the combiner's interval tuning below: that governs stage 1,
    # this governs stage 2 onward.
    n_partitions: Optional[int] = None
    overwrite: bool = False
    output_prefix: str = "cohort"

    # GVCF-combiner tuning (stage 1). Hail derives ONE PARTITION PER INTERVAL, so the
    # interval size caps how many cores can do useful work in the combine stage: with
    # Hail's 1.2 Mb genome default a single chromosome yields few partitions
    # (chr20 -> 54, chr1 -> 208) and any cores beyond that idle.
    import_interval_size: Optional[int] = None
    use_exome_default_intervals: bool = False
    gvcf_batch_size: Optional[int] = None
    branch_factor: Optional[int] = None

    # Stage skip flags
    skip_combine_gvcfs: bool = False
    skip_vds_to_mt: bool = False
    skip_compute_sample_qc: bool = False
    skip_compute_variant_qc: bool = False
    skip_export_pvcf: bool = False

    # Skip the biallelic audit + repair in the VDS -> MT stage. `hvantk hgc vds2mt` has exposed
    # this since forever; the pipeline hard-coded False, so users of the *recommended* entry point
    # could not opt out at all.
    skip_validation: bool = False

    # Path overrides (for resuming from intermediate stages)
    vds_path: Optional[str] = None
    mt_path: Optional[str] = None

    # QC configuration
    min_sample_call_rate: float = 0.85
    min_variant_call_rate: float = 0.85
    apply_qc_filters: bool = False

    # Output options
    keep_intermediates: bool = True
    generate_qc_report: bool = False

    def validate(self) -> List[str]:
        """
        Validate configuration and return list of errors.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        # Check required paths
        if not Path(self.input_dir).exists():
            errors.append(f"Input directory does not exist: {self.input_dir}")

        # Check skip flag dependencies
        if self.skip_combine_gvcfs and not self.vds_path:
            errors.append("--vds-path required when using --skip-combine-gvcfs")

        if self.skip_vds_to_mt and not self.mt_path:
            errors.append("--mt-path required when using --skip-vds-to-mt")

        # Check path overrides exist if provided
        if self.vds_path and not Path(self.vds_path).exists():
            errors.append(f"Specified VDS path does not exist: {self.vds_path}")

        if self.mt_path and not Path(self.mt_path).exists():
            errors.append(f"Specified MT path does not exist: {self.mt_path}")

        # Validate QC thresholds
        if not 0 <= self.min_sample_call_rate <= 1:
            errors.append("min_sample_call_rate must be between 0 and 1")

        if not 0 <= self.min_variant_call_rate <= 1:
            errors.append("min_variant_call_rate must be between 0 and 1")

        # Validate GVCF-combiner tuning. Hail rejects these too, but only once the
        # combiner is constructed - well after Hail init and gVCF validation - and the
        # error is then wrapped in generic "check your Spark/GVCFs" guidance.
        if self.import_interval_size is not None and self.import_interval_size < 1:
            errors.append("import_interval_size must be at least 1 bp")

        if self.import_interval_size is not None and self.use_exome_default_intervals:
            errors.append(
                "import_interval_size and use_exome_default_intervals are mutually exclusive"
            )

        if self.gvcf_batch_size is not None and self.gvcf_batch_size < 1:
            errors.append("gvcf_batch_size must be at least 1")

        if self.branch_factor is not None and self.branch_factor < 2:
            errors.append("branch_factor must be at least 2")

        # Checked here, not only in convert_vds_to_mt, for the same reason as the three
        # knobs above: this is the gate that runs before Hail init and before stage 1.
        # convert_vds_to_mt does reject < 1, but stage 2 is reached only after the gVCF
        # combine has run to completion -- hours on a real cohort -- so a value that was
        # knowably wrong before any work started would cost the whole combine first.
        if self.n_partitions is not None and self.n_partitions < 1:
            errors.append("n_partitions must be at least 1")

        return errors

    def combiner_kwargs(self) -> Dict[str, Any]:
        """Build the kwargs forwarded to ``combine_gvcfs`` for stage 1.

        Only keys the user actually set are included; an empty dict preserves Hail's
        default partitioning (``use_genome_default_intervals``, 1.2 Mb).
        """
        kwargs: Dict[str, Any] = {}
        if self.import_interval_size is not None:
            kwargs["import_interval_size"] = self.import_interval_size
        if self.use_exome_default_intervals:
            kwargs["use_exome_default_intervals"] = True
        if self.gvcf_batch_size is not None:
            kwargs["gvcf_batch_size"] = self.gvcf_batch_size
        if self.branch_factor is not None:
            kwargs["branch_factor"] = self.branch_factor
        return kwargs


@dataclass
class PipelineState:
    """Tracks the state of pipeline execution."""

    config: PipelineConfig
    current_stage: Optional[PipelineStage] = None
    completed_stages: List[str] = field(default_factory=list)
    outputs: Dict[str, str] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    start_time: Optional[str] = None
    end_time: Optional[str] = None

    def save(self, path: Path):
        """
        Save pipeline state to JSON file.

        Args:
            path: Path to save state file
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
    def load(cls, path: Path) -> "PipelineState":
        """
        Load pipeline state from JSON file.

        Args:
            path: Path to state file

        Returns:
            PipelineState object
        """
        with open(path, "r") as f:
            state_dict = json.load(f)

        # Reconstruct config
        config = PipelineConfig(**state_dict["config"])

        # Reconstruct state
        state = cls(
            config=config,
            current_stage=(
                PipelineStage(state_dict["current_stage"])
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
        self, stage: PipelineStage, output_path: Optional[str] = None
    ):
        """Mark a stage as completed."""
        stage_name = stage.value
        if stage_name not in self.completed_stages:
            self.completed_stages.append(stage_name)

        if output_path:
            self.outputs[stage_name] = output_path

        logger.info(f"Stage {stage_name} marked as complete")

    def is_stage_complete(self, stage: PipelineStage) -> bool:
        """Check if a stage has been completed."""
        return stage.value in self.completed_stages


class PipelineRunner:
    """Orchestrates the HGC pipeline workflow."""

    def __init__(self, config: PipelineConfig):
        """
        Initialize pipeline runner.

        Args:
            config: Pipeline configuration
        """
        self.config = config
        self.state = PipelineState(config=config)
        self.logger = logging.getLogger(__name__)
        self._setup_output_paths()
        self._initialize_hail()

    def _setup_output_paths(self):
        """Initialize output paths based on configuration."""
        output_dir = Path(self.config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix = self.config.output_prefix
        qc_dir = output_dir / "qc"
        qc_dir.mkdir(exist_ok=True)

        logs_dir = output_dir / "logs"
        logs_dir.mkdir(exist_ok=True)

        self.paths = {
            "vds": self.config.vds_path or str(output_dir / f"{prefix}.vds"),
            "mt": self.config.mt_path or str(output_dir / f"{prefix}.mt"),
            "sample_qc": str(qc_dir / f"{prefix}_sample_qc.ht"),
            "variant_qc": str(qc_dir / f"{prefix}_variant_qc.ht"),
            "sample_qc_csv": str(qc_dir / f"{prefix}_sample_qc.csv"),
            "variant_qc_csv": str(qc_dir / f"{prefix}_variant_qc.csv"),
            "pvcf": str(output_dir / f"{prefix}.vcf.bgz"),
            "qc_report": str(qc_dir / f"{prefix}_qc_report.html"),
            "state": str(output_dir / ".pipeline_state.json"),
            "log": str(
                logs_dir / f'pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            ),
        }

        self.logger.info(f"Output paths configured: {self.paths}")

    def _initialize_hail(self):
        """Initialize Hail if not already initialized.

        Delegates to the managed, idempotent, thread-safe ``init_hail`` so
        repeated pipeline invocations in one process do not trigger a second
        ``hl.init`` (which Hail rejects), and ``HVANTK_SKIP_HAIL_INIT`` is honored.
        """
        tmp_dir = self.config.tmp_dir or tempfile.gettempdir()
        init_hail(tmp_dir=tmp_dir, quiet=False)
        self.logger.info(f"Hail initialized with tmp_dir: {tmp_dir}")

    def show_plan(self):
        """Display the execution plan without running."""
        print("\n" + "=" * 70)
        print("HGC PIPELINE EXECUTION PLAN")
        print("=" * 70)

        print("\n📁 Input/Output:")
        print(f"  Input directory:  {self.config.input_dir}")
        print(f"  Output directory: {self.config.output_dir}")
        print(f"  Output prefix:    {self.config.output_prefix}")

        print("\n🔧 Configuration:")
        print(f"  Reference genome: {self.config.reference_genome}")
        # Reaches convert_vds_to_mt as of #208. Before that this line printed a setting
        # no stage read -- keep it truthful, and say which stage it governs.
        print(
            f"  Partitions (MT):  "
            f"{_shown_partitions(self.config.n_partitions, skipped=self.config.skip_vds_to_mt)}"
        )
        print(f"  Overwrite:        {self.config.overwrite}")
        print(
            f"  Combiner options: {self.config.combiner_kwargs() or 'defaults (genome 1.2 Mb intervals)'}"
        )

        print("\n📊 Stages to execute:")
        stage_num = 1

        if not self.config.skip_combine_gvcfs:
            print(f"  [{stage_num}] Combine gVCFs → VDS")
            print(f"      Output: {self.paths['vds']}")
            stage_num += 1
        else:
            print(
                f"  [SKIP] Combine gVCFs (using existing VDS: {self.config.vds_path})"
            )

        if not self.config.skip_vds_to_mt:
            print(f"  [{stage_num}] Convert VDS → MatrixTable")
            print(f"      Output: {self.paths['mt']}")
            stage_num += 1
        else:
            print(f"  [SKIP] VDS to MT (using existing MT: {self.config.mt_path})")

        if not self.config.skip_compute_sample_qc:
            print(f"  [{stage_num}] Compute Sample QC")
            print(f"      Output: {self.paths['sample_qc']}")
            stage_num += 1
        else:
            print(f"  [SKIP] Compute Sample QC")

        if not self.config.skip_compute_variant_qc:
            print(f"  [{stage_num}] Compute Variant QC")
            print(f"      Output: {self.paths['variant_qc']}")
            stage_num += 1
        else:
            print(f"  [SKIP] Compute Variant QC")

        if not self.config.skip_export_pvcf:
            print(f"  [{stage_num}] Export cohort VCF")
            print(f"      Output: {self.paths['pvcf']}")
            if self.config.apply_qc_filters:
                print(f"      Apply QC filters: Yes")
                print(
                    f"        Min sample call rate: {self.config.min_sample_call_rate}"
                )
                print(
                    f"        Min variant call rate: {self.config.min_variant_call_rate}"
                )
            stage_num += 1
        else:
            print(f"  [SKIP] Export pVCF")

        if self.config.generate_qc_report:
            print(f"  [{stage_num}] Generate QC Report")
            print(f"      Output: {self.paths['qc_report']}")

        print("\n" + "=" * 70 + "\n")

    def run(self) -> PipelineState:
        """
        Execute the pipeline.

        Returns:
            PipelineState object with execution results
        """
        self.state.start_time = datetime.now().isoformat()
        self.logger.info("Starting HGC pipeline execution")

        try:
            # Stage 1: Combine gVCFs
            if not self.config.skip_combine_gvcfs:
                self._run_stage(PipelineStage.COMBINE_GVCFS)
            else:
                self.logger.info("Skipping gVCF combination stage")
                self.state.outputs["vds"] = self.config.vds_path

            # Stage 2: VDS to MT
            if not self.config.skip_vds_to_mt:
                self._run_stage(PipelineStage.VDS_TO_MT)
            else:
                self.logger.info("Skipping VDS to MT conversion stage")
                self.state.outputs["mt"] = self.config.mt_path

            # Stage 3a: Compute Sample QC
            if not self.config.skip_compute_sample_qc:
                self._run_stage(PipelineStage.COMPUTE_SAMPLE_QC)
            else:
                self.logger.info("Skipping sample QC computation stage")

            # Stage 3b: Compute Variant QC
            if not self.config.skip_compute_variant_qc:
                self._run_stage(PipelineStage.COMPUTE_VARIANT_QC)
            else:
                self.logger.info("Skipping variant QC computation stage")

            # Stage 4: Export pVCF
            if not self.config.skip_export_pvcf:
                self._run_stage(PipelineStage.EXPORT_PVCF)
            else:
                self.logger.info("Skipping pVCF export stage")

            # Generate QC report if requested
            if self.config.generate_qc_report:
                self._generate_qc_report()

            self.state.end_time = datetime.now().isoformat()
            self.logger.info("Pipeline execution completed successfully")

        except Exception as e:
            self.state.errors.append(str(e))
            self.state.end_time = datetime.now().isoformat()
            self.logger.exception(f"Pipeline execution failed: {e}")
            raise

        finally:
            # Save final state
            self.state.save(Path(self.paths["state"]))

        return self.state

    def _run_stage(self, stage: PipelineStage) -> bool:
        """
        Execute a single pipeline stage.

        Args:
            stage: Pipeline stage to execute

        Returns:
            True if successful, False otherwise
        """
        self.state.current_stage = stage
        self.logger.info(f"Running stage: {stage.value}")

        try:
            if stage == PipelineStage.COMBINE_GVCFS:
                output = self._run_combine_gvcfs()
            elif stage == PipelineStage.VDS_TO_MT:
                output = self._run_vds_to_mt()
            elif stage == PipelineStage.COMPUTE_SAMPLE_QC:
                output = self._run_compute_sample_qc()
            elif stage == PipelineStage.COMPUTE_VARIANT_QC:
                output = self._run_compute_variant_qc()
            elif stage == PipelineStage.EXPORT_PVCF:
                output = self._run_export_pvcf()
            else:
                raise ValueError(f"Unknown stage: {stage}")

            self.state.mark_stage_complete(stage, output)
            return True

        except Exception as e:
            error_msg = f"Stage {stage.value} failed: {e}"
            self.state.errors.append(error_msg)
            self.logger.exception(error_msg)
            raise

    def _run_combine_gvcfs(self) -> str:
        """Stage 1: Combine gVCFs into VDS."""
        self.logger.info("🔄 [1/5] Combining gVCF files into VDS...")

        combiner_kwargs = self.config.combiner_kwargs()
        if combiner_kwargs:
            self.logger.info(f"   Combiner options: {combiner_kwargs}")

        combine_gvcfs(
            gvcf_dir=self.config.input_dir,
            vds_output_path=self.paths["vds"],
            tmp_path=self.config.tmp_dir or tempfile.mkdtemp(),
            save_path=f"{self.paths['vds']}.plan",
            vdses=[],
            kwargs=combiner_kwargs,
            reference_genome=self.config.reference_genome,
        )

        self.logger.info(f"   ✓ VDS created: {self.paths['vds']}")
        return self.paths["vds"]

    def _run_vds_to_mt(self) -> str:
        """Stage 2: Convert VDS to MatrixTable."""
        self.logger.info("🔄 [2/5] Converting VDS to MatrixTable...")

        vds_path = self.state.outputs.get("vds") or self.config.vds_path

        convert_vds_to_mt(
            vds_path=vds_path,
            output_path=self.paths["mt"],
            adjust_genotypes=True,
            skip_split_multi=False,
            skip_validation=self.config.skip_validation,
            skip_keying_by_cols=False,
            overwrite=self.config.overwrite,
            # #208: this is the consumer `n_partitions` never had. It was accepted from
            # the CLI and echoed back in the run plan while reaching no stage at all, so
            # the run plan was affirmatively telling the user a setting had taken effect
            # when it had not.
            n_partitions=self.config.n_partitions,
        )

        self.logger.info(f"   ✓ MatrixTable created: {self.paths['mt']}")
        return self.paths["mt"]

    def _run_compute_sample_qc(self) -> str:
        """Stage 3a: Compute sample QC metrics."""
        self.logger.info("🔄 [3/5] Computing sample QC metrics...")

        mt_path = self.state.outputs.get("mt") or self.config.mt_path
        mt = hl.read_matrix_table(mt_path)

        # Compute sample QC
        mt_qc = compute_sample_qc(mt)

        # Extract and save sample QC table
        sample_qc_ht = mt_qc.cols().select("sample_qc")
        sample_qc_ht.write(self.paths["sample_qc"], overwrite=self.config.overwrite)

        # Export to CSV for easy inspection
        sample_qc_ht.export(self.paths["sample_qc_csv"])

        self.logger.info(f"   ✓ Sample QC saved: {self.paths['sample_qc']}")
        self.logger.info(f"   ✓ Sample QC CSV: {self.paths['sample_qc_csv']}")

        return self.paths["sample_qc"]

    def _run_compute_variant_qc(self) -> str:
        """Stage 3b: Compute variant QC metrics."""
        self.logger.info("🔄 [4/5] Computing variant QC metrics...")

        mt_path = self.state.outputs.get("mt") or self.config.mt_path
        mt = hl.read_matrix_table(mt_path)

        # Compute variant QC
        mt_qc = compute_variant_qc(mt)

        # Extract and save variant QC table
        variant_qc_ht = mt_qc.rows().select("variant_qc")
        variant_qc_ht.write(self.paths["variant_qc"], overwrite=self.config.overwrite)

        # Export sample to CSV (limit rows for large datasets)
        variant_qc_sample = variant_qc_ht.head(10000)
        variant_qc_sample.export(self.paths["variant_qc_csv"])

        self.logger.info(f"   ✓ Variant QC saved: {self.paths['variant_qc']}")
        self.logger.info(
            f"   ✓ Variant QC CSV (sample): {self.paths['variant_qc_csv']}"
        )

        return self.paths["variant_qc"]

    def _run_export_pvcf(self) -> str:
        """Stage 5: Export to project VCF."""
        self.logger.info("🔄 [5/5] Exporting cohort VCF...")

        mt_path = self.state.outputs.get("mt") or self.config.mt_path

        # Apply QC filters if requested
        if self.config.apply_qc_filters:
            self.logger.info("   Applying QC filters...")

            # Load the MT
            mt = hl.read_matrix_table(mt_path)

            # Filter samples
            if not self.config.skip_compute_sample_qc:
                mt = filter_samples_by_qc(
                    mt, min_call_rate=self.config.min_sample_call_rate
                )
                self.logger.info(
                    f"   ✓ Filtered samples by call rate >= {self.config.min_sample_call_rate}"
                )

            # Filter variants
            if not self.config.skip_compute_variant_qc:
                mt = filter_variants_by_qc(
                    mt, min_call_rate=self.config.min_variant_call_rate
                )
                self.logger.info(
                    f"   ✓ Filtered variants by call rate >= {self.config.min_variant_call_rate}"
                )

            # Save filtered MT
            filtered_mt_path = str(
                Path(self.paths["mt"]).parent
                / f"{self.config.output_prefix}_filtered.mt"
            )
            mt.write(filtered_mt_path, overwrite=self.config.overwrite)
            self.logger.info(f"   ✓ Filtered MT saved: {filtered_mt_path}")

            # Use the filtered MT for export
            mt_path_to_export = filtered_mt_path
        else:
            # Use the original MT
            mt_path_to_export = mt_path

        # Export to VCF
        convert_mt_to_multi_sample_vcf(
            mt_path=mt_path_to_export,
            vcf_path=self.paths["pvcf"],
            filter_adj_genotypes=True,
            min_ac=1,
            split_multi=True,
        )

        self.logger.info(f"   ✓ pVCF exported: {self.paths['pvcf']}")
        return self.paths["pvcf"]

    def _generate_qc_report(self):
        """Generate HTML QC report."""
        self.logger.info("📝 Generating HTML QC report...")

        mt_path = self.state.outputs.get("mt") or self.config.mt_path
        mt = hl.read_matrix_table(mt_path)

        # Check for QC annotations
        has_sample_qc = "sample_qc" in mt.col
        has_variant_qc = "variant_qc" in mt.row

        if not has_sample_qc and not has_variant_qc:
            self.logger.warning("No QC annotations found, skipping report generation")
            return

        # Create QCMetrics object
        sample_qc = mt.cols().select("sample_qc") if has_sample_qc else None
        variant_qc = mt.rows().select("variant_qc") if has_variant_qc else None

        qc_results = QCMetrics(mt=mt, sample_qc=sample_qc, variant_qc=variant_qc)

        # Generate report
        report_path = qc_results.generate_html_report(
            output_path=self.paths["qc_report"],
            title=f"QC Report - {self.config.output_prefix}",
        )

        self.logger.info(f"   ✓ QC report generated: {report_path}")
