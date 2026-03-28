"""
Stage-based pipeline orchestration for QTL cascade analysis.

Follows the same patterns as :mod:`hvantk.psroc.pipeline` and
:mod:`hvantk.enrichex.pipeline`:

* Configuration dataclass with ``validate()``
* Enumerated pipeline stages with per-stage error handling
* ``run()`` for single-tissue execution
* ``run_collection()`` for multi-tissue iteration
"""

import logging
import re
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from hvantk.qtlcascade.constants import (
    DEFAULT_COLOC_H4_THRESHOLD,
    DEFAULT_COLOC_P1,
    DEFAULT_COLOC_P2,
    DEFAULT_COLOC_P12,
    DEFAULT_COLOC_W,
    DEFAULT_COLOC_WINDOW_KB,
    DEFAULT_EQTL_P_THRESHOLD,
    DEFAULT_PQTL_P_THRESHOLD,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class CascadeConfig:
    """Configuration for the QTL cascade pipeline."""

    # Required paths
    eqtl_ht: str = ""
    pqtl_ht: str = ""
    output_dir: str = ""

    # Optional allpairs for coloc
    eqtl_allpairs_ht: str = ""
    pqtl_allpairs_ht: str = ""

    # Optional overlays
    constraint_ht: str = ""
    disease_genes_ht: str = ""

    # Tissue configuration
    tissues: List[str] = field(default_factory=list)

    # Thresholds
    eqtl_p_threshold: float = DEFAULT_EQTL_P_THRESHOLD
    pqtl_p_threshold: float = DEFAULT_PQTL_P_THRESHOLD
    coloc_window_kb: int = DEFAULT_COLOC_WINDOW_KB
    coloc_p1: float = DEFAULT_COLOC_P1
    coloc_p2: float = DEFAULT_COLOC_P2
    coloc_p12: float = DEFAULT_COLOC_P12
    coloc_W: float = DEFAULT_COLOC_W

    # Output options
    generate_plots: bool = True
    generate_report: bool = True
    overwrite: bool = False

    def validate(self) -> List[str]:
        """Return a list of validation errors (empty if valid)."""
        errors: List[str] = []
        if not self.eqtl_ht:
            errors.append("eqtl_ht is required.")
        if not self.pqtl_ht:
            errors.append("pqtl_ht is required.")
        if not self.output_dir:
            errors.append("output_dir is required.")
        # Coloc needs both allpairs
        if bool(self.eqtl_allpairs_ht) != bool(self.pqtl_allpairs_ht):
            errors.append(
                "Both eqtl_allpairs_ht and pqtl_allpairs_ht are required "
                "for colocalization (or omit both to skip coloc)."
            )
        return errors


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------


class CascadeStage(Enum):
    """Pipeline stages executed in order."""
    BUILD_CASCADE = "build_cascade"
    BUILD_GENE_SUMMARY = "build_gene_summary"
    RUN_COLOC = "run_coloc"
    GENERATE_OUTPUTS = "generate_outputs"


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class CascadeResult:
    """Per-tissue cascade results."""
    tissue: str
    cascade_ht_path: str = ""
    gene_summary_ht_path: str = ""
    class_counts: Optional[dict] = None
    n_cascade_genes: int = 0
    coloc_df: Optional[pd.DataFrame] = None


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


class CascadePipeline:
    """QTL cascade pipeline with stage-based execution.

    Parameters
    ----------
    config : CascadeConfig
        Validated configuration.
    """

    def __init__(self, config: CascadeConfig) -> None:
        errors = config.validate()
        if errors:
            raise ValueError(
                "CascadeConfig validation failed:\n"
                + "\n".join(f"  - {e}" for e in errors)
            )
        self.config = config
        self._output_dir = Path(config.output_dir)
        self._per_tissue_dir = self._output_dir / "per_tissue"
        self._plots_dir = self._output_dir / "plots"
        self._output_dir.mkdir(parents=True, exist_ok=True)
        self._per_tissue_dir.mkdir(exist_ok=True)
        if config.generate_plots:
            self._plots_dir.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show_plan(self) -> str:
        """Display the execution plan (no Hail required)."""
        cfg = self.config
        tissues = cfg.tissues or ["(all — no tissue filter)"]
        has_coloc = bool(cfg.eqtl_allpairs_ht and cfg.pqtl_allpairs_ht)
        stages = list(CascadeStage)
        if not has_coloc:
            stages = [s for s in stages if s != CascadeStage.RUN_COLOC]

        lines = [
            "=" * 70,
            "QTL CASCADE PIPELINE — EXECUTION PLAN",
            "=" * 70,
            "",
            f"Tissues:          {', '.join(tissues)}",
            f"eQTL table:       {cfg.eqtl_ht}",
            f"pQTL table:       {cfg.pqtl_ht}",
            f"Coloc:            {'yes' if has_coloc else 'no'}",
            f"Constraint:       {cfg.constraint_ht or '(none)'}",
            f"Disease genes:    {cfg.disease_genes_ht or '(none)'}",
            f"Output directory: {cfg.output_dir}",
            "",
            "Stages:",
        ]
        for s in stages:
            lines.append(f"  {s.value}")
        lines.append("")
        lines.append(f"Total runs: {len(cfg.tissues) if cfg.tissues else 1}")
        lines.append("=" * 70)

        plan = "\n".join(lines)
        logger.info("\n%s", plan)
        return plan

    def run(self, tissue: Optional[str] = None) -> CascadeResult:
        """Execute the cascade pipeline for a single tissue."""
        t0 = time.time()
        tissue_label = tissue or "all"
        logger.info("=" * 70)
        logger.info("QTL CASCADE — START (tissue: %s)", tissue_label)
        logger.info("=" * 70)

        safe_name = _sanitize(tissue_label)
        tissue_dir = str(self._per_tissue_dir / safe_name)

        result = CascadeResult(tissue=tissue_label)

        try:
            self._initialize_hail()

            # Stage 1: Build cascade
            result = self._stage_build_cascade(result, tissue, tissue_dir)

            # Stage 2: Gene summary
            result = self._stage_gene_summary(result, tissue_dir)

            # Stage 3: Coloc (optional)
            if self.config.eqtl_allpairs_ht and self.config.pqtl_allpairs_ht:
                result = self._stage_coloc(result, tissue)

            # Stage 4: Outputs
            result = self._stage_outputs(result, safe_name)

        except Exception:
            logger.exception("Pipeline failed for tissue %s", tissue_label)
            raise

        elapsed = time.time() - t0
        logger.info(
            "QTL CASCADE — COMPLETE (tissue: %s, %.1fs)", tissue_label, elapsed
        )
        return result

    def run_collection(self) -> Dict[str, CascadeResult]:
        """Execute the pipeline for each tissue in ``config.tissues``.

        Returns a dict mapping tissue name → :class:`CascadeResult`.
        """
        if not self.config.tissues:
            raise ValueError(
                "config.tissues is empty. Use run() for a single execution."
            )

        t0 = time.time()
        results: Dict[str, CascadeResult] = {}
        failed: List[str] = []

        for i, tissue in enumerate(sorted(self.config.tissues), 1):
            logger.info(
                "[%d/%d] Running cascade for tissue: %s",
                i, len(self.config.tissues), tissue,
            )
            try:
                results[tissue] = self.run(tissue=tissue)
            except Exception as exc:
                logger.warning("Tissue %s failed: %s", tissue, exc)
                failed.append(tissue)

        if not results:
            raise RuntimeError(
                f"All tissues failed: {', '.join(failed)}"
            )
        if failed:
            logger.warning(
                "%d tissue(s) failed: %s", len(failed), ", ".join(failed)
            )

        # Cross-tissue outputs
        if self.config.generate_plots and len(results) >= 2:
            self._generate_cross_tissue_outputs(results)

        if self.config.generate_report:
            self._generate_collection_report(results)

        elapsed = time.time() - t0
        logger.info(
            "Collection complete: %d/%d tissues (%.1fs)",
            len(results), len(self.config.tissues), elapsed,
        )
        return results

    # ------------------------------------------------------------------
    # Stage implementations
    # ------------------------------------------------------------------

    def _stage_build_cascade(self, result, tissue, tissue_dir):
        logger.info("Stage: %s", CascadeStage.BUILD_CASCADE.value)
        from hvantk.qtlcascade.cascade import build_cascade

        cascade_path = f"{tissue_dir}/cascade.ht"
        ht = build_cascade(
            eqtl_ht_path=self.config.eqtl_ht,
            pqtl_ht_path=self.config.pqtl_ht,
            output_path=cascade_path,
            eqtl_p_threshold=self.config.eqtl_p_threshold,
            pqtl_p_threshold=self.config.pqtl_p_threshold,
            tissue=tissue,
            overwrite=self.config.overwrite,
        )
        result.cascade_ht_path = cascade_path

        # Collect class counts
        import hail as hl
        counts = dict(
            ht.group_by(ht.cascade_class)
            .aggregate(n=hl.agg.count())
            .to_pandas()
            .set_index("cascade_class")["n"]
        )
        result.class_counts = counts
        logger.info("Cascade class counts: %s", counts)
        return result

    def _stage_gene_summary(self, result, tissue_dir):
        logger.info("Stage: %s", CascadeStage.BUILD_GENE_SUMMARY.value)
        from hvantk.qtlcascade.gene_summary import build_cascade_gene_summary

        gene_path = f"{tissue_dir}/gene_summary.ht"
        gene_ht = build_cascade_gene_summary(
            cascade_ht_path=result.cascade_ht_path,
            output_path=gene_path,
            constraint_ht_path=self.config.constraint_ht or None,
            disease_genes_ht_path=self.config.disease_genes_ht or None,
            coloc_df=result.coloc_df,
            overwrite=self.config.overwrite,
        )
        result.gene_summary_ht_path = gene_path
        result.n_cascade_genes = gene_ht.count()
        logger.info("Gene summary: %d genes", result.n_cascade_genes)

        # Export TSV for downstream use
        tsv_path = f"{tissue_dir}/gene_summary.tsv"
        gene_ht.export(tsv_path)
        logger.info("Gene summary exported: %s", tsv_path)
        return result

    def _stage_coloc(self, result, tissue):
        logger.info("Stage: %s", CascadeStage.RUN_COLOC.value)
        from hvantk.qtlcascade.coloc import run_coloc_per_gene
        import hail as hl

        # Get cascade gene IDs (genes with both eQTL and pQTL)
        ht = hl.read_table(result.cascade_ht_path)
        cascade_genes = list(
            ht.filter(ht.cascade_class == "eqtl_mediated")
            .aggregate(hl.agg.collect_as_set(ht.gene_id))
        )
        if not cascade_genes:
            logger.info("No eqtl_mediated genes — skipping coloc")
            return result

        logger.info("Running coloc for %d cascade genes", len(cascade_genes))
        coloc_df = run_coloc_per_gene(
            eqtl_allpairs_ht_path=self.config.eqtl_allpairs_ht,
            pqtl_allpairs_ht_path=self.config.pqtl_allpairs_ht,
            cascade_genes=cascade_genes,
            tissue=tissue,
            window_kb=self.config.coloc_window_kb,
            p1=self.config.coloc_p1,
            p2=self.config.coloc_p2,
            p12=self.config.coloc_p12,
            W=self.config.coloc_W,
        )
        result.coloc_df = coloc_df

        n_pass = (coloc_df["H4"] > DEFAULT_COLOC_H4_THRESHOLD).sum()
        logger.info(
            "Coloc: %d/%d genes with P(H4) > %.1f",
            n_pass, len(coloc_df), DEFAULT_COLOC_H4_THRESHOLD,
        )

        # Save coloc results
        safe = _sanitize(tissue or "all")
        coloc_path = self._per_tissue_dir / safe / "coloc_results.tsv"
        coloc_path.parent.mkdir(parents=True, exist_ok=True)
        coloc_df.to_csv(str(coloc_path), sep="\t", index=False)
        return result

    def _stage_outputs(self, result, safe_name):
        logger.info("Stage: %s", CascadeStage.GENERATE_OUTPUTS.value)

        if self.config.generate_plots:
            self._generate_tissue_plots(result, safe_name)

        return result

    # ------------------------------------------------------------------
    # Plot generation
    # ------------------------------------------------------------------

    def _generate_tissue_plots(self, result, safe_name):
        from hvantk.qtlcascade import plot as qtl_plot

        prefix = str(self._plots_dir / safe_name)

        if result.class_counts:
            try:
                qtl_plot.plot_cascade_classes(
                    result.class_counts,
                    output_path=f"{prefix}_cascade_classes.png",
                    title=f"Cascade Classes — {result.tissue}",
                )
            except Exception as exc:
                logger.warning("Cascade class plot failed: %s", exc)

        if result.coloc_df is not None and not result.coloc_df.empty:
            try:
                qtl_plot.plot_coloc_posteriors(
                    result.coloc_df,
                    output_path=f"{prefix}_coloc_posteriors.png",
                    title=f"Coloc P(H4) — {result.tissue}",
                )
            except Exception as exc:
                logger.warning("Coloc posterior plot failed: %s", exc)

    def _generate_cross_tissue_outputs(self, results):
        """Cross-tissue heatmap from collection results."""
        from hvantk.qtlcascade import plot as qtl_plot

        rows = []
        for tissue, res in results.items():
            if res.class_counts:
                for cls, cnt in res.class_counts.items():
                    rows.append({
                        "gene_id": cls,
                        "tissue": tissue,
                        "n_concordant": cnt,
                    })
        if not rows:
            return

        try:
            df = pd.DataFrame(rows)
            qtl_plot.plot_cross_tissue_heatmap(
                df,
                output_path=str(self._plots_dir / "cross_tissue_heatmap.png"),
                title="Cascade Counts Across Tissues",
            )
        except Exception as exc:
            logger.warning("Cross-tissue heatmap failed: %s", exc)

    def _generate_collection_report(self, results):
        from hvantk.qtlcascade.report import generate_report

        # Combine gene summaries
        import hail as hl
        gene_dfs = []
        for tissue, res in results.items():
            if res.gene_summary_ht_path:
                try:
                    df = hl.read_table(res.gene_summary_ht_path).to_pandas()
                    df["tissue"] = tissue
                    gene_dfs.append(df)
                except Exception:
                    pass
        gene_df = pd.concat(gene_dfs, ignore_index=True) if gene_dfs else None

        # Combine coloc
        coloc_dfs = [
            res.coloc_df for res in results.values()
            if res.coloc_df is not None and not res.coloc_df.empty
        ]
        coloc_df = pd.concat(coloc_dfs, ignore_index=True) if coloc_dfs else None

        # Aggregate class counts
        combined_counts: dict = {}
        for res in results.values():
            if res.class_counts:
                for cls, cnt in res.class_counts.items():
                    combined_counts[cls] = combined_counts.get(cls, 0) + cnt

        # Collect plot paths
        plot_paths = {}
        for p in self._plots_dir.glob("*.png"):
            plot_paths[p.stem] = str(p)

        generate_report(
            output_path=str(self._output_dir / "qtlcascade_report.html"),
            gene_summary_df=gene_df,
            coloc_df=coloc_df,
            class_counts=combined_counts,
            plot_paths=plot_paths,
            tissues=list(results.keys()),
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _initialize_hail(self):
        try:
            import hail as hl
            hl.current_backend()
            logger.info("Hail already initialised")
        except Exception:
            from hvantk.core.hail_context import init_hail
            init_hail()


def _sanitize(name: str) -> str:
    """Convert a tissue/group name to a filesystem-safe slug."""
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", name)
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "unknown"
