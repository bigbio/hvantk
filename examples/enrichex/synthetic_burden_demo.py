#!/usr/bin/env python3
"""
End-to-end EnrichEx demo on a small synthetic dataset.

Exercises the full pipeline built across Phases 1–5.1:
  1. Generate synthetic cohort (Phase 5.1)
  2. Run BurdenPipeline with variant-class stratification (Phases 1–3)
  3. Generate plots: heatmap, volcano, forest (Phase 4)
  4. Check type-I error calibration on synonymous (negative control)
"""

import json
import logging
import sys
import tempfile
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s [%(levelname)s] %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger("enrichex_demo")


def main():
    import hail as hl
    from hvantk.core.utils.hail_context import init_hail

    # ── 0. Init Hail ──────────────────────────────────────────────────
    init_hail(app_name="enrichex-demo")

    from hvantk.algorithms.enrichex.burden import VariantFilter
    from hvantk.algorithms.enrichex.pipeline import BurdenConfig, BurdenPipeline
    from hvantk.algorithms.enrichex.plot import (
        plot_burden_volcano,
        plot_celltype_burden_heatmap,
        plot_celltype_forest,
    )
    from hvantk.algorithms.enrichex.simulation import (
        check_type_i_error,
        generate_synthetic_burden_cohort,
    )

    with tempfile.TemporaryDirectory(prefix="enrichex_demo_") as tmp:
        tmp_dir = Path(tmp)
        logger.info("Working directory: %s", tmp_dir)

        # ── 1. Generate synthetic data ────────────────────────────────
        logger.info("=" * 60)
        logger.info("STEP 1: Generating synthetic cohort")
        logger.info("=" * 60)

        # Two signal sets with different effect sizes
        signal_sets = {
            "cardiomyocytes": [f"GENE_{i}" for i in range(8)],
            "fibroblasts": [f"GENE_{i}" for i in range(15, 22)],
        }
        effect_sizes = {
            "cardiomyocytes": 3.0,  # strong signal
            "fibroblasts": 2.0,  # moderate signal
        }

        # Also add a null set (no planted signal) for calibration
        null_set = {
            "neurons": [f"GENE_{i}" for i in range(30, 38)],
        }
        all_gene_sets = {**signal_sets, **null_set}

        mt, pheno_ht, _ = generate_synthetic_burden_cohort(
            n_cases=200,
            n_controls=200,
            n_genes=50,
            variants_per_gene=8,
            signal_gene_sets=signal_sets,
            effect_sizes=effect_sizes,
            baseline_carrier_rate=0.06,
            seed=42,
        )

        # Write to disk for the pipeline
        mt_path = str(tmp_dir / "cohort.mt")
        pheno_path = str(tmp_dir / "pheno.ht")
        mt.write(mt_path, overwrite=True)
        pheno_ht.write(pheno_path, overwrite=True)

        n_rows, n_cols = mt.count()
        logger.info("Cohort: %d variants, %d samples", n_rows, n_cols)

        # Write gene set collections as JSON (one per "tissue")
        gs_dir = tmp_dir / "gene_sets"
        gs_dir.mkdir()
        for tissue_name, gs_dict in [
            (
                "heart",
                {
                    "cardiomyocytes": signal_sets["cardiomyocytes"],
                    "fibroblasts": signal_sets["fibroblasts"],
                },
            ),
            ("brain", {"neurons": null_set["neurons"]}),
        ]:
            coll_data = {
                "gene_sets": {
                    name: {"name": name, "genes": genes}
                    for name, genes in gs_dict.items()
                },
                "background_genes": [f"GENE_{i}" for i in range(50)],
            }
            gs_path = gs_dir / f"{tissue_name}.json"
            gs_path.write_text(json.dumps(coll_data, indent=2))
            logger.info("Wrote gene set collection: %s", gs_path)

        # ── 2. Run BurdenPipeline ─────────────────────────────────────
        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 2: Running BurdenPipeline")
        logger.info("=" * 60)

        output_dir = str(tmp_dir / "results")
        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            phenotype_field="is_case",
            phenotype_type="binary",
            gene_set_collections={
                "heart": str(gs_dir / "heart.json"),
                "brain": str(gs_dir / "brain.json"),
            },
            variant_classes={
                "lof": VariantFilter(
                    consequences=["stop_gained"],
                    max_af=1.0,  # no AF filtering on synthetic data
                    pass_only=False,
                    min_gq=0,
                    min_dp=0,
                    min_score=None,
                ),
                "missense": VariantFilter(
                    consequences=["missense_variant"],
                    max_af=1.0,
                    pass_only=False,
                    min_gq=0,
                    min_dp=0,
                    min_score=None,
                ),
                "synonymous": VariantFilter(
                    consequences=["synonymous_variant"],
                    max_af=1.0,
                    pass_only=False,
                    min_gq=0,
                    min_dp=0,
                    min_score=None,
                ),
            },
            min_carriers=0,
            correction_method="benjamini-hochberg",
            alpha=0.05,
            output_dir=output_dir,
            generate_report=False,
        )

        pipeline = BurdenPipeline(config)
        pipeline.show_plan()
        combined_df = pipeline.run()

        logger.info("")
        logger.info("Combined results shape: %s", combined_df.shape)
        if not combined_df.empty:
            logger.info("\nTop results by p-value:")
            cols = [
                "gene_set_name",
                "variant_class",
                "collection",
                "p_value",
                "p_adjusted",
                "significant",
            ]
            display_cols = [c for c in cols if c in combined_df.columns]
            print(
                combined_df.sort_values("p_value")[display_cols]
                .head(15)
                .to_string(index=False)
            )

        # ── 3. Generate plots ─────────────────────────────────────────
        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 3: Generating plots")
        logger.info("=" * 60)

        plots_dir = tmp_dir / "plots"
        plots_dir.mkdir()

        if not combined_df.empty:
            # 3a. Heatmap
            fig = plot_celltype_burden_heatmap(
                combined_df,
                str(plots_dir / "heatmap.png"),
                title="Synthetic Burden Heatmap",
            )
            import matplotlib.pyplot as plt

            plt.close(fig)
            logger.info("Wrote heatmap: %s", plots_dir / "heatmap.png")

            # 3b. Volcano
            if "beta" in combined_df.columns:
                effect = "beta"
            elif "odds_ratio" in combined_df.columns:
                effect = "odds_ratio"
            else:
                effect = None

            if effect:
                fig = plot_burden_volcano(
                    combined_df,
                    str(plots_dir / "volcano.png"),
                    effect_col=effect,
                    color_by="variant_class",
                    title="Synthetic Burden Volcano",
                )
                plt.close(fig)
                logger.info("Wrote volcano: %s", plots_dir / "volcano.png")

            # 3c. Forest for top cell type
            top_gs = combined_df.sort_values("p_value").iloc[0]["gene_set_name"]
            if "ci_lower" in combined_df.columns and effect:
                fig = plot_celltype_forest(
                    combined_df,
                    str(plots_dir / "forest.png"),
                    cell_type=top_gs,
                    effect_col=effect,
                    title=f"Forest: {top_gs}",
                )
                plt.close(fig)
                logger.info("Wrote forest: %s", plots_dir / "forest.png")
            else:
                logger.info("Skipping forest plot (missing ci_lower/ci_upper columns)")

        # ── 4. Type-I error check on synonymous ───────────────────────
        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 4: Synonymous calibration check")
        logger.info("=" * 60)

        if not combined_df.empty:
            syn = combined_df[combined_df["variant_class"] == "synonymous"]
            if not syn.empty:
                syn_pvals = syn["p_value"].dropna().tolist()
                cal = check_type_i_error(syn_pvals)
                logger.info("Synonymous p-values: %d tests", cal["n_tests"])
                logger.info(
                    "  Rejection rate: %.3f (expected: %.3f)",
                    cal["rejection_rate"],
                    cal["expected_rate"],
                )
                logger.info("  KS p-value:     %.4f", cal["ks_pvalue"])
                logger.info("  Calibrated:      %s", cal["calibrated"])
            else:
                logger.info("No synonymous results to check.")

            lof = combined_df[combined_df["variant_class"] == "lof"]
            if not lof.empty:
                logger.info("")
                logger.info(
                    "LoF results (expect signal in cardiomyocytes/fibroblasts):"
                )
                lof_display = lof.sort_values("p_value")[display_cols].head(5)
                print(lof_display.to_string(index=False))

        # ── 5. List output files ──────────────────────────────────────
        logger.info("")
        logger.info("=" * 60)
        logger.info("OUTPUT FILES")
        logger.info("=" * 60)
        for f in sorted(Path(output_dir).rglob("*")):
            if f.is_file():
                logger.info("  %s (%s)", f.relative_to(tmp_dir), f.stat().st_size)
        for f in sorted(plots_dir.rglob("*")):
            if f.is_file():
                logger.info("  %s (%s)", f.relative_to(tmp_dir), f.stat().st_size)

        logger.info("")
        logger.info("Demo complete.")


if __name__ == "__main__":
    main()
