#!/usr/bin/env python
"""
Ancestry Inference Example: End-to-End Workflow

This example demonstrates how to use the ancestry inference pipeline
to predict genetic ancestry for query samples using synthetic data
that simulates the 1000 Genomes reference panel.

The example covers:
1. Generating synthetic query and reference MatrixTables
2. Running the ancestry inference pipeline via Python API
3. Interpreting results and viewing predictions
4. Generating visualizations (PCA, confusion matrix, etc.)
5. Creating an HTML report

Requirements:
    - hvantk installed with Hail: poetry install
    - Run from the repository root

Usage:
    python examples/ancestry/run_ancestry_example.py [--output-dir DIR]

Output:
    - predictions.ht: Hail Table with ancestry predictions
    - predictions.tsv: TSV export of predictions
    - pc_scores.ht: Hail Table with PC scores
    - rf_model.pkl: Trained Random Forest model
    - pca_loadings.ht: PCA loadings for projection
    - pipeline_stats.json: Pipeline execution statistics
    - ancestry_report.html: Comprehensive HTML report
    - plots/: Directory with visualization images
"""

import argparse
import sys
from pathlib import Path


def create_synthetic_data(data_dir: Path) -> tuple:
    """Generate synthetic MatrixTables for testing.

    Uses Hail's Balding-Nichols population genetics model to create
    realistic synthetic genotype data with known population structure.

    The approach:
    1. Generate all samples (reference + query) from the same simulation
    2. Split into reference (with ancestry labels) and query (labels hidden)
    3. This ensures query and reference share the same genetic variants

    The reference panel simulates 5 continental populations (EUR, AFR, EAS, SAS, AMR)
    similar to the 1000 Genomes superpopulations.

    Parameters
    ----------
    data_dir : Path
        Directory to save synthetic MatrixTables.

    Returns
    -------
    tuple
        (query_mt_path, reference_mt_path)
    """
    import hail as hl

    print("Generating synthetic data with known population structure...")

    query_path = str(data_dir / "synthetic_query.mt")
    reference_path = str(data_dir / "synthetic_reference.mt")

    # =========================================================================
    # Generate all samples from the same simulation
    # =========================================================================
    print("\n  Creating combined dataset (5 populations, 300 samples)...")

    # We'll generate 300 samples: 200 for reference, 100 for query
    # All from the same population model to ensure genetic compatibility
    total_samples = 300
    n_reference = 200
    n_query = 100

    # Simulate 5 continental populations with realistic Fst values
    # Higher Fst = more differentiation between populations
    combined_mt = hl.balding_nichols_model(
        n_populations=5,
        n_samples=total_samples,
        n_variants=10000,  # More variants for better PCA separation
        pop_dist=[0.2, 0.2, 0.2, 0.2, 0.2],  # Equal distribution
        fst=[0.12, 0.15, 0.10, 0.08, 0.11],  # Realistic Fst values
        n_partitions=8,
    )

    # Map population indices to superpopulation labels
    pop_labels = hl.literal(["EUR", "AFR", "EAS", "SAS", "AMR"])
    combined_mt = combined_mt.annotate_cols(
        ancestry=pop_labels[combined_mt.pop],
    )
    combined_mt = combined_mt.key_rows_by("locus", "alleles")

    # =========================================================================
    # Split into reference panel (first 200 samples with labels)
    # =========================================================================
    print(f"\n  Creating reference panel ({n_reference} samples with labels)...")

    # Filter to first 200 samples for reference
    ref_mt = combined_mt.filter_cols(combined_mt.sample_idx < n_reference)
    ref_mt = ref_mt.annotate_cols(
        s=hl.str("ref_") + hl.str(ref_mt.sample_idx),
    )
    ref_mt = ref_mt.key_cols_by("s")

    # Print reference statistics
    ref_counts = ref_mt.aggregate_cols(hl.agg.counter(ref_mt.ancestry))
    print(f"  Reference populations: {dict(sorted(ref_counts.items()))}")

    ref_mt.write(reference_path, overwrite=True)
    print(f"  Saved reference: {reference_path}")

    # =========================================================================
    # Split into query cohort (last 100 samples, labels hidden)
    # =========================================================================
    print(f"\n  Creating query cohort ({n_query} samples, labels hidden)...")

    # Filter to last 100 samples for query
    query_mt = combined_mt.filter_cols(combined_mt.sample_idx >= n_reference)

    # Renumber sample indices for query
    query_mt = query_mt.annotate_cols(
        s=hl.str("query_") + hl.str(query_mt.sample_idx - n_reference),
        # Store true population for validation (hidden in real scenarios)
        _true_pop=query_mt.ancestry,
    )
    # Drop the 'ancestry' column from query since it's "unknown" in real scenarios
    query_mt = query_mt.drop("ancestry")
    query_mt = query_mt.key_cols_by("s")

    # Print query statistics (hidden in real scenarios)
    true_counts = query_mt.aggregate_cols(hl.agg.counter(query_mt._true_pop))
    print(f"  True query populations (hidden): {dict(sorted(true_counts.items()))}")

    query_mt.write(query_path, overwrite=True)
    print(f"  Saved query: {query_path}")

    return query_path, reference_path


def main(output_dir: str = "./examples/ancestry/results") -> int:
    """Run the ancestry inference example with synthetic data."""

    print("=" * 70)
    print("Ancestry Inference Example: End-to-End Workflow")
    print("=" * 70)
    print()

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    data_dir = output_path / "data"
    data_dir.mkdir(exist_ok=True)

    results_dir = output_path
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print(f"Output directory: {output_path}")
    print()

    # =========================================================================
    # Step 1: Import hvantk modules
    # =========================================================================
    print("Loading hvantk modules...")
    try:
        from hvantk.core.hail_context import init_hail
        from hvantk.ancestry import (
            run_ancestry_inference,
            PipelineConfig,
            SUPERPOP_COLORS,
        )
    except ImportError as e:
        print(f"ERROR: Failed to import hvantk modules: {e}")
        print("Make sure hvantk is installed: poetry install")
        return 1

    # =========================================================================
    # Step 2: Initialize Hail
    # =========================================================================
    print("Initializing Hail...")
    init_hail(quiet=True)
    import hail as hl

    # =========================================================================
    # Step 3: Generate Synthetic Data
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 1: Generating Synthetic Data")
    print("-" * 70)

    query_path, reference_path = create_synthetic_data(data_dir)

    # =========================================================================
    # Step 4: Configure Pipeline
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 2: Configuring Ancestry Inference Pipeline")
    print("-" * 70)

    config = PipelineConfig(
        # Variant filtering (adjusted for synthetic data)
        min_af=0.05,
        max_af=0.95,
        min_call_rate=0.95,
        # Skip LD pruning for synthetic data (already independent)
        skip_ld_pruning=True,
        # PCA settings
        n_pcs=15,
        n_pcs_classify=10,
        # Classification
        n_estimators=100,
        min_prob=0.70,  # Lower threshold for more assignments
        random_seed=42,
        # Validation
        validate_model=True,
        n_cv_folds=5,
        min_samples_per_pop=10,
        # Merge settings
        min_shared_variants=100,
        # Checkpointing (optional, for large datasets)
        checkpoint_path=str(output_path / "checkpoints"),
        overwrite_checkpoints=True,
    )

    print(f"  Min AF: {config.min_af}")
    print(f"  Min call rate: {config.min_call_rate}")
    print(f"  Number of PCs: {config.n_pcs}")
    print(f"  PCs for classification: {config.n_pcs_classify}")
    print(f"  Minimum probability: {config.min_prob}")
    print(f"  Cross-validation folds: {config.n_cv_folds}")

    # =========================================================================
    # Step 5: Load MatrixTables
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 3: Loading MatrixTables")
    print("-" * 70)

    query_mt = hl.read_matrix_table(query_path)
    reference_mt = hl.read_matrix_table(reference_path)

    print(f"  Query samples: {query_mt.count_cols()}")
    print(f"  Query variants: {query_mt.count_rows()}")
    print(f"  Reference samples: {reference_mt.count_cols()}")
    print(f"  Reference variants: {reference_mt.count_rows()}")

    # Store hidden true labels for validation later
    true_labels_df = query_mt.cols().select("_true_pop").to_pandas()
    true_labels = true_labels_df.set_index("s")["_true_pop"].to_dict()

    # =========================================================================
    # Step 6: Run Ancestry Inference
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 4: Running Ancestry Inference Pipeline")
    print("-" * 70)

    result = run_ancestry_inference(
        query_mt=query_mt,
        reference_mt=reference_mt,
        ancestry_col="ancestry",
        config=config,
    )

    # =========================================================================
    # Step 7: Display Results
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 5: Results Summary")
    print("-" * 70)

    # Get predictions
    predictions_df = result.get_predictions_df()
    query_predictions = result.get_query_predictions()

    print(f"\n  Total samples analyzed: {len(predictions_df)}")
    print(f"  Query samples: {len(query_predictions)}")

    # Ancestry distribution
    print("\n  Predicted ancestry distribution (query samples):")
    for ancestry, count in (
        query_predictions["predicted_ancestry"].value_counts().items()
    ):
        pct = 100 * count / len(query_predictions)
        print(f"    {ancestry}: {count} ({pct:.1f}%)")

    # Cross-validation accuracy
    accuracy = result.get_accuracy()
    if accuracy:
        print(f"\n  Cross-validation accuracy: {accuracy:.2%}")

    # Variance explained
    var_explained = result.variance_explained()
    print(f"\n  Variance explained by PC1: {var_explained[0]:.2%}")
    print(f"  Variance explained by PC1-5: {sum(var_explained[:5]):.2%}")

    # Pipeline statistics
    stats = result.pipeline_stats
    print("\n  Pipeline statistics:")
    print(f"    Shared variants: {stats['n_shared_variants']}")
    print(
        f"    Variants after filtering: {stats.get('n_variants_after_filter', 'N/A')}"
    )
    print(f"    Training samples: {stats['n_training_samples']}")
    print(f"    Populations: {', '.join(stats['populations'])}")

    # =========================================================================
    # Step 8: Validate Against Hidden True Labels
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 6: Validation Against True Labels (Synthetic Data)")
    print("-" * 70)

    # Compare predictions to hidden true labels
    query_pred = query_predictions.copy()
    # The sample IDs are in the 's' column, not the index
    query_pred["true_ancestry"] = query_pred["s"].map(
        lambda sample_id: true_labels.get(sample_id, "unknown")
    )

    # Calculate inference accuracy (excluding unassigned)
    assigned = query_pred[query_pred["predicted_ancestry"] != "unassigned"]
    if len(assigned) > 0:
        correct = (assigned["predicted_ancestry"] == assigned["true_ancestry"]).sum()
        inference_accuracy = correct / len(assigned)
        print(f"\n  Inference accuracy (assigned samples): {inference_accuracy:.2%}")
        print(
            f"  Samples assigned: {len(assigned)}/{len(query_pred)} ({100*len(assigned)/len(query_pred):.1f}%)"
        )

        # Per-population accuracy
        print("\n  Per-population accuracy:")
        for pop in sorted(assigned["true_ancestry"].unique()):
            pop_mask = assigned["true_ancestry"] == pop
            pop_correct = (assigned.loc[pop_mask, "predicted_ancestry"] == pop).sum()
            pop_total = pop_mask.sum()
            if pop_total > 0:
                pop_acc = pop_correct / pop_total
                print(f"    {pop}: {pop_acc:.2%} ({pop_correct}/{pop_total})")
    else:
        print("  No samples were assigned ancestry (all below probability threshold)")

    # =========================================================================
    # Step 9: Save Results
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 7: Saving Results")
    print("-" * 70)

    # Use the built-in save method
    saved_paths = result.save(
        output_path=results_dir,
        save_model=True,
        save_loadings=True,
    )

    for name, path in saved_paths.items():
        print(f"  {name}: {path}")

    # =========================================================================
    # Step 10: Generate Visualizations
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 8: Generating Visualizations")
    print("-" * 70)

    # PCA scatter plot (PC1 vs PC2)
    fig = result.plot_pca(pc_x=1, pc_y=2, colors=SUPERPOP_COLORS)
    pca_path = plots_dir / "pca_scatter_pc1_pc2.png"
    fig.savefig(pca_path, dpi=150, bbox_inches="tight")
    print(f"  PCA scatter plot: {pca_path}")
    import matplotlib.pyplot as plt

    plt.close(fig)

    # PCA scatter plot (PC1 vs PC3)
    fig = result.plot_pca(pc_x=1, pc_y=3, colors=SUPERPOP_COLORS)
    pca_path = plots_dir / "pca_scatter_pc1_pc3.png"
    fig.savefig(pca_path, dpi=150, bbox_inches="tight")
    print(f"  PCA scatter plot: {pca_path}")
    plt.close(fig)

    # Two-panel PCA plot
    fig = result.plot_pca_panel()
    panel_path = plots_dir / "pca_panel.png"
    fig.savefig(panel_path, dpi=150, bbox_inches="tight")
    print(f"  PCA panel plot: {panel_path}")
    plt.close(fig)

    # Ancestry proportions
    fig = result.plot_ancestry_proportions()
    props_path = plots_dir / "ancestry_proportions.png"
    fig.savefig(props_path, dpi=150, bbox_inches="tight")
    print(f"  Ancestry proportions: {props_path}")
    plt.close(fig)

    # Probability distribution
    fig = result.plot_probability_distribution()
    prob_path = plots_dir / "probability_distribution.png"
    fig.savefig(prob_path, dpi=150, bbox_inches="tight")
    print(f"  Probability distribution: {prob_path}")
    plt.close(fig)

    # Variance explained
    fig = result.plot_variance_explained(n_pcs=config.n_pcs)
    var_path = plots_dir / "variance_explained.png"
    fig.savefig(var_path, dpi=150, bbox_inches="tight")
    print(f"  Variance explained: {var_path}")
    plt.close(fig)

    # Confusion matrix (if validation was performed)
    fig = result.plot_confusion_matrix(normalize=True)
    if fig:
        cm_path = plots_dir / "confusion_matrix.png"
        fig.savefig(cm_path, dpi=150, bbox_inches="tight")
        print(f"  Confusion matrix: {cm_path}")
        plt.close(fig)

    # =========================================================================
    # Step 11: Generate HTML Report
    # =========================================================================
    print()
    print("-" * 70)
    print("Step 9: Generating HTML Report")
    print("-" * 70)

    report_path = results_dir / "ancestry_report.html"
    result.generate_report(
        output_path=report_path,
        title="Ancestry Inference Example Report",
    )
    print(f"  Report: {report_path}")

    # =========================================================================
    # Summary
    # =========================================================================
    print()
    print("=" * 70)
    print("ANCESTRY INFERENCE EXAMPLE COMPLETE")
    print("=" * 70)
    print()
    print("Generated Output Files:")
    print(f"  {results_dir}/")
    print("  ├── data/")
    print("  │   ├── synthetic_query.mt/")
    print("  │   └── synthetic_reference.mt/")
    print("  ├── predictions.ht/")
    print("  ├── predictions.tsv")
    print("  ├── pc_scores.ht/")
    print("  ├── rf_model.pkl")
    print("  ├── pca_loadings.ht/")
    print("  ├── pipeline_stats.json")
    print("  ├── ancestry_report.html")
    print("  ├── plots/")
    print("  │   ├── pca_scatter_pc1_pc2.png")
    print("  │   ├── pca_scatter_pc1_pc3.png")
    print("  │   ├── pca_panel.png")
    print("  │   ├── ancestry_proportions.png")
    print("  │   ├── probability_distribution.png")
    print("  │   ├── variance_explained.png")
    print("  │   └── confusion_matrix.png")
    print("  └── checkpoints/")
    print()
    print(f"View the report: open {report_path}")
    print()

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run ancestry inference example with synthetic data"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./examples/ancestry/results",
        help="Output directory for results (default: ./examples/ancestry/results)",
    )
    args = parser.parse_args()

    sys.exit(main(output_dir=args.output_dir))
