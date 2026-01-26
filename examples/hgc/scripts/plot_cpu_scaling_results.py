#!/usr/bin/env python3
"""
HGC CPU Scaling Results Plotter

Generate publication-quality plots showing HGC workflow CPU scaling (strong scaling):
- Runtime vs CPU Count
- Speedup vs CPU Count (with ideal linear speedup reference)
- Parallel Efficiency vs CPU Count
- Per-step breakdown of runtime scaling

Usage:
    python plot_cpu_scaling_results.py --results-dir ./cpu_scaling_results
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["font.size"] = 11


def load_cpu_scaling_data(results_dir: Path) -> pd.DataFrame:
    """Load CPU scaling results from CSV."""
    csv_path = results_dir / "cpu_scaling_summary.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CPU scaling summary not found: {csv_path}")

    df = pd.read_csv(csv_path)
    print(f"Loaded CPU scaling data: {len(df)} runs")
    print(f"CPU counts: {df['num_cpus'].tolist()}")
    return df


def plot_runtime_vs_cpus(df: pd.DataFrame, output_path: Path):
    """
    Plot total runtime vs CPU count.
    Shows how runtime decreases with more CPUs.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    cpu_counts = df["num_cpus"].values
    total_times = df["total_sec"].values

    # Plot actual data
    ax.plot(
        cpu_counts,
        total_times,
        "o-",
        linewidth=2.5,
        markersize=12,
        color="#2c3e50",
        label="Measured Runtime",
        zorder=3,
    )

    # Add ideal scaling curve (inverse relationship)
    baseline_cpus = cpu_counts[0]
    baseline_time = total_times[0]
    ideal_times = baseline_time * baseline_cpus / cpu_counts
    ax.plot(
        cpu_counts,
        ideal_times,
        "--",
        linewidth=2,
        color="#27ae60",
        alpha=0.7,
        label="Ideal Linear Scaling",
        zorder=2,
    )

    ax.set_xlabel("Number of CPU Cores", fontsize=14, fontweight="bold")
    ax.set_ylabel("Total Runtime (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Runtime vs CPU Count (Strong Scaling)",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper right", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add value labels
    for cpu, time in zip(cpu_counts, total_times):
        ax.annotate(
            f"{time:.0f}s",
            xy=(cpu, time),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3),
        )

    # Add sample size annotation
    sample_size = df["sample_size"].iloc[0]
    ax.text(
        0.02,
        0.98,
        f"Fixed Cohort: {sample_size} samples",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved runtime vs CPUs plot: {output_path}")
    plt.close()


def plot_speedup(df: pd.DataFrame, output_path: Path):
    """
    Plot speedup vs CPU count.
    Speedup = T(baseline) / T(n_cpus)
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    cpu_counts = df["num_cpus"].values
    speedups = df["speedup"].values

    # Calculate ideal linear speedup
    baseline_cpus = cpu_counts[0]
    ideal_speedup = cpu_counts / baseline_cpus

    # Plot actual speedup
    ax.plot(
        cpu_counts,
        speedups,
        "o-",
        linewidth=2.5,
        markersize=12,
        color="#e74c3c",
        label="Measured Speedup",
        zorder=3,
    )

    # Plot ideal linear speedup
    ax.plot(
        cpu_counts,
        ideal_speedup,
        "--",
        linewidth=2,
        color="#27ae60",
        alpha=0.7,
        label="Ideal Linear Speedup",
        zorder=2,
    )

    ax.set_xlabel("Number of CPU Cores", fontsize=14, fontweight="bold")
    ax.set_ylabel("Speedup (×)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Speedup vs CPU Count", fontsize=16, fontweight="bold", pad=20
    )
    ax.legend(loc="upper left", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add value labels
    for cpu, speedup in zip(cpu_counts, speedups):
        ax.annotate(
            f"{speedup:.2f}×",
            xy=(cpu, speedup),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3),
        )

    # Add sample size annotation
    sample_size = df["sample_size"].iloc[0]
    ax.text(
        0.98,
        0.02,
        f"Fixed Cohort: {sample_size} samples\nBaseline: {baseline_cpus} CPUs",
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved speedup plot: {output_path}")
    plt.close()


def plot_efficiency(df: pd.DataFrame, output_path: Path):
    """
    Plot parallel efficiency vs CPU count.
    Efficiency = Speedup / (CPUs / baseline_CPUs) * 100%
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    cpu_counts = df["num_cpus"].values
    efficiencies = df["efficiency"].values

    # Plot efficiency
    ax.plot(
        cpu_counts,
        efficiencies,
        "o-",
        linewidth=2.5,
        markersize=12,
        color="#3498db",
        label="Parallel Efficiency",
        zorder=3,
    )

    # Add 100% efficiency reference line
    ax.axhline(
        y=100,
        color="#27ae60",
        linestyle="--",
        linewidth=2,
        alpha=0.7,
        label="Ideal Efficiency (100%)",
        zorder=2,
    )

    # Add 50% efficiency reference line
    ax.axhline(
        y=50,
        color="#f39c12",
        linestyle=":",
        linewidth=1.5,
        alpha=0.5,
        label="50% Efficiency",
        zorder=1,
    )

    ax.set_xlabel("Number of CPU Cores", fontsize=14, fontweight="bold")
    ax.set_ylabel("Parallel Efficiency (%)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Parallel Efficiency vs CPU Count",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper right", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 110)

    # Add value labels
    for cpu, eff in zip(cpu_counts, efficiencies):
        ax.annotate(
            f"{eff:.1f}%",
            xy=(cpu, eff),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3),
        )

    # Add sample size annotation
    sample_size = df["sample_size"].iloc[0]
    ax.text(
        0.02,
        0.98,
        f"Fixed Cohort: {sample_size} samples",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved efficiency plot: {output_path}")
    plt.close()


def plot_step_breakdown(df: pd.DataFrame, output_path: Path):
    """
    Plot runtime breakdown by workflow step vs CPU count.
    Shows which steps scale better with more CPUs.
    """
    fig, ax = plt.subplots(figsize=(14, 8))

    cpu_counts = df["num_cpus"].values
    steps = ["gvcf_combine_sec", "vds_to_mt_sec", "compute_qc_sec", "mt_to_vcf_sec"]
    step_labels = ["GVCF → VDS", "VDS → MT", "Compute QC", "MT → VCF"]
    colors = ["#3498db", "#2ecc71", "#f39c12", "#e74c3c"]
    markers = ["o", "s", "^", "d"]

    # Plot each step
    for step, label, color, marker in zip(steps, step_labels, colors, markers):
        times = df[step].values
        ax.plot(
            cpu_counts,
            times,
            marker=marker,
            linewidth=2.5,
            markersize=10,
            color=color,
            label=label,
            alpha=0.85,
        )

    ax.set_xlabel("Number of CPU Cores", fontsize=14, fontweight="bold")
    ax.set_ylabel("Runtime (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Step-by-Step Scaling", fontsize=16, fontweight="bold", pad=20
    )
    ax.legend(loc="upper right", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add sample size annotation
    sample_size = df["sample_size"].iloc[0]
    ax.text(
        0.02,
        0.98,
        f"Fixed Cohort: {sample_size} samples",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved step breakdown plot: {output_path}")
    plt.close()


def plot_comprehensive_scaling(df: pd.DataFrame, output_path: Path):
    """
    Create a comprehensive 2x2 panel figure with all key metrics.
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    cpu_counts = df["num_cpus"].values
    total_times = df["total_sec"].values
    speedups = df["speedup"].values
    efficiencies = df["efficiency"].values

    baseline_cpus = cpu_counts[0]
    baseline_time = total_times[0]
    ideal_times = baseline_time * baseline_cpus / cpu_counts
    ideal_speedup = cpu_counts / baseline_cpus

    # Panel 1: Runtime vs CPUs
    ax1.plot(
        cpu_counts,
        total_times,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#2c3e50",
        label="Measured",
        zorder=3,
    )
    ax1.plot(
        cpu_counts,
        ideal_times,
        "--",
        linewidth=2,
        color="#27ae60",
        alpha=0.7,
        label="Ideal",
        zorder=2,
    )
    ax1.set_xlabel("CPU Cores", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Runtime (seconds)", fontsize=12, fontweight="bold")
    ax1.set_title("A) Runtime vs CPU Count", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Speedup
    ax2.plot(
        cpu_counts,
        speedups,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#e74c3c",
        label="Measured",
        zorder=3,
    )
    ax2.plot(
        cpu_counts,
        ideal_speedup,
        "--",
        linewidth=2,
        color="#27ae60",
        alpha=0.7,
        label="Ideal",
        zorder=2,
    )
    ax2.set_xlabel("CPU Cores", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Speedup (×)", fontsize=12, fontweight="bold")
    ax2.set_title("B) Speedup vs CPU Count", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Efficiency
    ax3.plot(
        cpu_counts,
        efficiencies,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#3498db",
        zorder=3,
    )
    ax3.axhline(y=100, color="#27ae60", linestyle="--", linewidth=2, alpha=0.7)
    ax3.axhline(y=50, color="#f39c12", linestyle=":", linewidth=1.5, alpha=0.5)
    ax3.set_xlabel("CPU Cores", fontsize=12, fontweight="bold")
    ax3.set_ylabel("Efficiency (%)", fontsize=12, fontweight="bold")
    ax3.set_title("C) Parallel Efficiency", fontsize=13, fontweight="bold")
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 110)

    # Panel 4: Step breakdown
    steps = ["gvcf_combine_sec", "vds_to_mt_sec", "compute_qc_sec", "mt_to_vcf_sec"]
    step_labels = ["GVCF→VDS", "VDS→MT", "QC", "MT→VCF"]
    colors = ["#3498db", "#2ecc71", "#f39c12", "#e74c3c"]
    markers = ["o", "s", "^", "d"]

    for step, label, color, marker in zip(steps, step_labels, colors, markers):
        times = df[step].values
        ax4.plot(
            cpu_counts,
            times,
            marker=marker,
            linewidth=2,
            markersize=8,
            color=color,
            label=label,
            alpha=0.85,
        )

    ax4.set_xlabel("CPU Cores", fontsize=12, fontweight="bold")
    ax4.set_ylabel("Runtime (seconds)", fontsize=12, fontweight="bold")
    ax4.set_title("D) Per-Step Scaling", fontsize=13, fontweight="bold")
    ax4.legend(fontsize=9, ncol=2)
    ax4.grid(True, alpha=0.3)

    # Add overall title
    sample_size = df["sample_size"].iloc[0]
    fig.suptitle(
        f"HGC Workflow CPU Scaling Analysis (Fixed Cohort: {sample_size} samples)",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved comprehensive scaling plot: {output_path}")
    plt.close()


def generate_summary_table(df: pd.DataFrame, output_path: Path):
    """Generate a formatted summary table as text file."""
    with open(output_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("HGC CPU SCALING BENCHMARK SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        sample_size = df["sample_size"].iloc[0]
        f.write(f"Fixed Cohort Size: {sample_size} samples\n")
        f.write(
            f"CPU Counts Tested: {', '.join(map(str, df['num_cpus'].tolist()))}\n\n"
        )

        f.write("-" * 80 + "\n")
        f.write(
            f"{'CPUs':>6} {'Total(s)':>10} {'Speedup':>10} {'Efficiency':>12} "
            f"{'GVCF→VDS':>12} {'VDS→MT':>10} {'QC':>10} {'MT→VCF':>10}\n"
        )
        f.write("-" * 80 + "\n")

        for _, row in df.iterrows():
            f.write(
                f"{row['num_cpus']:>6.0f} "
                f"{row['total_sec']:>10.1f} "
                f"{row['speedup']:>9.2f}× "
                f"{row['efficiency']:>11.1f}% "
                f"{row['gvcf_combine_sec']:>12.1f} "
                f"{row['vds_to_mt_sec']:>10.1f} "
                f"{row['compute_qc_sec']:>10.1f} "
                f"{row['mt_to_vcf_sec']:>10.1f}\n"
            )

        f.write("-" * 80 + "\n\n")

        # Calculate and report key metrics
        max_speedup = df["speedup"].max()
        max_speedup_cpus = df.loc[df["speedup"].idxmax(), "num_cpus"]
        best_efficiency = df["efficiency"].max()
        best_eff_cpus = df.loc[df["efficiency"].idxmax(), "num_cpus"]

        f.write("KEY METRICS:\n")
        f.write(
            f"  Maximum Speedup: {max_speedup:.2f}× at {max_speedup_cpus:.0f} CPUs\n"
        )
        f.write(
            f"  Best Efficiency: {best_efficiency:.1f}% at {best_eff_cpus:.0f} CPUs\n"
        )

        # Calculate time saved
        baseline_time = df.iloc[0]["total_sec"]
        fastest_time = df["total_sec"].min()
        time_saved = baseline_time - fastest_time
        percent_saved = (time_saved / baseline_time) * 100

        f.write(
            f"  Time Saved (baseline vs fastest): {time_saved:.0f}s ({percent_saved:.1f}%)\n"
        )

        f.write("\n" + "=" * 80 + "\n")

    print(f"Saved summary table: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot HGC CPU scaling benchmark results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing CPU scaling results",
    )

    args = parser.parse_args()

    if not args.results_dir.exists():
        print(f"ERROR: Results directory not found: {args.results_dir}")
        sys.exit(1)

    print(f"Loading CPU scaling results from: {args.results_dir}")

    try:
        # Load data
        df = load_cpu_scaling_data(args.results_dir)

        # Sort by CPU count
        df = df.sort_values("num_cpus")

        # Generate all plots
        print("\nGenerating plots...")

        plot_runtime_vs_cpus(df, args.results_dir / "cpu_scaling_runtime.png")
        plot_speedup(df, args.results_dir / "cpu_scaling_speedup.png")
        plot_efficiency(df, args.results_dir / "cpu_scaling_efficiency.png")
        plot_step_breakdown(df, args.results_dir / "cpu_scaling_step_breakdown.png")
        plot_comprehensive_scaling(
            df, args.results_dir / "cpu_scaling_comprehensive.png"
        )

        # Generate summary table
        generate_summary_table(df, args.results_dir / "cpu_scaling_summary.txt")

        print("\n" + "=" * 80)
        print("✓ All plots and summaries generated successfully!")
        print("=" * 80)
        print(f"\nOutput files in {args.results_dir}:")
        print("  - cpu_scaling_runtime.png")
        print("  - cpu_scaling_speedup.png")
        print("  - cpu_scaling_efficiency.png")
        print("  - cpu_scaling_step_breakdown.png")
        print("  - cpu_scaling_comprehensive.png (4-panel summary)")
        print("  - cpu_scaling_summary.txt (text summary)")

    except Exception as e:
        print(f"ERROR: Failed to generate plots: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
