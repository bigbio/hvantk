#!/usr/bin/env python3
"""
HGC Scalability Results Plotter

Generate publication-quality plots showing HGC workflow scalability:
- Runtime vs Sample Size (total and per-step breakdown)
- Memory usage vs Sample Size
- Scaling efficiency analysis

Usage:
    python plot_scalability_results.py --results-dir ./scalability_results
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["font.size"] = 11


def load_timing_data(results_dir: Path) -> pd.DataFrame:
    """Load timing results from CSV."""
    timing_csv = results_dir / "timings.csv"
    if not timing_csv.exists():
        raise FileNotFoundError(f"Timing CSV not found: {timing_csv}")

    df = pd.read_csv(timing_csv)
    print(f"Loaded timing data: {len(df)} runs")
    print(f"Sample sizes: {df['sample_size'].tolist()}")
    return df


def load_memory_data(results_dir: Path) -> pd.DataFrame:
    """Load memory usage results from CSV."""
    memory_csv = results_dir / "memory_usage.csv"
    if not memory_csv.exists():
        raise FileNotFoundError(f"Memory CSV not found: {memory_csv}")

    df = pd.read_csv(memory_csv)
    # Handle N/A values
    df["peak_memory_mb"] = pd.to_numeric(df["peak_memory_mb"], errors="coerce")
    print(f"Loaded memory data: {len(df)} runs")
    return df


def plot_runtime_breakdown(timing_df: pd.DataFrame, output_path: Path):
    """Create stacked bar chart showing runtime breakdown by workflow step."""
    fig, ax = plt.subplots(figsize=(14, 8))

    # Prepare data for stacked bars
    sample_sizes = timing_df["sample_size"]
    steps = ["gvcf_combine_sec", "vds_to_mt_sec", "compute_qc_sec", "mt_to_vcf_sec"]
    step_labels = ["GVCF → VDS", "VDS → MT", "Compute QC", "MT → VCF"]
    colors = ["#3498db", "#2ecc71", "#f39c12", "#e74c3c"]

    # Create stacked bars
    bottom = np.zeros(len(sample_sizes))

    for step, label, color in zip(steps, step_labels, colors):
        values = timing_df[step].values
        ax.bar(
            sample_sizes, values, bottom=bottom, label=label, color=color, alpha=0.85
        )
        bottom += values

    ax.set_xlabel("Sample Size", fontsize=14, fontweight="bold")
    ax.set_ylabel("Runtime (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Runtime Breakdown by Sample Size",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper left", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3, axis="y")

    # Add value labels on bars
    for i, (size, total) in enumerate(zip(sample_sizes, timing_df["total_sec"])):
        ax.text(
            size,
            total + total * 0.02,
            f"{total:.0f}s",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved runtime breakdown plot: {output_path}")
    plt.close()


def plot_total_runtime(timing_df: pd.DataFrame, output_path: Path):
    """Create line plot showing total runtime vs sample size with trend analysis."""
    fig, ax = plt.subplots(figsize=(12, 8))

    sample_sizes = timing_df["sample_size"].values
    total_times = timing_df["total_sec"].values

    # Plot actual data
    ax.plot(
        sample_sizes,
        total_times,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#2c3e50",
        label="Actual Runtime",
        zorder=3,
    )

    # Fit linear trend
    z_linear = np.polyfit(sample_sizes, total_times, 1)
    p_linear = np.poly1d(z_linear)
    ax.plot(
        sample_sizes,
        p_linear(sample_sizes),
        "--",
        linewidth=2,
        color="#e74c3c",
        alpha=0.7,
        label=f"Linear Fit: y = {z_linear[0]:.2f}x + {z_linear[1]:.0f}",
    )

    # Calculate R-squared
    residuals = total_times - p_linear(sample_sizes)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((total_times - np.mean(total_times)) ** 2)
    if ss_tot > 0:
        r_squared = 1 - (ss_res / ss_tot)
    else:
        r_squared = 1.0

    ax.set_xlabel("Sample Size", fontsize=14, fontweight="bold")
    ax.set_ylabel("Total Runtime (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Total Runtime vs Sample Size",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper left", fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add R-squared annotation
    ax.text(
        0.98,
        0.05,
        f"$R^2$ = {r_squared:.4f}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # Add value labels
    for size, time in zip(sample_sizes, total_times):
        ax.annotate(
            f"{time:.0f}s",
            xy=(size, time),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved total runtime plot: {output_path}")
    plt.close()


def plot_memory_usage(memory_df: pd.DataFrame, output_path: Path):
    """Create line plot showing memory usage vs sample size."""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Filter out N/A values
    valid_data = memory_df[memory_df["peak_memory_mb"].notna()]

    if len(valid_data) == 0:
        print("WARNING: No valid memory data available for plotting")
        return

    sample_sizes = valid_data["sample_size"].values
    memory_mb = valid_data["peak_memory_mb"].values
    memory_gb = memory_mb / 1024

    # Plot memory usage
    ax.plot(
        sample_sizes,
        memory_gb,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#9b59b6",
        label="Peak Memory Usage",
        zorder=3,
    )

    # Fit linear trend
    if len(sample_sizes) > 1:
        z = np.polyfit(sample_sizes, memory_gb, 1)
        p = np.poly1d(z)
        ax.plot(
            sample_sizes,
            p(sample_sizes),
            "--",
            linewidth=2,
            color="#e74c3c",
            alpha=0.7,
            label=f"Linear Fit: y = {z[0]:.4f}x + {z[1]:.2f}",
        )

    ax.set_xlabel("Sample Size", fontsize=14, fontweight="bold")
    ax.set_ylabel("Peak Memory Usage (GB)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Memory Usage vs Sample Size",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper left", fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add value labels
    for size, mem in zip(sample_sizes, memory_gb):
        ax.annotate(
            f"{mem:.1f}GB",
            xy=(size, mem),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved memory usage plot: {output_path}")
    plt.close()


def plot_step_comparison(timing_df: pd.DataFrame, output_path: Path):
    """Create line plot comparing runtime of each workflow step."""
    fig, ax = plt.subplots(figsize=(14, 8))

    sample_sizes = timing_df["sample_size"]
    steps = {
        "GVCF → VDS": ("gvcf_combine_sec", "#3498db"),
        "VDS → MT": ("vds_to_mt_sec", "#2ecc71"),
        "Compute QC": ("compute_qc_sec", "#f39c12"),
        "MT → VCF": ("mt_to_vcf_sec", "#e74c3c"),
    }

    for label, (col, color) in steps.items():
        ax.plot(
            sample_sizes,
            timing_df[col],
            "o-",
            linewidth=2,
            markersize=8,
            color=color,
            label=label,
            alpha=0.85,
        )

    ax.set_xlabel("Sample Size", fontsize=14, fontweight="bold")
    ax.set_ylabel("Runtime (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Step-by-Step Runtime Comparison",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper left", fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved step comparison plot: {output_path}")
    plt.close()


def plot_scaling_efficiency(timing_df: pd.DataFrame, output_path: Path):
    """Plot scaling efficiency (time per sample)."""
    fig, ax = plt.subplots(figsize=(12, 8))

    sample_sizes = timing_df["sample_size"].values
    total_times = timing_df["total_sec"].values

    # Validate that sample sizes are greater than zero
    if (sample_sizes == 0).any():
        raise ValueError("Sample sizes must be greater than zero")

    time_per_sample = total_times / sample_sizes

    ax.plot(
        sample_sizes,
        time_per_sample,
        "o-",
        linewidth=2.5,
        markersize=10,
        color="#16a085",
        label="Time per Sample",
        zorder=3,
    )

    ax.set_xlabel("Sample Size", fontsize=14, fontweight="bold")
    ax.set_ylabel("Time per Sample (seconds)", fontsize=14, fontweight="bold")
    ax.set_title(
        "HGC Workflow Scaling Efficiency", fontsize=16, fontweight="bold", pad=20
    )
    ax.legend(loc="best", fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)

    # Add value labels
    for size, tps in zip(sample_sizes, time_per_sample):
        ax.annotate(
            f"{tps:.2f}s",
            xy=(size, tps),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )

    # Add horizontal line for reference (ideal would be constant)
    mean_tps = np.mean(time_per_sample)
    ax.axhline(
        y=mean_tps,
        color="red",
        linestyle="--",
        alpha=0.5,
        label=f"Mean: {mean_tps:.2f}s/sample",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved scaling efficiency plot: {output_path}")
    plt.close()


def generate_summary_report(
    timing_df: pd.DataFrame, memory_df: pd.DataFrame, output_path: Path
):
    """Generate a text summary report with key statistics."""
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("HGC SCALABILITY BENCHMARK - SUMMARY REPORT")
    report_lines.append("=" * 80)
    report_lines.append("")

    # Runtime statistics
    report_lines.append("RUNTIME STATISTICS:")
    report_lines.append("-" * 80)
    for _, row in timing_df.iterrows():
        size = int(row["sample_size"])
        total = row["total_sec"]
        report_lines.append(
            f"Sample Size: {size:4d} | Total Runtime: {total:8.1f}s ({total/60:6.1f}min)"
        )

    report_lines.append("")
    report_lines.append("STEP BREAKDOWN (% of total time):")
    report_lines.append("-" * 80)
    steps = {
        "GVCF → VDS": "gvcf_combine_sec",
        "VDS → MT": "vds_to_mt_sec",
        "Compute QC": "compute_qc_sec",
        "MT → VCF": "mt_to_vcf_sec",
    }

    for step_name, col in steps.items():
        avg_pct = (timing_df[col] / timing_df["total_sec"] * 100).mean()
        report_lines.append(f"{step_name:15s}: {avg_pct:5.1f}% average")

    report_lines.append("")

    # Scaling analysis
    sample_sizes = timing_df["sample_size"].values
    total_times = timing_df["total_sec"].values

    report_lines.append("SCALING ANALYSIS:")
    report_lines.append("-" * 80)

    # Linear fit
    z = np.polyfit(sample_sizes, total_times, 1)
    p = np.poly1d(z)
    residuals = total_times - p(sample_sizes)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((total_times - np.mean(total_times)) ** 2)
    if ss_tot > 0:
        r_squared = 1 - (ss_res / ss_tot)
    else:
        r_squared = 1.0

    report_lines.append(f"Linear fit: y = {z[0]:.4f}x + {z[1]:.2f}")
    report_lines.append(f"R-squared: {r_squared:.4f}")
    report_lines.append(f"Interpretation: Adding 1 sample adds ~{z[0]:.2f} seconds")
    report_lines.append("")

    # Doubling analysis
    report_lines.append("DOUBLING ANALYSIS:")
    report_lines.append("-" * 80)
    for i in range(len(sample_sizes) - 1):
        size_ratio = sample_sizes[i + 1] / sample_sizes[i]
        time_ratio = total_times[i + 1] / total_times[i]
        report_lines.append(
            f"{sample_sizes[i]:4.0f} → {sample_sizes[i+1]:4.0f} samples "
            f"({size_ratio:4.2f}x): {time_ratio:4.2f}x time increase"
        )

    report_lines.append("")

    # Memory statistics
    valid_memory = memory_df[memory_df["peak_memory_mb"].notna()]
    if len(valid_memory) > 0:
        report_lines.append("MEMORY USAGE:")
        report_lines.append("-" * 80)
        for _, row in valid_memory.iterrows():
            size = int(row["sample_size"])
            mem_gb = row["peak_memory_mb"] / 1024
            mem_per_sample = mem_gb / size * 1024  # MB per sample
            report_lines.append(
                f"Sample Size: {size:4d} | Peak Memory: {mem_gb:6.1f}GB "
                f"({mem_per_sample:5.1f}MB/sample)"
            )
        report_lines.append("")

    report_lines.append("=" * 80)

    # Write report
    with open(output_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"Saved summary report: {output_path}")

    # Also print to console
    print("\n" + "\n".join(report_lines))


def main():
    parser = argparse.ArgumentParser(
        description="Generate scalability plots from HGC benchmark results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing benchmark results (timings.csv, memory_usage.csv)",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="hgc_scalability",
        help="Prefix for output plot files (default: hgc_scalability)",
    )

    args = parser.parse_args()

    if not args.results_dir.exists():
        print(f"ERROR: Results directory not found: {args.results_dir}")
        return 1

    print("=" * 80)
    print("HGC Scalability Results Plotter")
    print("=" * 80)
    print(f"Results directory: {args.results_dir}")
    print("")

    # Load data
    try:
        timing_df = load_timing_data(args.results_dir)
        memory_df = load_memory_data(args.results_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return 1

    # Create plots directory
    plots_dir = args.results_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print("")
    print("Generating plots...")
    print("-" * 80)

    # Generate all plots
    plot_runtime_breakdown(
        timing_df, plots_dir / f"{args.output_prefix}_runtime_breakdown.png"
    )
    plot_total_runtime(timing_df, plots_dir / f"{args.output_prefix}_total_runtime.png")
    plot_step_comparison(
        timing_df, plots_dir / f"{args.output_prefix}_step_comparison.png"
    )
    plot_scaling_efficiency(
        timing_df, plots_dir / f"{args.output_prefix}_scaling_efficiency.png"
    )

    if memory_df["peak_memory_mb"].notna().any():
        plot_memory_usage(
            memory_df, plots_dir / f"{args.output_prefix}_memory_usage.png"
        )

    # Generate summary report
    report_path = args.results_dir / f"{args.output_prefix}_report.txt"
    generate_summary_report(timing_df, memory_df, report_path)

    print("")
    print("=" * 80)
    print("All plots generated successfully!")
    print(f"Plots saved to: {plots_dir}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    exit(main())
