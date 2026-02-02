#!/usr/bin/env python3
"""
Test script to validate HGC scalability benchmark setup.

This script checks that all required dependencies are available
and the configuration is correct.
"""

import argparse
import os
import sys
from pathlib import Path


def check_import(module_name: str, package_name: str = None) -> bool:
    """Check if a module can be imported."""
    try:
        __import__(module_name)
        print(f"✓ {package_name or module_name} is available")
        return True
    except ImportError:
        print(f"✗ {package_name or module_name} is NOT available")
        return False


def check_hvantk_hgc() -> bool:
    """Check if hvantk.hgc functions are available."""
    try:
        from hvantk.hgc import (
            combine_gvcfs,
            convert_vds_to_mt,
            compute_full_qc,
            convert_mt_to_multi_sample_vcf,
        )

        print("✓ hvantk.hgc functions are available")
        return True
    except ImportError as e:
        print(f"✗ hvantk.hgc import failed: {e}")
        return False


def check_hail() -> bool:
    """Check if Hail is available and can be initialized."""
    try:
        import hail as hl

        print(f"✓ Hail is available (version: {hl.__version__})")
        return True
    except ImportError:
        print("✗ Hail is NOT available")
        return False


def check_gvcf_directory(gvcf_dir: str) -> bool:
    """Check if GVCF directory exists and contains files."""
    path = Path(gvcf_dir)
    if not path.exists():
        print(f"✗ GVCF directory does not exist: {gvcf_dir}")
        return False

    if not path.is_dir():
        print(f"✗ GVCF path is not a directory: {gvcf_dir}")
        return False

    # Check for GVCF files
    gvcf_patterns = ["*.g.vcf.gz", "*.gvcf.gz", "*.g.vcf", "*.gvcf"]
    gvcf_files = []
    for pattern in gvcf_patterns:
        gvcf_files.extend(list(path.glob(pattern)))

    if not gvcf_files:
        print(f"✗ No GVCF files found in: {gvcf_dir}")
        return False

    print(f"✓ GVCF directory exists with {len(gvcf_files)} files: {gvcf_dir}")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate HGC scalability benchmark setup"
    )
    parser.add_argument(
        "--gvcf-dir",
        default=os.environ.get("GVCF_DIR", "/path/to/your/gvcf_directory"),
        help="Path to GVCF directory (overrides env var GVCF_DIR)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    print("=" * 80)
    print("HGC Scalability Benchmark - Setup Validation")
    print("=" * 80)
    print()

    all_ok = True

    print("Checking Python version...")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        print(f"✓ Python {version.major}.{version.minor}.{version.micro}")
    else:
        print(
            f"✗ Python {version.major}.{version.minor}.{version.micro} (requires 3.8+)"
        )
        all_ok = False
    print()

    print("Checking required Python packages...")
    all_ok &= check_import("hail")
    all_ok &= check_import("pandas")
    all_ok &= check_import("numpy")
    all_ok &= check_import("matplotlib")
    all_ok &= check_import("seaborn")
    print()

    print("Checking hvantk library...")
    all_ok &= check_hvantk_hgc()
    print()

    print("Checking GVCF directory...")
    gvcf_dir = args.gvcf_dir
    if gvcf_dir == "/path/to/your/gvcf_directory":
        print(
            "  WARNING: GVCF_DIR is still a placeholder. Set --gvcf-dir or export GVCF_DIR."
        )
    gvcf_ok = check_gvcf_directory(gvcf_dir)
    if not gvcf_ok:
        print("  Note: You can specify a different directory with --gvcf-dir")
    all_ok &= gvcf_ok
    print()

    print("Checking benchmark scripts...")
    common_dir = Path(__file__).parent
    hgc_dir = common_dir.parent
    scalability_dir = hgc_dir / "scalability"
    cpu_scaling_dir = hgc_dir / "cpu_scaling"

    # Check scalability benchmark scripts
    print("\nScalability benchmark:")
    scalability_scripts = [
        ("benchmark.py", scalability_dir),
        ("benchmark.sh", scalability_dir),
        ("plot_results.py", scalability_dir),
    ]
    for script, location in scalability_scripts:
        script_path = location / script
        if script_path.exists():
            print(f"  ✓ {script} exists")
        else:
            print(f"  ✗ {script} NOT found")
            all_ok = False

    # Check CPU scaling benchmark scripts
    print("\nCPU scaling benchmark:")
    cpu_scaling_scripts = [
        ("benchmark.py", cpu_scaling_dir),
        ("benchmark.sh", cpu_scaling_dir),
        ("plot_results.py", cpu_scaling_dir),
    ]
    for script, location in cpu_scaling_scripts:
        script_path = location / script
        if script_path.exists():
            print(f"  ✓ {script} exists")
        else:
            print(f"  ✗ {script} NOT found")
            all_ok = False
    print()

    print("=" * 80)
    if all_ok:
        print("✓ All checks passed! Ready to run benchmarks.")
        print()
        print("To run scalability benchmark:")
        print(f"  cd {scalability_dir}")
        print("  bash benchmark.sh --gvcf-dir <path>")
        print()
        print("To run CPU scaling benchmark:")
        print(f"  cd {cpu_scaling_dir}")
        print("  bash benchmark.sh --gvcf-list <file> --output-dir <dir>")
    else:
        print("✗ Some checks failed. Please resolve the issues above.")
        print()
        print("Common fixes:")
        print(
            "  - Install missing packages: pip install hail pandas numpy matplotlib seaborn"
        )
        print("  - Check that hvantk is installed: pip install -e .")
        print(
            "  - Verify GVCF directory path is correct (or set --gvcf-dir / GVCF_DIR)"
        )
    print("=" * 80)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
