#!/usr/bin/env python3
"""
Test script to validate HGC scalability benchmark setup.

This script checks that all required dependencies are available
and the configuration is correct.
"""

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
            convert_mt_to_multi_sample_vcf
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
    gvcf_patterns = ['*.g.vcf.gz', '*.gvcf.gz', '*.g.vcf', '*.gvcf']
    gvcf_files = []
    for pattern in gvcf_patterns:
        gvcf_files.extend(list(path.glob(pattern)))

    if not gvcf_files:
        print(f"✗ No GVCF files found in: {gvcf_dir}")
        return False

    print(f"✓ GVCF directory exists with {len(gvcf_files)} files: {gvcf_dir}")
    return True

def main():
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
        print(f"✗ Python {version.major}.{version.minor}.{version.micro} (requires 3.8+)")
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
    gvcf_dir = "/mnt/nfs/KOL_UOL/projects/CHD_1000WGS/variant_calling/split_vcfs/chr20"
    gvcf_ok = check_gvcf_directory(gvcf_dir)
    if not gvcf_ok:
        print(f"  Note: You can specify a different directory with --gvcf-dir")
    all_ok &= gvcf_ok
    print()

    print("Checking benchmark scripts...")
    script_dir = Path(__file__).parent
    scripts = [
        "hgc_scalability_benchmark.sh",
        "hgc_scalability_benchmark.py",
        "plot_scalability_results.py"
    ]
    for script in scripts:
        script_path = script_dir / script
        if script_path.exists():
            print(f"✓ {script} exists")
        else:
            print(f"✗ {script} NOT found")
            all_ok = False
    print()

    print("=" * 80)
    if all_ok:
        print("✓ All checks passed! Ready to run benchmark.")
        print()
        print("To start the benchmark, run:")
        print(f"  cd {script_dir}")
        print("  bash hgc_scalability_benchmark.sh")
    else:
        print("✗ Some checks failed. Please resolve the issues above.")
        print()
        print("Common fixes:")
        print("  - Install missing packages: pip install hail pandas numpy matplotlib seaborn")
        print("  - Check that hvantk is installed: pip install -e .")
        print("  - Verify GVCF directory path is correct")
    print("=" * 80)

    return 0 if all_ok else 1

if __name__ == '__main__':
    sys.exit(main())

