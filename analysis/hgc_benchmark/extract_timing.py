#!/usr/bin/env python3
"""
Extract timing data from JSON file and format as CSV row.

This script is called by hgc_scalability_benchmark.sh to extract timing data
from individual run JSON files and append them to the timings.csv file.

Usage:
    python3 extract_timing.py <timing_json_file> <sample_size>

Example:
    python3 extract_timing.py scalability_results/run_2/timing_2.json 2
    # Output: 2,45.2,12.3,8.5,15.7,81.7
"""

import json
import sys
from pathlib import Path


def main():
    """Extract timing data and format as CSV row."""
    if len(sys.argv) != 3:
        print("Usage: extract_timing.py <timing_json_file> <sample_size>", file=sys.stderr)
        print("", file=sys.stderr)
        print("Example:", file=sys.stderr)
        print("  python3 extract_timing.py scalability_results/run_2/timing_2.json 2", file=sys.stderr)
        sys.exit(1)

    timing_file = sys.argv[1]
    sample_size = sys.argv[2]

    # Required fields in the timing JSON
    required_fields = ["gvcf_combine", "vds_to_mt", "compute_qc", "mt_to_vcf", "total"]

    try:
        # Read and parse JSON file
        timing_path = Path(timing_file)
        if not timing_path.exists():
            print(f"ERROR: File not found: {timing_file}", file=sys.stderr)
            sys.exit(1)

        with open(timing_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Validate JSON structure
        if not isinstance(data, dict):
            print(f"ERROR: JSON file does not contain a dictionary: {timing_file}", file=sys.stderr)
            sys.exit(1)

        # Check for missing required fields
        missing = [k for k in required_fields if k not in data]
        if missing:
            print(f"ERROR: Missing required fields in {timing_file}: {missing}", file=sys.stderr)
            print(f"Available fields: {list(data.keys())}", file=sys.stderr)
            sys.exit(1)

        # Validate that values are numeric
        for field in required_fields:
            value = data[field]
            if not isinstance(value, (int, float)):
                print(f"ERROR: Field '{field}' has non-numeric value: {value} (type: {type(value).__name__})",
                      file=sys.stderr)
                sys.exit(1)

        # Format as CSV row: sample_size,gvcf_combine,vds_to_mt,compute_qc,mt_to_vcf,total
        row = f"{sample_size},{data['gvcf_combine']:.1f},{data['vds_to_mt']:.1f},{data['compute_qc']:.1f},{data['mt_to_vcf']:.1f},{data['total']:.1f}"

        # Print to stdout (will be captured by bash script)
        print(row)
        sys.exit(0)

    except FileNotFoundError:
        print(f"ERROR: File not found: {timing_file}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in {timing_file}: {e}", file=sys.stderr)
        print(f"File contents:", file=sys.stderr)
        try:
            with open(timing_file, 'r') as f:
                print(f.read()[:500], file=sys.stderr)  # Show first 500 chars
        except:
            pass
        sys.exit(1)
    except KeyError as e:
        print(f"ERROR: Missing key in JSON: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error processing {timing_file}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

