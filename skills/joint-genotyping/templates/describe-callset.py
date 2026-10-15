"""Describe a finished HGC pipeline run: the sample set, the MatrixTable row count and
the VCF record count. Runs inside the image (Hail). Writes <out>/counts.json.

    python describe-callset.py <output dir> [--prefix cohort]
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path


def vcf_records(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--prefix", default="cohort")
    args = ap.parse_args()
    out = Path(args.out)

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    mt = hl.read_matrix_table(str(out / f"{args.prefix}.mt"))
    samples = sorted(mt.s.collect())
    counts = {
        "samples": samples,
        "n_samples": len(samples),
        "mt_rows": mt.count_rows(),
        "vcf_records": vcf_records(out / f"{args.prefix}.vcf.bgz"),
        "vds_success": all(
            (out / f"{args.prefix}.vds" / part / "_SUCCESS").exists()
            for part in ("reference_data", "variant_data")
        ),
        "sample_qc_csv": (out / "qc" / f"{args.prefix}_sample_qc.csv").is_file(),
    }
    (out / "counts.json").write_text(json.dumps(counts, indent=2))
    for k, v in counts.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
