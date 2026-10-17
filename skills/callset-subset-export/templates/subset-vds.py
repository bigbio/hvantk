"""Subset a lossless VDS by sample list, or list its samples. Runs inside the image.

    python subset-vds.py --list-samples <callset.vds>
    python subset-vds.py --vds <callset.vds> --samples <list.txt> --out <subset.vds>

The subset keeps exactly the listed samples and drops alleles none of them carries, so
it equals a fresh joint call of those samples. Refuses to overwrite an existing output
and verifies the written sample set before exiting.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_list(path: Path) -> list[str]:
    names = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            names.append(s)
    return sorted(set(names))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list-samples", metavar="VDS")
    ap.add_argument("--vds")
    ap.add_argument("--samples")
    ap.add_argument("--out")
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    if args.list_samples:
        vds = hl.vds.read_vds(args.list_samples)
        for s in sorted(vds.variant_data.s.collect()):
            print(s)
        return 0
    if not (args.vds and args.samples and args.out):
        ap.error("--vds, --samples and --out are required (or --list-samples)")
    out = Path(args.out)
    if out.exists():
        print(f"ERROR: refusing to overwrite {out}", file=sys.stderr)
        return 1
    wanted = read_list(Path(args.samples))
    vds = hl.vds.read_vds(args.vds)
    present = set(vds.variant_data.s.collect())
    missing = [s for s in wanted if s not in present]
    if missing:
        print(
            f"ERROR: {len(missing)} listed sample(s) are not in the callset, e.g. {missing[:5]}; "
            "list the callset's spelling with --list-samples and map your list onto it",
            file=sys.stderr,
        )
        return 1
    sub = hl.vds.filter_samples(vds, wanted, keep=True, remove_dead_alleles=True)
    sub.write(str(out), overwrite=False)
    written = sorted(hl.vds.read_vds(str(out)).variant_data.s.collect())
    if written != wanted:
        print(
            f"ERROR: written sample set differs from the list: {written}",
            file=sys.stderr,
        )
        return 1
    print(f"subset written: {out} ({len(written)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
