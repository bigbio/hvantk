#!/usr/bin/env python3
"""Keep one contig's rows of a tab-separated text file.

Used by ``stage.sh`` to cut the ml-features eval's inputs down to chromosome 22: the
ClinVar VCF, the dbNSFP variant file and the INSIDER BED. A line is kept when it is a
header or track line (it starts with ``#``, ``track`` or ``browser``) or when its first
tab-separated field is one of ``--contig``. Kept lines are copied byte for byte.

An output named ``*.bgz`` is written through Hail's ``hadoop_open``, so it is
block-gzipped, which ``hl.import_table(force_bgz=True)`` and ``hl.import_vcf`` read in
parallel (the ClinVar and dbNSFP builders expect it); any other name is written as plain
text (the INSIDER builder reads a plain ``.bed``). Runs inside the image, as a job (the
dbNSFP file is tens of GB).

    python slice_text.py --input IN[.gz|.bgz] --output OUT[.bgz] --contig 22 --contig chr22
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys

_HEADER_PREFIXES = (b"#", b"track", b"browser")


def _open_input(path: str):
    """A binary line stream over ``path``, through ``gzip -dc`` when it is compressed."""
    with open(path, "rb") as fh:
        compressed = fh.read(2) == b"\x1f\x8b"
    if not compressed:
        return open(path, "rb"), None
    gzip_bin = shutil.which("gzip")
    if gzip_bin:
        proc = subprocess.Popen(  # noqa: S603 -- a fixed program on a path we were given
            [gzip_bin, "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 20
        )
        return proc.stdout, proc
    return gzip.open(path, "rb"), None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True, help="*.bgz for BGZF, else plain text")
    ap.add_argument("--contig", action="append", required=True)
    args = ap.parse_args()

    if args.output.endswith(".bgz"):
        from hvantk.core.utils.hail_context import init_hail

        init_hail(quiet=True)
        import hail as hl

        sink = hl.hadoop_open(args.output, "wb")
    else:
        sink = open(args.output, "wb")
    contigs = {c.encode() for c in args.contig}
    stream, proc = _open_input(args.input)
    data = 0
    with stream, sink as out:
        for line in stream:
            if line.startswith(_HEADER_PREFIXES):
                out.write(line)
            elif line.split(b"\t", 1)[0] in contigs:
                out.write(line)
                data += 1
    if proc is not None and proc.wait() != 0:
        print(f"ERROR: gzip -dc {args.input} exited {proc.returncode}", file=sys.stderr)
        return 1
    print(f"{args.output}: {data} data rows on {sorted(args.contig)}")
    if data == 0:
        print("ERROR: no data rows kept; wrong contig spelling?", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
