"""Preflight check for a cohort MatrixTable before spending a `hvantk cohort burden`
job on it: the clean-input contract (split bi-allelic rows, a defined gene field, a
route field, a boolean arm column defined for every sample -- see
`hvantk/algorithms/burden/aggregate.py`'s `assert_clean_mt`), plus descriptive counts
that help pick `--carrier-mode` and `--key`. Runs inside the image (Hail).

    python inspect-cohort-mt.py --mt M.mt --gene-col gene_id --route-col route \
        --arm-col is_case [--score-field revel]

Prints row/column counts, every field and its type, multiallelic row count, rows with
an undefined gene, a rough shape breakdown of the sampled gene ids (unversioned ENSG /
versioned ENSG / symbol-like), per-arm sample counts, and route value counts. Exits 1
on a missing field, a non-boolean or partly undefined arm column, or any multiallelic
row -- the same contract `hvantk cohort burden` itself enforces, checked here before a
real job is spent on an input that would just be rejected.
"""

from __future__ import annotations

import argparse
import re
import sys

_ENSG_UNVERSIONED = re.compile(r"^ENSG\d+$")
_ENSG_VERSIONED = re.compile(r"^ENSG\d+\.\d+$")
_SAMPLE_GENE_IDS = 2000


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mt", required=True)
    ap.add_argument("--gene-col", required=True)
    ap.add_argument("--route-col", required=True)
    ap.add_argument("--arm-col", required=True)
    ap.add_argument(
        "--score-field", default=None, help="optional; checked for presence only"
    )
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    mt = hl.read_matrix_table(args.mt)
    row_fields = list(mt.row)
    col_fields = list(mt.col)
    entry_fields = list(mt.entry)

    missing = []
    for req in ("locus", "alleles"):
        if req not in row_fields:
            missing.append(f"row field {req!r} not found")
    if args.gene_col not in row_fields:
        missing.append(
            f"gene column {args.gene_col!r} not found in row fields {row_fields}"
        )
    if args.route_col not in row_fields:
        missing.append(
            f"route column {args.route_col!r} not found in row fields {row_fields}"
        )
    if "GT" not in entry_fields:
        missing.append(f"entry field 'GT' not found in entry fields {entry_fields}")
    if args.arm_col not in col_fields:
        missing.append(
            f"arm column {args.arm_col!r} not found in column fields {col_fields}"
        )
    if args.score_field and args.score_field not in row_fields:
        missing.append(
            f"score field {args.score_field!r} not found in row fields {row_fields}"
        )
    if missing:
        print("MISSING FIELD(S):")
        for m in missing:
            print(f"  - {m}")
        return 1

    n_rows = mt.count_rows()
    n_cols = mt.count_cols()
    print(f"rows (variants): {n_rows}")
    print(f"cols (samples):  {n_cols}")
    print("row fields:   " + ", ".join(f"{f}:{mt[f].dtype}" for f in row_fields))
    print("col fields:   " + ", ".join(f"{f}:{mt[f].dtype}" for f in col_fields))
    print("entry fields: " + ", ".join(f"{f}:{mt[f].dtype}" for f in entry_fields))

    n_multi = mt.filter_rows(hl.len(mt.alleles) != 2).count_rows()
    print(f"multiallelic rows (len(alleles) != 2): {n_multi}")

    n_undef_gene = mt.filter_rows(hl.is_missing(mt[args.gene_col])).count_rows()
    print(f"rows with an undefined gene ({args.gene_col}): {n_undef_gene}")

    sampled = mt.aggregate_rows(
        hl.agg.filter(
            hl.is_defined(mt[args.gene_col]),
            hl.agg.take(mt[args.gene_col], _SAMPLE_GENE_IDS),
        )
    )
    n_unversioned = sum(1 for g in sampled if _ENSG_UNVERSIONED.match(g))
    n_versioned = sum(1 for g in sampled if _ENSG_VERSIONED.match(g))
    n_symbol_like = len(sampled) - n_unversioned - n_versioned
    print(
        f"gene id shape (of {len(sampled)} sampled values): "
        f"unversioned ENSG={n_unversioned} versioned ENSG={n_versioned} "
        f"symbol-like={n_symbol_like}"
    )

    arm = mt.aggregate_cols(
        hl.struct(
            n=hl.agg.count(),
            n_defined=hl.agg.count_where(hl.is_defined(mt[args.arm_col])),
            n_true=hl.agg.count_where(mt[args.arm_col]),
        )
    )
    n_undefined_arm = arm.n - arm.n_defined
    n_control = arm.n_defined - arm.n_true
    print(
        f"arm column {args.arm_col!r}: {arm.n_true} case, {n_control} control, "
        f"{n_undefined_arm} undefined (of {arm.n} samples)"
    )

    if str(mt[args.route_col].dtype).startswith("array<"):
        route_counts = mt.aggregate_rows(
            hl.agg.explode(lambda r: hl.agg.counter(r), mt[args.route_col])
        )
    else:
        route_counts = mt.aggregate_rows(hl.agg.counter(mt[args.route_col]))
    print(f"route ({args.route_col}) value counts: {dict(route_counts)}")

    problems = []
    if n_multi:
        problems.append(
            f"{n_multi} multiallelic row(s); run hl.split_multi_hts upstream"
        )
    if n_undefined_arm:
        problems.append(
            f"arm column {args.arm_col!r} is undefined for {n_undefined_arm} "
            "sample(s); define it for every sample before burden"
        )
    if arm.n_true in (0, arm.n):
        problems.append(
            f"arm column {args.arm_col!r} is not binary: {arm.n_true} case / "
            f"{n_control} control"
        )
    if problems:
        print("FAIL:")
        for p in problems:
            print(f"  - {p}")
        return 1

    print("OK: clean-input contract satisfied")
    return 0


if __name__ == "__main__":
    sys.exit(main())
