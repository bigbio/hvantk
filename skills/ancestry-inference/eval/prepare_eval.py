"""Prepare the ancestry eval inputs from one 1000 Genomes chromosome: a labelled
reference, an unlabelled query of held-out samples, the held-out truth, and a reference
with permuted labels for the negative control. Deterministic. Runs inside the image.

    python prepare_eval.py --vcf <chr.vcf.gz> --samples-table <samples.ht> --out <stage dir>
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent / "templates" / "prepare-1kg-reference.py"


def load_template():
    spec = importlib.util.spec_from_file_location("prepare_1kg_reference", TEMPLATE)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--samples-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target", type=int, default=20000)
    args = ap.parse_args()
    out = Path(args.out)

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    full = load_template().build_reference(
        hl, [args.vcf], args.samples_table, "GRCh38", 0.05, args.target, "super_pop"
    )
    full = full.add_col_index("_cidx")
    query = full.filter_cols(full._cidx % 5 == 0).drop("_cidx")
    reference = full.filter_cols(full._cidx % 5 != 0).drop("_cidx")

    reference.write(str(out / "reference.mt"), overwrite=True)
    query.select_cols().write(str(out / "query.mt"), overwrite=True)
    query.cols().select("super_pop").export(str(out / "truth.tsv"))

    # Permuted labels for the negative control: shift the label list by half its length
    # over the sorted sample order, which breaks the sample-to-label link while keeping
    # the label counts.
    ref = hl.read_matrix_table(str(out / "reference.mt"))
    pairs = sorted(ref.cols().select("super_pop").collect(), key=lambda r: r.s)
    labels = [p.super_pop for p in pairs]
    shifted = labels[len(labels) // 2 :] + labels[: len(labels) // 2]
    mapping = hl.literal({p.s: lab for p, lab in zip(pairs, shifted)})
    ref.annotate_cols(super_pop=mapping[ref.s]).write(
        str(out / "reference_shuffled.mt"), overwrite=True
    )

    r = hl.read_matrix_table(str(out / "reference.mt"))
    q = hl.read_matrix_table(str(out / "query.mt"))
    print(f"reference: samples={r.count_cols()} variants={r.count_rows()}")
    print(f"query:     samples={q.count_cols()} variants={q.count_rows()}")
    print(r.aggregate_cols(hl.agg.counter(r.super_pop)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
