"""Build a labelled reference MatrixTable from 1000 Genomes chromosome VCFs and the
sample table. Runs inside the image (Hail).

    python prepare-1kg-reference.py --vcf <chr.vcf.gz> [--vcf ...] --samples-table <samples.ht> \
        --out <reference.mt> [--build GRCh38] [--min-maf 0.05] [--target 20000] [--label-col super_pop]

Keeps biallelic SNPs, joins the labels (samples without one are dropped, which removes
the related samples of the 3202-sample release), filters on minor allele frequency among
the labelled samples, thins deterministically to about --target variants spread along
the genome, and writes the reference with the label as a flat column.
"""

from __future__ import annotations

import argparse


def build_reference(hl, vcfs, samples_table, build, min_maf, target, label_col):
    mt = hl.import_vcf(
        vcfs, force_bgz=True, reference_genome=build, array_elements_required=False
    )
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    ht = hl.read_table(samples_table)
    mt = mt.annotate_cols(**{label_col: ht[mt.s][label_col]})
    mt = mt.filter_cols(hl.is_defined(mt[label_col]))
    mt = mt.select_cols(label_col).select_entries("GT").select_rows()
    mt = hl.variant_qc(mt)
    af = mt.variant_qc.AF[1]
    mt = mt.filter_rows((af >= min_maf) & (af <= 1 - min_maf)).select_rows()
    n = mt.count_rows()
    k = max(1, n // target)
    mt = mt.add_row_index("_idx")
    mt = mt.filter_rows(mt._idx % k == 0).drop("_idx")
    return mt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", action="append", required=True)
    ap.add_argument("--samples-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--min-maf", type=float, default=0.05)
    ap.add_argument("--target", type=int, default=20000)
    ap.add_argument("--label-col", default="super_pop")
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    mt = build_reference(
        hl,
        args.vcf,
        args.samples_table,
        args.build,
        args.min_maf,
        args.target,
        args.label_col,
    )
    mt.write(args.out, overwrite=False)
    mt = hl.read_matrix_table(args.out)
    print(
        f"reference written: {args.out} samples={mt.count_cols()} variants={mt.count_rows()}"
    )
    print(mt.aggregate_cols(hl.agg.counter(mt[args.label_col])))


if __name__ == "__main__":
    main()
