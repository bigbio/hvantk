"""Build the cohort-attach-burden eval's SYNTHETIC case/control cohort MatrixTable from
one 1000 Genomes high-coverage chromosome VCF and the staged gene index: a fixed
case/control split of 1,000 samples, rare (AC>0, AF<0.01 AMONG those 1,000 samples)
split-biallelic variants, each annotated with the smallest-id gene whose body it falls
inside and a snv/indel route. A staging helper next to the eval's five fixed files --
not read by the task, not graded directly by check.py.

    python prepare_cohort.py --vcf <chr22.vcf.gz> --gene-index <gene-index.ht> \
        --out <cohort.mt> --summary <cohort.summary.json>

Runs inside the image (Hail). Deterministic given the same VCF and gene index: the
case/control split is a fixed rule over the sorted sample ids (rank parity), not a
random draw.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


# --- canonical digest (COMMON.md's shared convention; duplicated here, not imported,
# --- because this staging helper and the skill's eval/check.py are graded/read
# --- independently and must each stay self-contained) -------------------------------


def norm(value):
    """JSON-stable form for a digest: floats become a fixed-precision string (stable
    across platforms, finite for nan/inf), lists/tuples map elementwise, sets become
    sorted lists, dicts/Hail Structs become dicts, Locus/Interval render via str()."""
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, float):
        return format(value, ".10g")
    if isinstance(value, (int, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [norm(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(norm(v) for v in value)
    if isinstance(value, dict):
        return {k: norm(v) for k, v in value.items()}
    try:
        import hail as hl

        if isinstance(value, hl.Struct):
            return {k: norm(v) for k, v in dict(value).items()}
        if isinstance(value, (hl.Locus, hl.Interval)):
            return str(value)
    except ImportError:
        pass
    return str(value)


def rows_digest(rows) -> str:
    strings = sorted(
        json.dumps(norm(dict(r)), sort_keys=True, separators=(",", ":")) for r in rows
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--gene-index", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--n-samples", type=int, default=1000)
    ap.add_argument("--max-af", type=float, default=0.01)
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    mt = hl.import_vcf(
        args.vcf,
        force_bgz=True,
        reference_genome="GRCh38",
        array_elements_required=False,
    )

    # Fixed case/control split of the first N samples by sorted sample id: deterministic
    # and reproducible across stagings, with no external label file or randomness.
    ids = sorted(mt.s.collect())
    if len(ids) < args.n_samples:
        raise SystemExit(
            f"VCF has only {len(ids)} samples, fewer than --n-samples {args.n_samples}"
        )
    kept = ids[: args.n_samples]
    is_case_of = {s: (i % 2 == 0) for i, s in enumerate(kept)}
    mt = mt.filter_cols(hl.literal(set(kept), hl.tset(hl.tstr)).contains(mt.s))
    mt = mt.annotate_cols(
        is_case=hl.literal(is_case_of, hl.tdict(hl.tstr, hl.tbool))[mt.s]
    )

    mt = hl.split_multi_hts(mt)

    # Rare among the kept 1,000 samples only: call_stats is computed AFTER the column
    # filter above, so AC/AF reflect exactly this cohort, not the full 1000-Genomes
    # high-coverage release's 3,202 samples.
    mt = mt.annotate_rows(_cs=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.filter_rows((mt._cs.AC[1] > 0) & (mt._cs.AF[1] < args.max_af))

    mt = mt.annotate_rows(
        route=hl.if_else(hl.is_snp(mt.alleles[0], mt.alleles[1]), "snv", "indel")
    )

    # Gene-body overlap against the staged gene index. The gene index's `chromosome`
    # values carry no "chr" prefix (Ensembl GTF convention --
    # hvantk/skills/ensembl_gene/structure/parse.py: chromosome is the GTF's own column
    # 1, verbatim); the 1000 Genomes high-coverage VCF's contigs do ("chr22"). Handled
    # both ways so a differently-built gene index still matches. Pre-filtered to this
    # VCF's own chromosome(s) and broadcast as a literal: chr22 alone has a few hundred
    # protein-coding genes, cheap to scan per variant on the driver; collecting the
    # whole ~20,000-gene index would not be.
    def with_chr(chrom: str) -> str:
        return chrom if chrom.startswith("chr") else "chr" + chrom

    vcf_contigs = {
        with_chr(c) for c in mt.aggregate_rows(hl.agg.collect_as_set(mt.locus.contig))
    }
    gi = hl.read_table(args.gene_index).select("chromosome", "gene_start", "gene_end")
    genes = [g for g in gi.collect() if with_chr(g.chromosome) in vcf_contigs]
    if not genes:
        raise SystemExit(
            f"no gene-index genes overlap this VCF's chromosome(s) {sorted(vcf_contigs)}"
        )
    gene_struct_type = hl.tstruct(
        gene_id=hl.tstr, chromosome=hl.tstr, gene_start=hl.tint32, gene_end=hl.tint32
    )
    gene_structs = [
        hl.Struct(
            gene_id=g.gene_id,
            chromosome=with_chr(g.chromosome),
            gene_start=g.gene_start,
            gene_end=g.gene_end,
        )
        for g in genes
    ]
    gl = hl.literal(gene_structs, dtype=hl.tarray(gene_struct_type))
    matches = gl.filter(
        lambda g: (
            (g.chromosome == mt.locus.contig)
            & (g.gene_start <= mt.locus.position)
            & (mt.locus.position <= g.gene_end)
        )
    )
    mt = mt.annotate_rows(_gene_ids=matches.map(lambda g: g.gene_id))
    mt = mt.annotate_rows(
        gene_id=hl.if_else(
            hl.len(mt._gene_ids) > 0, hl.sorted(mt._gene_ids)[0], hl.missing(hl.tstr)
        )
    )
    mt = mt.filter_rows(hl.is_defined(mt.gene_id))

    mt = mt.select_rows("gene_id", "route")
    mt = mt.select_cols("is_case")
    mt = mt.select_entries("GT")
    mt = mt.naive_coalesce(32)
    mt.write(args.out, overwrite=True)

    mt = hl.read_matrix_table(args.out)
    n_rows = mt.count_rows()
    n_cols = mt.count_cols()
    arm = mt.aggregate_cols(
        hl.struct(
            n_case=hl.agg.count_where(mt.is_case),
            n_control=hl.agg.count_where(~mt.is_case),
        )
    )
    genes_covered = mt.aggregate_rows(hl.agg.collect_as_set(mt.gene_id))
    digest = rows_digest(mt.rows().select("gene_id", "route").collect())
    summary = {
        "n_rows": n_rows,
        "n_cols": n_cols,
        "n_case": arm.n_case,
        "n_control": arm.n_control,
        "n_genes": len(genes_covered),
        "rows_sha256": digest,
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"cohort.mt written: {args.out}")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
