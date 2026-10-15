"""Stage the qtl-cascade eval's eQTL and pQTL raw files from one real, public eQTL
Catalogue window (GTEx liver gene expression, dataset QTD000266, study QTS000015,
GRCh38 chr1:108,500,000-110,500,000 — about 2 Mb, picked only for being a quiet,
reasonably gene-dense region). Runs on the LOGIN NODE, through `run-in-image.sh` (it
needs egress; it does not need Hail): it calls hvantk's own
``hvantk.algorithms.qtlcascade.gwas_coloc.fetch_eqtl_region``, the same remote-tabix
code path ``qtlcascade gwas-coloc`` uses. pysam's remote tabix downloads the
remote index (a ``.tbi``) into the CURRENT WORKING DIRECTORY, not a temp
location -- ``stage.sh`` runs this script from a scratch directory for exactly
that reason and removes the stray file afterwards; this script itself does not
need to care where it is run from.

Writes two files:

- ``--out-eqtl``: every variant-gene pair in the window, in the real GTEx v8 allpairs
  TSV layout (``gene_id variant_id tss_distance ma_samples ma_count maf pval_nominal
  slope slope_se``), gzipped. ``fetch_eqtl_region`` returns no ``maf``/``ma_samples``/
  ``ma_count``/``tss_distance`` (and the GTEx v8 builder does not read them anyway),
  so those four columns are written as Hail's missing token, ``NA``.
- ``--out-pqtl``: a SYNTHETIC pQTL table in the real Fang et al. allpairs layout
  (``gene_name SNP CHR BP A1 NMISS BETA STAT P``, space-delimited, gzipped), derived
  from the eQTL rows by a fixed, deterministic rule (no randomness) so the two tables
  exercise every cascade class without depending on any non-public pQTL source. The
  construction below is internal to this staging script — nothing it writes states or
  implies that the pQTL file is anything other than a plain Fang-layout input.

Gene symbols for the pQTL file come from the staged HGNC raw TSV (``--hgnc``,
``hgnc_id``/``symbol``/``ensembl_gene_id``/``status`` columns, approved rows only) —
the real mapping ``pqtl:metrics --plugin-arg hgnc_ht=...`` needs at build time, run
here in reverse (Ensembl ID -> symbol) to produce symbol-keyed Fang rows that build
back into the *same* Ensembl IDs.

    python stage_qtl.py --hgnc RAW/hgnc/hgnc_complete_set.txt \
        --out-eqtl RAW/gtex-eqtl/Liver.v8.allpairs.txt.gz \
        --out-pqtl RAW/pqtl/Liver.allpairs_nobsGE72.txt.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from scipy.stats import norm

DATASET = "QTD000266"  # GTEx liver, gene expression
STUDY = "QTS000015"  # GTEx, in the eQTL Catalogue
CHROM = "1"
START, END = 108_500_000, 110_500_000  # GRCh38

EQTL_HEADER = (
    "gene_id",
    "variant_id",
    "tss_distance",
    "ma_samples",
    "ma_count",
    "maf",
    "pval_nominal",
    "slope",
    "slope_se",
)
PQTL_HEADER = ("gene_name", "SNP", "CHR", "BP", "A1", "NMISS", "BETA", "STAT", "P")
SHIFT = 25
NMISS = 200


def _variant_id(pos: int, ref: str, alt: str) -> str:
    return f"chr{CHROM}_{pos}_{ref}_{alt}_b38"


def _fmt(x: float) -> str:
    return format(x, ".10g")


def fetch_window():
    """``{gene_id: [((pos, ref, alt), beta, se, pval), ...]}`` for the staged window.

    ``fetch_eqtl_region`` already strips the gene-id version suffix and drops rows
    with ``se <= 0``; there is no ``maf`` in this payload.
    """
    from hvantk.algorithms.qtlcascade.gwas_coloc import fetch_eqtl_region

    return fetch_eqtl_region(DATASET, CHROM, START, END, study=STUDY)


def write_eqtl_raw(data: dict, out_path: Path) -> int:
    rows = []
    for gene_id, variants in data.items():
        for (pos, ref, alt), beta, se, pval in variants:
            rows.append((gene_id, pos, ref, alt, beta, se, pval))
    rows.sort(key=lambda r: (r[0], r[1]))  # gene, then position

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as f:
        f.write("\t".join(EQTL_HEADER) + "\n")
        for gene_id, pos, ref, alt, beta, se, pval in rows:
            f.write(
                "\t".join(
                    (
                        gene_id,
                        _variant_id(pos, ref, alt),
                        "NA",  # tss_distance
                        "NA",  # ma_samples
                        "NA",  # ma_count
                        "NA",  # maf
                        _fmt(pval),
                        _fmt(beta),  # slope
                        _fmt(se),  # slope_se
                    )
                )
                + "\n"
            )
    return len(rows)


def load_ensembl_to_symbol(hgnc_path: Path) -> dict:
    """``{ensembl_gene_id: symbol}`` from approved HGNC rows (first occurrence wins;
    file order is fixed, so this is deterministic)."""
    mapping: dict = {}
    with open(hgnc_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("status") != "Approved":
                continue
            ensembl_id = (row.get("ensembl_gene_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            if ensembl_id and symbol and ensembl_id not in mapping:
                mapping[ensembl_id] = symbol
    return mapping


def build_pqtl_rows(data: dict, ensembl_to_symbol: dict):
    """Synthetic Fang-layout rows; returns ``(rows, class_counts)``.

    Genes with a symbol, sorted by gene id, take class ``i % 4``:

    - 0 concordant: ``BETA = 0.8 * slope``, same SE as the eQTL row.
    - 1 discordant: ``BETA = -0.8 * slope``, same SE as the eQTL row.
    - 2: no pQTL rows at all for this gene (eQTL-only by construction).
    - 3 shifted: each variant's BETA/SE are the eQTL slope/SE of the variant
      ``SHIFT`` positions later in that gene's position-sorted list (wrapping) —
      same physical SNP (SNP/CHR/BP/A1 unchanged), unrelated effect size.
    """
    matched_genes = sorted(g for g in data if g in ensembl_to_symbol)
    class_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    rows = []
    for i, gene_id in enumerate(matched_genes):
        cls = i % 4
        class_counts[cls] += 1
        if cls == 2:
            continue
        symbol = ensembl_to_symbol[gene_id]
        variants = sorted(data[gene_id], key=lambda t: t[0][0])
        n = len(variants)
        for j, ((pos, ref, alt), beta, se, _pval) in enumerate(variants):
            if cls == 0:
                eff_beta, eff_se = 0.8 * beta, se
            elif cls == 1:
                eff_beta, eff_se = -0.8 * beta, se
            else:  # cls == 3
                _, shift_beta, shift_se, _ = variants[(j + SHIFT) % n]
                eff_beta, eff_se = shift_beta, shift_se
            stat = eff_beta / eff_se
            p_value = 2 * norm.sf(abs(stat))
            rows.append(
                (
                    symbol,
                    _variant_id(pos, ref, alt),
                    CHROM,
                    pos,
                    alt,
                    NMISS,
                    eff_beta,
                    stat,
                    p_value,
                )
            )
    rows.sort(key=lambda r: (r[0], r[3]))  # gene_name, then position
    return rows, class_counts


def write_pqtl_raw(rows, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as f:
        f.write(" ".join(PQTL_HEADER) + "\n")
        for gene_name, snp, chrom, bp, a1, nmiss, beta, stat, p_value in rows:
            f.write(
                " ".join(
                    (
                        gene_name,
                        snp,
                        chrom,
                        str(bp),
                        a1,
                        str(nmiss),
                        _fmt(beta),
                        _fmt(stat),
                        _fmt(p_value),
                    )
                )
                + "\n"
            )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True, help="staged hgnc_complete_set.txt")
    ap.add_argument("--out-eqtl", required=True)
    ap.add_argument("--out-pqtl", required=True)
    args = ap.parse_args()

    print(f"fetching {DATASET} (study {STUDY}) chr{CHROM}:{START}-{END} ...")
    data = fetch_window()
    n_genes = len(data)
    n_eqtl_rows = write_eqtl_raw(data, Path(args.out_eqtl))
    print(f"eQTL window: {n_genes} genes, {n_eqtl_rows} variant-gene pairs")

    ensembl_to_symbol = load_ensembl_to_symbol(Path(args.hgnc))
    pqtl_rows, class_counts = build_pqtl_rows(data, ensembl_to_symbol)
    write_pqtl_raw(pqtl_rows, Path(args.out_pqtl))
    matched = sum(class_counts.values())
    print(
        f"pQTL genes: {matched} of {n_genes} matched an approved HGNC symbol; "
        f"class counts (0=concordant,1=discordant,2=none,3=shifted): {class_counts}"
    )
    print(f"pQTL rows written: {len(pqtl_rows)}")
    print(f"wrote {args.out_eqtl}")
    print(f"wrote {args.out_pqtl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
