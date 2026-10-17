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

The eQTL Catalogue's all-variant files repeat a (gene, variant) row once per extra
rsid, with identical statistics (a fixed hvantk defect: ``fetch_eqtl_region`` in
hvantk 0.5.0 and earlier returned every copy). This script de-duplicates on
``(gene_id, pos, ref, alt)`` itself, first row wins, before deriving anything from
``data`` — defensive, since a checkout predating the fix would otherwise double-count
some pairs in both the staged eQTL file and the pQTL rows derived from it.

Writes two files:

- ``--out-eqtl``: every (de-duplicated) variant-gene pair in the window, in the real
  GTEx v8 allpairs TSV layout (``gene_id variant_id tss_distance ma_samples ma_count
  maf pval_nominal slope slope_se``), gzipped. ``fetch_eqtl_region`` returns no
  ``maf``/``ma_samples``/``ma_count``/``tss_distance`` (and the GTEx v8 builder does
  not read them anyway), so those four columns are written as Hail's missing token,
  ``NA``.
- ``--out-pqtl``: a SYNTHETIC pQTL table in the real Fang et al. allpairs layout
  (``gene_name SNP CHR BP A1 NMISS BETA STAT P``, space-delimited, gzipped), derived
  from the eQTL rows by a fixed, deterministic rule (no randomness) so the two tables
  exercise every cascade class without depending on any non-public pQTL source. The
  construction below is internal to this staging script — nothing it writes states or
  implies that the pQTL file is anything other than a plain Fang-layout input.

Both files are written deterministically byte-for-byte run to run: rows are sorted
into a fixed order, gzipped with ``mtime=0`` (gzip's header otherwise embeds the
write time, which would change ``staged.sha256`` on every re-stage of identical
content), and written to a ``.part`` path then renamed into place.

Before writing anything, this script also simulates cascade.py's own outer-join
classification at p <= 5e-8 on both sides directly from the in-memory data (no Hail
needed) and prints the expected class counts and the number of coloc-candidate genes
(genes with a matching significant pair in both tables — the same gene set
``coloc_gene_ids()`` would return). It exits 1 if any class would be empty or there
is no coloc candidate, so a bad window is caught here, in seconds, rather than after
three cluster jobs.

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
import sys
from pathlib import Path

from scipy.stats import norm

DATASET = "QTD000266"  # GTEx liver, gene expression
STUDY = "QTS000015"  # GTEx, in the eQTL Catalogue
CHROM = "1"
START, END = 108_500_000, 110_500_000  # GRCh38
CASCADE_P_THRESHOLD = 5e-8  # cascade.py's own default, both sides

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


def _atomic_gzip_write(path: Path, text: str) -> None:
    """Write *text* to *path* gzipped, via a ``.part`` file renamed into place.

    Identical content gives identical bytes run to run: ``mtime=0`` drops the
    write-time timestamp gzip otherwise embeds in its header, and ``filename=""``
    (with the destination handed in as ``fileobj=`` instead of as the
    ``GzipFile``'s own ``filename=``) drops the embedded FNAME field, which would
    otherwise bake in the ``.part`` suffix of the intermediate path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.parent / (path.name + ".part")
    with (
        open(part, "wb") as raw,
        gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9
        ) as gz,
    ):
        gz.write(text.encode("utf-8"))
    part.replace(path)


def fetch_window():
    """``{gene_id: [((pos, ref, alt), beta, se, pval), ...]}`` for the staged window.

    ``fetch_eqtl_region`` already strips the gene-id version suffix and drops rows
    with ``se <= 0``; there is no ``maf`` in this payload.
    """
    from hvantk.algorithms.qtlcascade.gwas_coloc import fetch_eqtl_region

    return fetch_eqtl_region(DATASET, CHROM, START, END, study=STUDY)


def dedupe_window(data: dict) -> dict:
    """De-duplicate on ``(gene_id, pos, ref, alt)``, first row wins.

    Defensive: the eQTL Catalogue's all-variant files repeat a (gene, variant) row
    once per extra rsid, with identical statistics. ``fetch_eqtl_region`` itself
    de-duplicates this in hvantk after the fix (issue: a 2 Mb GTEx liver window
    carried 6,787 such repeats in 205,386 rows), but this script does it again so a
    checkout predating that fix still stages clean, non-doubled data.
    """
    out: dict = {}
    for gene_id, variants in data.items():
        seen: set = set()
        kept = []
        for key, beta, se, pval in variants:
            if key in seen:
                continue
            seen.add(key)
            kept.append((key, beta, se, pval))
        out[gene_id] = kept
    return out


def write_eqtl_raw(data: dict, out_path: Path) -> int:
    rows = []
    for gene_id, variants in data.items():
        for (pos, ref, alt), beta, se, pval in variants:
            rows.append((gene_id, pos, ref, alt, beta, se, pval))
    rows.sort(key=lambda r: (r[0], r[1]))  # gene, then position

    lines = ["\t".join(EQTL_HEADER)]
    for gene_id, pos, ref, alt, beta, se, pval in rows:
        lines.append(
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
        )
    _atomic_gzip_write(out_path, "\n".join(lines) + "\n")
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
    """Synthetic Fang-layout rows; returns ``(rows, class_counts, by_key)``.

    Genes with a symbol, sorted by gene id, take class ``i % 4``:

    - 0 concordant: ``BETA = 0.8 * slope``, same SE as the eQTL row.
    - 1 discordant: ``BETA = -0.8 * slope``, same SE as the eQTL row.
    - 2: no pQTL rows at all for this gene (eQTL-only by construction).
    - 3 shifted: each variant's BETA/SE are the eQTL slope/SE of the variant
      ``SHIFT`` positions later in that gene's position-sorted list (wrapping) —
      same physical SNP (SNP/CHR/BP/A1 unchanged), unrelated effect size.

    ``by_key`` maps ``(gene_id, pos, ref, alt)`` (the same key space as the eQTL
    side) to ``(BETA, STAT, P)``, for ``expected_cascade_classes`` below — the
    written rows themselves carry the gene *symbol*, not ``gene_id``, so this is
    computed here rather than reconstructed from them afterwards.
    """
    matched_genes = sorted(g for g in data if g in ensembl_to_symbol)
    class_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    rows = []
    by_key: dict = {}
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
            by_key[(gene_id, pos, ref, alt)] = (eff_beta, stat, p_value)
    rows.sort(key=lambda r: (r[0], r[3]))  # gene_name, then position
    return rows, class_counts, by_key


def write_pqtl_raw(rows, out_path: Path) -> None:
    lines = [" ".join(PQTL_HEADER)]
    for gene_name, snp, chrom, bp, a1, nmiss, beta, stat, p_value in rows:
        lines.append(
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
        )
    _atomic_gzip_write(out_path, "\n".join(lines) + "\n")


def expected_cascade_classes(
    data: dict, pqtl_by_key: dict, p_threshold: float = CASCADE_P_THRESHOLD
):
    """Simulate ``cascade.py``'s own outer-join classification at *p_threshold* on
    both sides, directly from this script's in-memory data (no Hail): a fast,
    Hail-free sanity check that the staged inputs will actually exercise every
    cascade class, and give coloc a candidate gene, before a cluster job runs them
    for real. Returns ``(class_counts, coloc_gene_ids)``.
    """
    eqtl_sig = {}
    for gene_id, variants in data.items():
        for (pos, ref, alt), beta, se, pval in variants:
            if pval <= p_threshold:
                eqtl_sig[(gene_id, pos, ref, alt)] = beta

    pqtl_sig = {
        key: beta
        for key, (beta, stat, p_value) in pqtl_by_key.items()
        # the pqtl builder drops STAT == 0 rows (SE = |BETA/STAT| undefined)
        # before anything downstream ever sees them.
        if stat != 0.0 and p_value <= p_threshold
    }

    class_counts = {"eqtl_mediated": 0, "discordant": 0, "eqtl_only": 0, "pqtl_only": 0}
    coloc_gene_ids: set = set()
    for key in set(eqtl_sig) | set(pqtl_sig):
        has_eqtl = key in eqtl_sig
        has_pqtl = key in pqtl_sig
        if has_eqtl and has_pqtl:
            cls = "eqtl_mediated" if eqtl_sig[key] * pqtl_sig[key] > 0 else "discordant"
            coloc_gene_ids.add(key[0])
        elif has_eqtl:
            cls = "eqtl_only"
        else:
            cls = "pqtl_only"
        class_counts[cls] += 1
    return class_counts, coloc_gene_ids


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True, help="staged hgnc_complete_set.txt")
    ap.add_argument("--out-eqtl", required=True)
    ap.add_argument("--out-pqtl", required=True)
    args = ap.parse_args()

    print(f"fetching {DATASET} (study {STUDY}) chr{CHROM}:{START}-{END} ...")
    raw_data = fetch_window()
    data = dedupe_window(raw_data)
    n_raw_rows = sum(len(v) for v in raw_data.values())
    n_rows = sum(len(v) for v in data.values())
    if n_rows != n_raw_rows:
        print(
            f"de-duplicated {n_raw_rows - n_rows} repeated (gene_id, pos, ref, alt) "
            "row(s) from the eQTL Catalogue fetch",
        )
    n_genes = len(data)

    n_eqtl_rows = write_eqtl_raw(data, Path(args.out_eqtl))
    print(f"eQTL window: {n_genes} genes, {n_eqtl_rows} variant-gene pairs")

    ensembl_to_symbol = load_ensembl_to_symbol(Path(args.hgnc))
    pqtl_rows, class_counts, pqtl_by_key = build_pqtl_rows(data, ensembl_to_symbol)
    write_pqtl_raw(pqtl_rows, Path(args.out_pqtl))
    matched = sum(class_counts.values())
    print(
        f"pQTL genes: {matched} of {n_genes} matched an approved HGNC symbol; "
        f"class counts (0=concordant,1=discordant,2=none,3=shifted): {class_counts}"
    )
    print(f"pQTL rows written: {len(pqtl_rows)}")
    print(f"wrote {args.out_eqtl}")
    print(f"wrote {args.out_pqtl}")

    expected_counts, coloc_gene_ids = expected_cascade_classes(data, pqtl_by_key)
    print(
        f"expected cascade class counts at p<={CASCADE_P_THRESHOLD:g}: {expected_counts}"
    )
    print(
        f"expected coloc-candidate genes (a significant, matching pair in both "
        f"tables): {len(coloc_gene_ids)}"
    )
    empty_classes = sorted(c for c, n in expected_counts.items() if n == 0)
    if empty_classes:
        print(
            f"ERROR: expected cascade class count is 0 for: {empty_classes} -- this "
            "window will not exercise every cascade class",
            file=sys.stderr,
        )
        return 1
    if not coloc_gene_ids:
        print(
            "ERROR: no coloc-candidate gene (a significant, matching pair in both "
            "tables) -- colocalisation would have nothing to test",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
