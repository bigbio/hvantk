"""Summarise one tissue's QTL cascade outputs under a `qtlcascade run`/`run_collection`
output directory: the cascade class counts, the cascade gene count, how many genes have
a complete cascade, how many genes were coloc-tested and colocalised (H4 > 0.8), the
top five genes by H4, the sha256 of `gene_summary.tsv` and `coloc_results.tsv`, and
whether the HTML report is present. Prints each value and writes the same as JSON to
`--receipt`. Runs inside the image: Hail to read `cascade.ht`, pandas for the TSVs.

    python summarize-cascade.py --out-dir OUT_DIR --tissue Liver --receipt OUT_DIR/summary.json

Exits 1 (after printing which path is missing) when `cascade.ht`, `gene_summary.tsv`,
or the report is absent — those three always exist after any completed run.
`coloc_results.tsv` is optional: a cascade run without `--eqtl-allpairs`/
`--pqtl-allpairs` never writes it, so its absence is reported, not an error.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

H4_THRESHOLD = 0.8


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _as_bool(value) -> bool:
    """`has_complete_cascade` arrives from a Hail TSV export as the literal string
    ``"true"``/``"false"`` (lowercase), not Python's ``"True"``/``"False"``."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def cascade_class_counts(cascade_ht_path: Path) -> dict:
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    ht = hl.read_table(str(cascade_ht_path))
    counts = dict(
        ht.group_by(ht.cascade_class)
        .aggregate(n=hl.agg.count())
        .to_pandas()
        .set_index("cascade_class")["n"]
    )
    return {str(k): int(v) for k, v in counts.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="a qtlcascade run's -o directory")
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--receipt", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    tissue_dir = out_dir / "per_tissue" / args.tissue
    cascade_ht = tissue_dir / "cascade.ht"
    gene_summary_tsv = tissue_dir / "gene_summary.tsv"
    coloc_tsv = tissue_dir / "coloc_results.tsv"
    report_html = out_dir / "qtlcascade_report.html"

    missing = [
        p
        for p in (cascade_ht / "_SUCCESS", gene_summary_tsv, report_html)
        if not p.exists()
    ]
    if missing:
        for p in missing:
            print(f"MISSING: {p}", file=sys.stderr)
        return 1

    class_counts = cascade_class_counts(cascade_ht)
    n_cascade_rows = sum(class_counts.values())

    gene_df = pd.read_csv(gene_summary_tsv, sep="\t")
    n_cascade_genes = len(gene_df)
    n_complete_cascade = (
        int(gene_df["has_complete_cascade"].map(_as_bool).sum())
        if "has_complete_cascade" in gene_df.columns
        else 0
    )
    gene_summary_sha256 = _sha256(gene_summary_tsv)

    coloc: dict = {"tested": None, "colocalised": None, "top5": None}
    coloc_results_sha256 = None
    if coloc_tsv.is_file():
        coloc_df = pd.read_csv(coloc_tsv, sep="\t")
        coloc["tested"] = len(coloc_df)
        coloc["colocalised"] = int((coloc_df["H4"] > H4_THRESHOLD).sum())
        top5 = coloc_df.sort_values("H4", ascending=False).head(5)
        coloc["top5"] = [
            {"gene_id": str(r["gene_id"]), "H4": float(r["H4"])}
            for _, r in top5.iterrows()
        ]
        coloc_results_sha256 = _sha256(coloc_tsv)

    summary = {
        "tissue": args.tissue,
        "cascade_class_counts": class_counts,
        "n_cascade_rows": n_cascade_rows,
        "n_cascade_genes": n_cascade_genes,
        "n_complete_cascade": n_complete_cascade,
        "coloc": coloc,
        "gene_summary_sha256": gene_summary_sha256,
        "coloc_results_sha256": coloc_results_sha256,
        "report_present": report_html.is_file(),
    }
    Path(args.receipt).write_text(json.dumps(summary, indent=2))

    print(f"tissue: {args.tissue}")
    print(f"cascade class counts: {class_counts} (sum {n_cascade_rows})")
    print(f"cascade genes: {n_cascade_genes}")
    print(f"genes with a complete cascade: {n_complete_cascade}")
    if coloc["tested"] is None:
        print("coloc: not run")
    else:
        print(f"coloc genes tested: {coloc['tested']}")
        print(f"colocalised (H4 > {H4_THRESHOLD}): {coloc['colocalised']}")
        print(f"top 5 genes by H4: {coloc['top5']}")
    print(f"gene_summary.tsv sha256: {gene_summary_sha256}")
    print(f"coloc_results.tsv sha256: {coloc_results_sha256}")
    print(f"report present: {summary['report_present']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
