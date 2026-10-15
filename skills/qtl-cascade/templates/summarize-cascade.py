"""Summarise one tissue's QTL cascade outputs under a `qtlcascade run`/`run_collection`
output directory: the cascade class counts, the cascade gene count, how many genes have
a complete cascade, how many genes were coloc-tested and colocalised (H4 > 0.8), the
top five genes by H4, the sha256 of `gene_summary.tsv` and `coloc_results.tsv`, and
whether the HTML report is present. Prints each value and writes the same as JSON to
`--receipt`. Runs inside the image: Hail to read `cascade.ht`, pandas for the TSVs.

    python summarize-cascade.py --out-dir OUT_DIR --tissue Liver --receipt OUT_DIR/summary.json \
        [--eqtl-ht E.ht --pqtl-ht P.ht --hgnc-ht H.ht]

The per-tissue directory is resolved with the pipeline's own tissue-name sanitiser
(`hvantk.algorithms.qtlcascade.pipeline._sanitize`): `qtlcascade run --tissues
<tissue>` writes under `per_tissue/<_sanitize(tissue)>`, but a run given no
`--tissues` filter at all writes under `per_tissue/all` instead (identical numbers,
for a tissue whose eQTL/pQTL tables are single-tissue to begin with) -- both are
accepted, trying the named one first.

Exits 1 (after printing which path is missing) when `cascade.ht` or
`gene_summary.tsv` is absent -- those always exist after any completed run. The
report is NOT required: a run with `--no-report` never writes one, so its absence
is recorded in the receipt (`report_present: false`), not an exit 1.
`coloc_results.tsv` is likewise optional: a cascade run without `--eqtl-allpairs`/
`--pqtl-allpairs` never writes it, so its absence is reported too, not an error.

`--eqtl-ht`/`--pqtl-ht`/`--hgnc-ht` are optional: when given, the receipt also
carries the sha256 of each table's `.provenance.json` sidecar (the same file
`dataset-build-at-scale`'s own verification ties a build to), so the summary can be
cross-checked against exactly which inputs produced it.
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


def resolve_tissue_dir(out_dir: Path, tissue: str) -> tuple[Path, str]:
    """The per-tissue directory ``qtlcascade run`` actually wrote to, and the
    label that run used for it. See the module docstring."""
    from hvantk.algorithms.qtlcascade.pipeline import _sanitize

    named = out_dir / "per_tissue" / _sanitize(tissue)
    if named.exists():
        return named, tissue
    return out_dir / "per_tissue" / "all", "all"


def cascade_class_counts(cascade_ht_path: Path) -> tuple[dict, int]:
    """Returns ``(class_counts, n_rows)``. ``n_rows`` is ``ht.count()`` -- the
    table's own row count -- rather than the sum of the grouped class counts, so a
    row whose class somehow failed to group (it never should; cascade.py's
    classification is total over every joined row) is never silently dropped from
    the total."""
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    ht = hl.read_table(str(cascade_ht_path))
    n_rows = ht.count()
    counts = dict(
        ht.group_by(ht.cascade_class)
        .aggregate(n=hl.agg.count())
        .to_pandas()
        .set_index("cascade_class")["n"]
    )
    return {str(k): int(v) for k, v in counts.items()}, n_rows


def _provenance_sha256(ht_path: str | None) -> str | None:
    if not ht_path:
        return None
    sidecar = Path(ht_path).with_name(Path(ht_path).name + ".provenance.json")
    return _sha256(sidecar) if sidecar.is_file() else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="a qtlcascade run's -o directory")
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--receipt", required=True)
    ap.add_argument(
        "--eqtl-ht", default=None, help="record its .provenance.json sha256"
    )
    ap.add_argument(
        "--pqtl-ht", default=None, help="record its .provenance.json sha256"
    )
    ap.add_argument(
        "--hgnc-ht", default=None, help="record its .provenance.json sha256"
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    tissue_dir, tissue_label = resolve_tissue_dir(out_dir, args.tissue)
    cascade_ht = tissue_dir / "cascade.ht"
    gene_summary_tsv = tissue_dir / "gene_summary.tsv"
    coloc_tsv = tissue_dir / "coloc_results.tsv"
    report_html = out_dir / "qtlcascade_report.html"

    missing = [p for p in (cascade_ht / "_SUCCESS", gene_summary_tsv) if not p.exists()]
    if missing:
        for p in missing:
            print(f"MISSING: {p}", file=sys.stderr)
        return 1

    class_counts, n_cascade_rows = cascade_class_counts(cascade_ht)

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

    report_present = report_html.is_file()

    summary = {
        "tissue": args.tissue,
        "tissue_label": tissue_label,
        "cascade_class_counts": class_counts,
        "n_cascade_rows": n_cascade_rows,
        "n_cascade_genes": n_cascade_genes,
        "n_complete_cascade": n_complete_cascade,
        "coloc": coloc,
        "gene_summary_sha256": gene_summary_sha256,
        "coloc_results_sha256": coloc_results_sha256,
        "report_present": report_present,
        "inputs": {
            "eqtl_ht_provenance_sha256": _provenance_sha256(args.eqtl_ht),
            "pqtl_ht_provenance_sha256": _provenance_sha256(args.pqtl_ht),
            "hgnc_ht_provenance_sha256": _provenance_sha256(args.hgnc_ht),
        },
    }
    Path(args.receipt).write_text(json.dumps(summary, indent=2))

    print(f"tissue: {args.tissue} (label: {tissue_label})")
    print(f"cascade class counts: {class_counts} (table rows: {n_cascade_rows})")
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
    print(f"report present: {report_present}")
    for key, value in summary["inputs"].items():
        if value is not None:
            print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
