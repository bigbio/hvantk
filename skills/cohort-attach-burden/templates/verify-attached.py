"""Post-hoc verification for a `cohort-attach-burden` run: reads the burden TSV, the
cohort manifest, the attached table and its coverage report, prints a summary, and
writes a receipt. Runs inside the image (Hail).

    python verify-attached.py --attached ATTACHED.ht --report ATTACHED.ht.report.json \
        --burden gene_burden.tsv --cohort cohort.yaml --receipt RECEIPT.json

Records: `n_genes`/`n_tested`/`mapping_rate`/`prior_non_null_rate` (from the attach
report), the burden table's own gene count and column list, its minimum `minp` and how
many rows have `p_adj < 0.05`, the attached table's key and which declared columns are
present, `burden_sha256` (the burden TSV's own file hash) and
`attached_content_sha256` (a digest over the *tested* rows' declared columns only, so
it is stable across reruns regardless of untested-gene ordering). Exits 1 if the
attached table has no `_SUCCESS` marker, its key is not `["gene_id"]`, a declared
cohort column is missing from it, or no gene was tested (`n_tested == 0`).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def read_burden(path: Path) -> list:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _norm(value):
    """Fixed-precision string for a float (stable across platforms); everything else
    passes through unchanged. Mirrors the shared digest convention for floats only --
    this receipt's content is plain scalars (gene_id plus the declared columns), never
    structs, locus or interval values."""
    return format(value, ".10g") if isinstance(value, float) else value


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--attached", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--burden", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--receipt", required=True)
    args = ap.parse_args()

    attached_path = Path(args.attached)
    report_path = Path(args.report)
    burden_path = Path(args.burden)

    problems = []
    if not (attached_path / "_SUCCESS").exists():
        problems.append(f"{attached_path} has no _SUCCESS marker")
    if not report_path.is_file():
        problems.append(f"report not found: {report_path}")
    if not burden_path.is_file():
        problems.append(f"burden table not found: {burden_path}")
    if problems:
        print("FAIL:")
        for p in problems:
            print(f"  - {p}")
        return 1

    from hvantk.algorithms.cohort.spec import load_cohort
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    manifest = load_cohort(args.cohort)
    report = json.loads(report_path.read_text())
    burden_rows = read_burden(burden_path)
    burden_columns = list(burden_rows[0].keys()) if burden_rows else []
    minp_values = [
        float(r["minp"]) for r in burden_rows if r.get("minp") not in (None, "")
    ]
    n_p_adj_sig = sum(
        1
        for r in burden_rows
        if r.get("p_adj") not in (None, "") and float(r["p_adj"]) < 0.05
    )

    ht = hl.read_table(str(attached_path))
    declared = list(manifest.declared_columns())
    missing_declared = [c for c in declared if c not in list(ht.row)]

    n_genes = ht.count()
    tested = ht.filter(ht.cohort_tested)
    n_tested = tested.count()
    content_strings = []
    if not missing_declared:
        for row in tested.select(*declared).collect():
            content_strings.append(
                json.dumps(
                    {k: _norm(v) for k, v in dict(row).items()},
                    sort_keys=True,
                    separators=(",", ":"),
                )
            )
    attached_content_sha256 = hashlib.sha256(
        "\n".join(sorted(content_strings)).encode()
    ).hexdigest()

    summary = {
        "n_genes": n_genes,
        "n_tested": n_tested,
        "mapping_rate": report.get("mapping_rate"),
        "prior_non_null_rate": report.get("prior", {}).get("non_null_rate"),
        "burden_n_genes": len(burden_rows),
        "burden_columns": burden_columns,
        "burden_min_minp": min(minp_values) if minp_values else None,
        "n_p_adj_sig": n_p_adj_sig,
        "attached_key": list(ht.key),
        "declared_columns_present": [c for c in declared if c not in missing_declared],
        "burden_sha256": hashlib.sha256(burden_path.read_bytes()).hexdigest(),
        "attached_content_sha256": attached_content_sha256,
    }
    print(json.dumps(summary, indent=2, default=str))

    fail = []
    if list(ht.key) != ["gene_id"]:
        fail.append(f"attached table key is {list(ht.key)}, expected ['gene_id']")
    if missing_declared:
        fail.append(
            f"declared column(s) missing from the attached table: {missing_declared}"
        )
    if n_tested == 0:
        fail.append("no gene was tested (n_tested == 0)")

    Path(args.receipt).write_text(
        json.dumps({**summary, "ok": not fail, "problems": fail}, indent=2, default=str)
    )

    if fail:
        print("FAIL:")
        for p in fail:
            print(f"  - {p}")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
