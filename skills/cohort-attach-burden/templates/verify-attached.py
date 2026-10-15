"""Post-hoc verification for a `cohort-attach-burden` run: reads the burden/gene
table, the cohort manifest, the attached table and its coverage report, prints a
summary, and writes a receipt -- in every exit path, including a failing one, so the
receipt always reflects whatever could be read. Runs inside the image (Hail).

    python verify-attached.py --attached ATTACHED.ht --report ATTACHED.ht.report.json \
        --burden gene_burden.tsv --cohort cohort.yaml --layer1 LAYER1.ht \
        --receipt RECEIPT.json

Records: `n_genes`/`n_tested`/`mapping_rate`/`prior_non_null_rate` (from the attach
report), the burden/gene table's own row count and column list, the manifest's own
prior column and its minimum value, how many rows have `p_adj < 0.05` (burden-specific;
null when the column is absent, e.g. for an externally-supplied table), the attached
table's key, which declared columns are present with their Hail type (and which of
those types are non-numeric -- a declared column such as `conc`/`driver_af` typed
`str` instead of `float64` is exactly the symptom of `cohort attach` importing a
cohort table with a blank cell that `hl.import_table`'s default `missing` list
(`["NA"]`) does not recognise as missing, in hvantk 0.5.0 and earlier -- see
`SKILL.md`'s Traps), `--layer1`'s own path (for traceability -- nothing the attached
table carries names which Layer-1 input it came from), `burden_sha256` (the burden/gene
table's own file hash) and `attached_content_sha256` (a digest over the *tested* rows'
declared columns only, so it is stable across reruns regardless of untested-gene
ordering).

Exits 1 if the attached table has no `_SUCCESS` marker, its key is not `["gene_id"]`,
a declared cohort column is missing from it, or no gene was tested (`n_tested == 0`).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def read_table_rows(path: Path) -> list:
    """Read a delimited gene/burden table with the same delimiter-sniffing rule
    `hvantk cohort validate`/`attach` use (a literal tab in the header line means TSV,
    else CSV; see `hvantk.algorithms.cohort.checks.detect_delimiter`), so this script
    never disagrees with them about how to parse the same file."""
    from hvantk.algorithms.cohort.checks import detect_delimiter

    delimiter = detect_delimiter(str(path))
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def _norm(value):
    """Fixed-precision string for a float (stable across platforms); everything else
    passes through unchanged. Mirrors the shared digest convention for floats only --
    this receipt's content is plain scalars (gene_id plus the declared columns), never
    structs, locus or interval values."""
    return format(value, ".10g") if isinstance(value, float) else value


def write_receipt(receipt_path: Path, summary: dict, fail: list) -> None:
    Path(receipt_path).write_text(
        json.dumps({**summary, "ok": not fail, "problems": fail}, indent=2, default=str)
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--attached", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--burden", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument(
        "--layer1",
        required=True,
        help="path to the Layer-1 table attached onto; recorded in the receipt only "
        "(the attached table itself carries no record of which Layer-1 input it came "
        "from)",
    )
    ap.add_argument("--receipt", required=True)
    args = ap.parse_args()

    attached_path = Path(args.attached)
    report_path = Path(args.report)
    burden_path = Path(args.burden)
    summary: dict = {"layer1": args.layer1}

    fail = []
    if not (attached_path / "_SUCCESS").exists():
        fail.append(f"{attached_path} has no _SUCCESS marker")
    if not report_path.is_file():
        fail.append(f"report not found: {report_path}")
    if not burden_path.is_file():
        fail.append(f"burden/gene table not found: {burden_path}")
    if fail:
        write_receipt(args.receipt, summary, fail)
        print("FAIL:")
        for p in fail:
            print(f"  - {p}")
        return 1

    from hvantk.algorithms.cohort.spec import load_cohort
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    manifest = load_cohort(args.cohort)
    report = json.loads(report_path.read_text())
    burden_rows = read_table_rows(burden_path)
    burden_columns = list(burden_rows[0].keys()) if burden_rows else []

    # The manifest's OWN declared prior column, not a hardcoded "minp": this script
    # also covers the second entry point (an externally-supplied gene table whose
    # prior is named whatever that cohort calls it, e.g. "p_dee").
    prior_col = manifest.prior.column
    prior_values = [
        float(r[prior_col])
        for r in burden_rows
        if prior_col in r and r.get(prior_col) not in (None, "")
    ]
    # p_adj is burden-specific (hvantk cohort burden --mtc); an externally-supplied
    # table need not have it, so this is None rather than 0 when the column is absent.
    if any("p_adj" in r for r in burden_rows):
        n_p_adj_sig = sum(
            1
            for r in burden_rows
            if r.get("p_adj") not in (None, "") and float(r["p_adj"]) < 0.05
        )
    else:
        n_p_adj_sig = None

    ht = hl.read_table(str(attached_path))
    declared = list(manifest.declared_columns())
    missing_declared = [c for c in declared if c not in list(ht.row)]
    present_declared = [c for c in declared if c not in missing_declared]
    declared_column_types = {c: str(ht[c].dtype) for c in present_declared}
    non_numeric_declared_columns = [
        c for c in present_declared if not hl.is_numeric(ht[c].dtype)
    ]

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

    summary.update(
        {
            "n_genes": n_genes,
            "n_tested": n_tested,
            "mapping_rate": report.get("mapping_rate"),
            "prior_non_null_rate": report.get("prior", {}).get("non_null_rate"),
            "burden_n_genes": len(burden_rows),
            "burden_columns": burden_columns,
            "prior_column": prior_col,
            "burden_min_prior": min(prior_values) if prior_values else None,
            "n_p_adj_sig": n_p_adj_sig,
            "attached_key": list(ht.key),
            "declared_columns_present": present_declared,
            "declared_column_types": declared_column_types,
            "non_numeric_declared_columns": non_numeric_declared_columns,
            "burden_sha256": hashlib.sha256(burden_path.read_bytes()).hexdigest(),
            "attached_content_sha256": attached_content_sha256,
        }
    )
    print(json.dumps(summary, indent=2, default=str))
    if non_numeric_declared_columns:
        print(
            "WARNING: declared column(s) with a non-numeric type (expected numeric): "
            f"{non_numeric_declared_columns}"
        )

    if list(ht.key) != ["gene_id"]:
        fail.append(f"attached table key is {list(ht.key)}, expected ['gene_id']")
    if missing_declared:
        fail.append(
            f"declared column(s) missing from the attached table: {missing_declared}"
        )
    if n_tested == 0:
        fail.append("no gene was tested (n_tested == 0)")

    write_receipt(args.receipt, summary, fail)

    if fail:
        print("FAIL:")
        for p in fail:
            print(f"  - {p}")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
