"""Describe a built Hail table or matrix table and write a receipt.

Run inside the container (through submit-hail-job.sh). Prints the key, the row count
and the field names, and writes a receipt JSON tying those numbers to the SHA-256 of
the output's ``.provenance.json`` sidecar.

    python hail-describe.py <path.ht|path.mt> [--receipt out.json]
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--receipt", default=None)
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    path = Path(args.path)
    is_mt = path.suffix == ".mt"
    obj = hl.read_matrix_table(str(path)) if is_mt else hl.read_table(str(path))
    key = list(obj.row_key) if is_mt else list(obj.key)
    n_rows = obj.count_rows() if is_mt else obj.count()
    fields = list(obj.row)
    sidecar = path.with_name(path.name + ".provenance.json")
    prov_sha = (
        hashlib.sha256(sidecar.read_bytes()).hexdigest() if sidecar.exists() else None
    )
    schema_id = (
        json.loads(sidecar.read_text()).get("schema_id") if sidecar.exists() else None
    )

    print(f"path:      {path}")
    print(f"key:       {key}")
    print(f"rows:      {n_rows}")
    print(f"fields:    {fields}")
    print(f"schema_id: {schema_id}")
    print(f"sidecar:   {sidecar if sidecar.exists() else 'MISSING'}")

    if args.receipt:
        receipt = {
            "output": str(path),
            "rows": n_rows,
            "key": key,
            "fields": fields,
            "schema_id": schema_id,
            "provenance_sha256": prov_sha,
            "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
        Path(args.receipt).write_text(json.dumps(receipt, indent=2))
        print(f"receipt:   {args.receipt}")


if __name__ == "__main__":
    main()
