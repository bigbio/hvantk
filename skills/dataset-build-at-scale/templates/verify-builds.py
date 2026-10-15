"""Verify built datasets against their plugin's committed schema and sample-row
snapshots, and write a receipt. Run inside the container image (through
submit-hail-job.sh or as an array-task-free follow-up job): it needs Hail for
``.ht``/``.mt`` outputs and anndata for ``.h5ad`` ones.

    python verify-builds.py [--list FILE] [--build DATASET=OUTPUT ...] [--receipt PATH]

FILE is the same TSV ``submit-build-array.sh`` reads
(``<provider:dataset>\\t<raw dir>\\t<output path>[\\t<extra options>]``); only the
first and third fields are used here. ``--build`` is repeatable and may be combined
with ``--list``, for verifying a build that was not made from a task list.

For each dataset, resolves its spec from the plugin registry, reads the output by
artifact type (``hl.read_table`` for ``.ht``, ``hl.read_matrix_table`` for ``.mt``,
``anndata.read_h5ad`` for ``.h5ad``), compares its schema with the plugin's committed
``schema_snapshot``, checks that every key in the committed ``row_snapshot`` is present
in the output, counts rows, and hashes the ``.provenance.json`` sidecar. Prints one
line per dataset and writes the receipt JSON. Exits 1 if any dataset is missing, has no
``_SUCCESS``, fails its committed schema, or is missing one of its committed sample-row
keys.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path


def pairs_from_list(path: str) -> list[tuple[str, str]]:
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError(f"{path}: line is missing fields: {line!r}")
        out.append((fields[0], fields[2]))
    return out


def pairs_from_build(specs: list[str]) -> list[tuple[str, str]]:
    out = []
    for item in specs:
        if "=" not in item:
            raise ValueError(f"--build must be DATASET=OUTPUT, got: {item!r}")
        ds, out_path = item.split("=", 1)
        out.append((ds, out_path))
    return out


def read_artifact(path: Path):
    """Return (schema_dict, row_count, collect_fn) for one built dataset.

    ``collect_fn(keys)`` collects sample rows for the given committed keys, or is
    ``None`` when the output's format has no sample-row collector (yet).
    """
    if path.suffix == ".ht":
        import hail as hl

        from hvantk.tests._snapshot_utils import (
            collect_sample_rows,
            hail_schema_to_dict,
        )

        ht = hl.read_table(str(path))
        return (
            hail_schema_to_dict(ht),
            ht.count(),
            lambda keys: collect_sample_rows(ht, keys),
        )
    if path.suffix == ".mt":
        import hail as hl

        from hvantk.tests._snapshot_utils import (
            collect_sample_rows,
            hail_schema_to_dict,
        )

        mt = hl.read_matrix_table(str(path))
        # The row table shares the MatrixTable's row key, which is what a VariantMatrix
        # snapshot's keys are expressed in.
        return (
            hail_schema_to_dict(mt),
            mt.count_rows(),
            lambda keys: collect_sample_rows(mt.rows(), keys),
        )
    if path.suffix == ".h5ad":
        import anndata as ad

        from hvantk.tests._snapshot_utils import anndata_schema_to_dict

        adata = ad.read_h5ad(str(path))
        return anndata_schema_to_dict(adata), int(adata.n_obs), None
    raise ValueError(f"no verifier for output extension {path.suffix!r}: {path}")


def verify_one(dataset: str, output: str) -> dict:
    from hvantk.core.plugin.loader import get_registry

    out = Path(output)
    details: dict = {"dataset": dataset, "output": str(out)}
    if not out.exists():
        details["error"] = "output not found"
        return details
    details["success_marker"] = (
        (out / "_SUCCESS").exists() if out.is_dir() else out.is_file()
    )
    sidecar = out.with_name(out.name + ".provenance.json")
    details["sidecar"] = sidecar.is_file()
    if not details["success_marker"] or not details["sidecar"]:
        details["error"] = "missing _SUCCESS or provenance sidecar"
        return details
    prov = json.loads(sidecar.read_text())
    details["schema_id"] = prov.get("schema_id")
    details["provenance_sha256"] = hashlib.sha256(sidecar.read_bytes()).hexdigest()

    spec = get_registry().get_dataset(dataset)
    schema, row_count, collect_fn = read_artifact(out)
    schema_snapshot = Path(spec.test_paths.schema_snapshot)
    expected_schema = json.loads(schema_snapshot.read_text())
    details["schema_matches"] = schema == expected_schema
    details["row_count"] = row_count

    row_snapshot = Path(spec.test_paths.row_snapshot)
    sample_rows: list = []
    sample_keys_present = True
    if collect_fn is not None and row_snapshot.is_file():
        keys = [row["key"] for row in json.loads(row_snapshot.read_text())]
        try:
            sample_rows = collect_fn(keys) if keys else []
        except KeyError:
            sample_keys_present = False
    details["sample_keys_present"] = sample_keys_present
    details["sample_rows"] = sample_rows
    return details


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", dest="list_path", default=None, help="task-list TSV")
    ap.add_argument(
        "--build",
        dest="builds",
        action="append",
        default=[],
        metavar="DATASET=OUTPUT",
        help="verify one dataset's output directly; repeatable",
    )
    ap.add_argument("--receipt", default=None, help="write the receipt JSON here")
    args = ap.parse_args()

    pairs = []
    if args.list_path:
        pairs.extend(pairs_from_list(args.list_path))
    pairs.extend(pairs_from_build(args.builds))
    if not pairs:
        ap.error("give --list, --build, or both")

    receipt = {
        "datasets": {},
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    ok = True
    for dataset, output in pairs:
        details = verify_one(dataset, output)
        receipt["datasets"][dataset] = details
        passed = (
            "error" not in details
            and details.get("schema_matches") is True
            and details.get("sample_keys_present") is True
        )
        ok = ok and passed
        status = "OK" if passed else "FAIL"
        print(
            f"{status}  {dataset:<24} rows={details.get('row_count', '?')!s:<10} "
            f"schema_id={details.get('schema_id')} "
            f"schema_matches={details.get('schema_matches')} "
            f"sample_keys_present={details.get('sample_keys_present')}"
            + (f" error={details['error']}" if "error" in details else "")
        )

    if args.receipt:
        Path(args.receipt).write_text(json.dumps(receipt, indent=2, default=str))
        print(f"receipt: {args.receipt}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
