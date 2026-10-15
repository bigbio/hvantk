"""Verify built datasets against their plugin's committed schema CONTRACT (not exact
equality) and write a receipt. Run inside the container image (through
submit-hail-job.sh, or through run-in-image.sh for a parquet-only dataset with no
Hail involved): it needs Hail for ``.ht``/``.mt`` outputs, anndata for ``.h5ad``, and
pandas for ``.parquet``.

    python verify-builds.py [--list FILE] [--build DATASET=OUTPUT ...] [--receipt PATH]

FILE is the same TSV ``submit-build-array.sh`` reads
(``<provider:dataset>\\t<raw dir>\\t<output path>[\\t<extra options>]``); only the
first and third fields are used here. ``--build`` is repeatable and may be combined
with ``--list``, for verifying a build that was not made from a task list.

For each dataset, resolves its spec from the plugin registry, reads the output by
artifact type (``hl.read_table`` for ``.ht``, ``hl.read_matrix_table`` for ``.mt``,
``anndata.read_h5ad`` for ``.h5ad``, ``pandas.read_parquet`` for ``.parquet``), and
checks it against the plugin's committed ``schema_snapshot`` as a CONTRACT, not exact
equality: a snapshot is fixture-derived, and a full build routinely carries more
columns than its fixture did (GenCC: 31 upstream columns, 13 in its fixture) while
still satisfying the plugin's contract. The contract passes when every field the
snapshot declares is present in the built schema with the same type; extra built
fields are allowed and listed (``extra_fields``). The committed ``row_snapshot``'s keys
are checked too, but only informationally (``snapshot_keys_present``): a hand-authored
fixture (ClinGen's invents its MONDO ids) can declare sample keys that do not exist in
a real build, so a missing one is reported, never a failure.

One dataset's exception never loses the others' results: each is wrapped in its own
try/except, and a receipt is written for whatever was verified. Prints one line per
dataset and writes the receipt JSON, keyed ``"<dataset>=<output>"`` (so verifying the
same dataset's output twice, e.g. before and after a rebuild, does not collide). Exits
1 if any dataset's output is missing, has no ``_SUCCESS``/sidecar, fails its schema
contract, or has a sidecar ``schema_id`` that does not match the registry's.
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


# --- schema contract (floor, not equality) -----------------------------------------


def _fields_contract(
    built: dict, snapshot: dict, prefix: str = ""
) -> tuple[list[str], dict, list[str]]:
    """(missing_fields, type_mismatches, extra_fields) for one field/type-dict level.

    A snapshot value that is itself a dict (a field whose type is a bare struct) is
    checked recursively with a dotted path name; any other value is compared for
    exact string equality (Hail/AnnData/pandas type labels already render as stable
    strings, so string equality is the right notion of "same type").
    """
    missing: list[str] = []
    mismatches: dict = {}
    for name, expected in snapshot.items():
        path = f"{prefix}{name}"
        if name not in built:
            missing.append(path)
            continue
        got = built[name]
        if isinstance(expected, dict):
            if not isinstance(got, dict):
                mismatches[path] = {"expected": expected, "got": got}
                continue
            sub_missing, sub_mismatches, _ = _fields_contract(
                got, expected, prefix=f"{path}."
            )
            missing.extend(sub_missing)
            mismatches.update(sub_mismatches)
        elif got != expected:
            mismatches[path] = {"expected": expected, "got": got}
    extra = sorted(set(built) - set(snapshot))
    return missing, mismatches, extra


def schema_contract(built: dict, snapshot: dict) -> dict:
    """Does *built* satisfy *snapshot* as a floor? Extra built fields are allowed.

    Dispatches on the snapshot's own shape: a Hail Table (``{"key", "row"}``), a Hail
    MatrixTable (``{"row_key", "col_key", "row", "col", "entry"}``), an AnnData
    (``{"obs_columns", "var_columns", "X_format", ...}`` -- ``n_obs``, ``n_vars``,
    ``X_dtype`` and ``layers`` are ignored), or a parquet-backed pandas table
    (``{"columns", "dtypes"}``). Returns
    ``{"ok", "missing_fields", "type_mismatches", "extra_fields"}``.
    """
    if "obs_columns" in snapshot:
        missing = sorted(
            [
                c
                for c in snapshot.get("obs_columns", [])
                if c not in built.get("obs_columns", [])
            ]
            + [
                c
                for c in snapshot.get("var_columns", [])
                if c not in built.get("var_columns", [])
            ]
        )
        mismatches = {}
        if built.get("X_format") != snapshot.get("X_format"):
            mismatches["X_format"] = {
                "expected": snapshot.get("X_format"),
                "got": built.get("X_format"),
            }
        extra = sorted(
            (set(built.get("obs_columns", [])) - set(snapshot.get("obs_columns", [])))
            | (set(built.get("var_columns", [])) - set(snapshot.get("var_columns", [])))
        )
        return {
            "ok": not missing and not mismatches,
            "missing_fields": missing,
            "type_mismatches": mismatches,
            "extra_fields": extra,
        }
    if "columns" in snapshot and "dtypes" in snapshot:
        missing = [c for c in snapshot["columns"] if c not in built.get("columns", [])]
        mismatches = {}
        for c in snapshot["columns"]:
            if c in built.get("columns", []):
                exp = snapshot.get("dtypes", {}).get(c)
                got = built.get("dtypes", {}).get(c)
                if got != exp:
                    mismatches[c] = {"expected": exp, "got": got}
        extra = sorted(set(built.get("columns", [])) - set(snapshot.get("columns", [])))
        return {
            "ok": not missing and not mismatches,
            "missing_fields": sorted(missing),
            "type_mismatches": mismatches,
            "extra_fields": extra,
        }
    if "row_key" in snapshot:
        missing = []
        mismatches = {}
        key_ok = built.get("row_key") == snapshot.get("row_key") and built.get(
            "col_key"
        ) == snapshot.get("col_key")
        if built.get("row_key") != snapshot.get("row_key"):
            mismatches["row_key"] = {
                "expected": snapshot.get("row_key"),
                "got": built.get("row_key"),
            }
        if built.get("col_key") != snapshot.get("col_key"):
            mismatches["col_key"] = {
                "expected": snapshot.get("col_key"),
                "got": built.get("col_key"),
            }
        extra = []
        for part in ("row", "col", "entry"):
            m, t, e = _fields_contract(
                built.get(part, {}), snapshot.get(part, {}), prefix=f"{part}."
            )
            missing.extend(m)
            mismatches.update(t)
            extra.extend(f"{part}.{x}" for x in e)
        ok = key_ok and not missing and not mismatches
        return {
            "ok": ok,
            "missing_fields": sorted(missing),
            "type_mismatches": mismatches,
            "extra_fields": sorted(extra),
        }
    # Hail Table: {"key": [...], "row": {...}}
    key_ok = built.get("key") == snapshot.get("key")
    mismatches = {}
    if not key_ok:
        mismatches["key"] = {"expected": snapshot.get("key"), "got": built.get("key")}
    missing, row_mismatches, extra = _fields_contract(
        built.get("row", {}), snapshot.get("row", {})
    )
    mismatches.update(row_mismatches)
    ok = key_ok and not missing and not mismatches
    return {
        "ok": ok,
        "missing_fields": sorted(missing),
        "type_mismatches": mismatches,
        "extra_fields": sorted(extra),
    }


# --- reading a built output ---------------------------------------------------------


def read_artifact(path: Path):
    """Return (schema_dict, row_count, keyed_table) for one built dataset.

    ``keyed_table`` is the Hail Table to semi-join sample keys against (the
    MatrixTable's own row table for a ``.mt``), or ``None`` when the format has no
    such concept here (AnnData, a flat pandas table) -- the informational row-snapshot
    check is skipped for those.
    """
    if path.suffix in (".ht", ".mt"):
        from hvantk.core.utils.hail_context import init_hail

        init_hail(quiet=True)
        import hail as hl

        from hvantk.tests._snapshot_utils import hail_schema_to_dict

        if path.suffix == ".ht":
            ht = hl.read_table(str(path))
            return hail_schema_to_dict(ht), ht.count(), ht
        mt = hl.read_matrix_table(str(path))
        return hail_schema_to_dict(mt), mt.count_rows(), mt.rows()
    if path.suffix == ".h5ad":
        import anndata as ad

        from hvantk.tests._snapshot_utils import anndata_schema_to_dict

        adata = ad.read_h5ad(str(path))
        return anndata_schema_to_dict(adata), int(adata.n_obs), None
    if path.suffix == ".parquet":
        import pandas as pd

        from hvantk.core.models.annotation_table import _normalize_dtype

        df = pd.read_parquet(path)
        schema = {
            "n_rows": int(len(df)),
            "columns": list(df.columns),
            "dtypes": {c: _normalize_dtype(t) for c, t in df.dtypes.items()},
        }
        return schema, int(len(df)), None
    raise ValueError(f"no verifier for output extension {path.suffix!r}: {path}")


def _declared_keys(
    row_snapshot_entries: list, key_field_names: list[str]
) -> list[dict]:
    """Key dicts from a committed row_snapshot.

    Most plugins' ``sample_rows.json`` is ``[{"key": {...}, "row": {...}}, ...]``
    (``collect_sample_rows``'s own shape). A few (``onek-genomes``) commit a flat list
    of row dicts with no ``"key"`` wrapper; for those, the key fields are picked out
    by name from the schema snapshot's own ``"key"`` list.
    """
    keys = []
    for entry in row_snapshot_entries:
        key = entry.get("key")
        if isinstance(key, dict):
            keys.append(key)
        else:
            keys.append(
                {name: entry[name] for name in key_field_names if name in entry}
            )
    return keys


def _present_key_count(
    table, keys: list[dict], key_field_names: list[str]
) -> tuple[int, int]:
    """(keys actually present in *table*, total keys declared), via semi_join.

    Never raises on a missing key, unlike ``collect_sample_rows``: the keys a
    committed row-snapshot declares may not exist in a build from the real upstream
    file (ClinGen's fixture invents its MONDO ids), so this is informational only.
    """
    if not keys or table is None:
        return 0, len(keys)
    import hail as hl

    from hvantk.tests._snapshot_utils import _jsonable_to_hail_python

    key_dtype = table.key.dtype
    key_rows = [
        {
            name: _jsonable_to_hail_python(k.get(name), key_dtype[name])
            for name in key_field_names
        }
        for k in keys
    ]
    requested = hl.Table.parallelize(key_rows, schema=key_dtype, key=key_field_names)
    found = table.semi_join(requested).count()
    return found, len(keys)


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
    details["schema_id_matches"] = details["schema_id"] == spec.schema_id

    schema, row_count, keyed_table = read_artifact(out)
    schema_snapshot = Path(spec.test_paths.schema_snapshot)
    expected_schema = json.loads(schema_snapshot.read_text())
    contract = schema_contract(schema, expected_schema)
    details["schema_contract_ok"] = contract["ok"]
    details["missing_fields"] = contract["missing_fields"]
    details["type_mismatches"] = contract["type_mismatches"]
    details["extra_fields"] = contract["extra_fields"]
    details["row_count"] = row_count

    row_snapshot = Path(spec.test_paths.row_snapshot)
    if keyed_table is not None and row_snapshot.is_file() and "key" in expected_schema:
        declared = _declared_keys(
            json.loads(row_snapshot.read_text()), expected_schema["key"]
        )
        found, total = _present_key_count(keyed_table, declared, expected_schema["key"])
        details["snapshot_keys_present"] = {"found": found, "total": total}
    else:
        details["snapshot_keys_present"] = None
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
        key = f"{dataset}={output}"
        try:
            details = verify_one(dataset, output)
        except Exception as exc:  # noqa: BLE001 -- one dataset's exception must not
            # cost the receipt for the others.
            details = {"dataset": dataset, "output": output, "error": str(exc)}
        receipt["datasets"][key] = details
        passed = (
            "error" not in details
            and details.get("schema_contract_ok") is True
            and details.get("schema_id_matches") is True
        )
        ok = ok and passed
        status = "OK" if passed else "FAIL"
        print(
            f"{status}  {dataset:<24} rows={details.get('row_count', '?')!s:<10} "
            f"schema_id={details.get('schema_id')} "
            f"schema_id_matches={details.get('schema_id_matches')} "
            f"schema_contract_ok={details.get('schema_contract_ok')} "
            f"extra_fields={details.get('extra_fields')} "
            f"snapshot_keys_present={details.get('snapshot_keys_present')}"
            + (f" error={details['error']}" if "error" in details else "")
        )

    if args.receipt:
        Path(args.receipt).write_text(json.dumps(receipt, indent=2, default=str))
        print(f"receipt: {args.receipt}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
