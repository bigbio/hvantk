"""Deterministic check for the dataset-build-at-scale eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json``, and exits 0 on pass. Runs inside the
container on a compute node (it needs Hail). Slurm accounting is read from a file
captured on the login node (``--sacct-file``), because the container has no Slurm
client; without a file it tries ``sacct`` and records M2 and M3 as null when that is
unavailable or fails.

M4 checks each output against its plugin's committed schema snapshot as a CONTRACT
(``schema_contract``, read from THIS checkout's own ``REPO_ROOT``, never from
``reference.json`` or whatever checkout produced the output): a snapshot is
fixture-derived and a full build may carry more columns than its fixture did (GenCC:
31 upstream columns, 13 in the fixture) while still satisfying the plugin's contract.
M5's ``content_ok`` is the stricter, exact check that this run REPRODUCES the
reference build: schema equal to the reference's stored schema, row counts equal, and
``rows_digest`` equal over a fixed, deterministically-chosen sample of keys -- never a
plugin's own ``SAMPLE_KEYS`` or committed ``sample_rows.json``, whose keys (ClinGen's
fixture invents its MONDO ids) are not guaranteed to exist in a real build.

Array-aware, unlike a single-job eval's check.py: the task asks for one array job
covering all three datasets. M3's "every listed job exists in accounting" allows a
notes id ``A`` to be satisfied by an accounting row ``A`` *or* ``A_<task>`` (an array
submission writes one ``job: <A>`` line; every finished task is only ever accounted
under ``<A>_<task>``, never bare ``A``). ``array_used`` is None when there is no
accounting, or no accounting row at all under this run's name prefix; otherwise it is
true only when some one array parent completed at least 2 of its tasks AND the
prefixed array-shaped tasks completed at least 3 in total -- true for one 3-task array
(or a 3-task array with one failure resubmitted as a 1-task array), false for three
separate one-task "arrays" or three plain jobs. ``array_parents`` is kept for
diagnostics.

    python check.py --out-dir DIR --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--verdict-dir DIR] [--expect-fail]
    python check.py --out-dir DIR --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

DIR holds the three fixed-name outputs (see DATASETS below). The accounting file is
the output of
``sacct -X -n -P -S <since> -o JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 30}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)
N_SAMPLE_KEYS = 12

#: Fixed output names under --out-dir, one per dataset the eval builds.
#: ``schema_snapshot`` is the plugin's own committed contract, read fresh from THIS
#: checkout every time -- never copied by value, never read from reference.json.
DATASETS = {
    "hgnc:lookup": {
        "output": "hgnc-lookup.ht",
        "schema_id": "hgnc-lookup-v1",
        "schema_snapshot": REPO_ROOT / "hvantk/skills/hgnc/tests/snapshots/schema.json",
    },
    "clingen:gene-disease": {
        "output": "clingen-gene-disease.ht",
        "schema_id": "clingen-gene-disease-v1",
        "schema_snapshot": REPO_ROOT
        / "hvantk/skills/clingen/tests/snapshots/schema.json",
    },
    "gencc:submissions": {
        "output": "gencc-submissions.ht",
        "schema_id": "gencc-submissions-v1",
        "schema_snapshot": REPO_ROOT
        / "hvantk/skills/gencc/tests/snapshots/schema.json",
    },
}

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")
_ARRAY_TASK_ID = re.compile(r"^\d+_\d+$")


# --- schema contract (floor, not equality); see templates/verify-builds.py for the
# --- twin copy this is deliberately kept in sync with by hand (each eval/check.py in
# --- this repository is self-contained; see e.g. hpc-hail-job's). ---------------------


def _fields_contract(
    built: dict, snapshot: dict, prefix: str = ""
) -> tuple[list[str], dict, list[str]]:
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

    All three eval datasets are Hail Tables (``{"key", "row"}``); the MatrixTable,
    AnnData and parquet branches are kept for parity with ``verify-builds.py`` and are
    unreachable here, but a future fourth eval dataset should not need them rewritten.
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
    # Hail Table: {"key": [...], "row": {...}} -- every eval dataset's own shape.
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


# --- deterministic sample keys and the order-insensitive row digest ----------------


def choose_sample_keys(ht, n_keys: int = N_SAMPLE_KEYS) -> list[dict]:
    """``n_keys`` keys chosen from *ht* by even position, not by any plugin's
    hand-picked (and possibly fixture-only) SAMPLE_KEYS."""
    from hvantk.tests._snapshot_utils import _to_jsonable

    import hail as hl

    n = ht.count()
    if n == 0:
        return []
    positions = sorted({(i * n) // n_keys for i in range(n_keys)})
    key_field_names = list(ht.key)
    indexed = ht.add_index("_i")
    chosen = indexed.filter(hl.literal(positions).contains(indexed["_i"])).select(
        *key_field_names
    )
    keys = [
        {name: _to_jsonable(row[name]) for name in key_field_names}
        for row in chosen.collect()
    ]
    keys.sort(key=lambda k: json.dumps(k, sort_keys=True))
    return keys


def rows_digest(ht, keys: list[dict]) -> str:
    """SHA-256 over every row of *ht* whose key is in *keys*, order-insensitive and
    duplicate-key-safe. A key absent from *ht* is simply absent from the result (via
    semi_join, which never raises) rather than a crash: a real content difference then
    shows up as a digest mismatch, which is what M5 is for."""
    if not keys:
        return hashlib.sha256(b"").hexdigest()
    from hvantk.tests._snapshot_utils import _jsonable_to_hail_python, _to_jsonable

    import hail as hl

    key_field_names = list(ht.key)
    key_dtype = ht.key.dtype
    key_rows = [
        {
            name: _jsonable_to_hail_python(k.get(name), key_dtype[name])
            for name in key_field_names
        }
        for k in keys
    ]
    requested = hl.Table.parallelize(key_rows, schema=key_dtype, key=key_field_names)
    collected = ht.semi_join(requested).collect()
    strings = sorted(
        json.dumps(_to_jsonable(dict(row)), sort_keys=True, separators=(",", ":"))
        for row in collected
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


# --- run notes and Slurm accounting (M1-M3) -----------------------------------------


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {"present": False, "commit": None, "container": None, "jobs": []}
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
    }


def mem_to_gb(s: str, cpus: int) -> float | None:
    """ReqMem as GB; ``c`` means per core, ``n`` per node. Binary units, as Slurm."""
    m = _MEM.match(s.strip())
    if not m:
        return None
    unit = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1.0, "T": 1024.0}[m.group(2)]
    n = float(m.group(1)) * unit
    return n * cpus if m.group(3) == "c" else n


def limit_to_minutes(s: str) -> float | None:
    """Timelimit as minutes; None for UNLIMITED, Partition_Limit or anything unparsed."""
    s = s.strip()
    days = 0
    try:
        if "-" in s:
            d, s = s.split("-", 1)
            days = int(d)
        parts = [int(p) for p in s.split(":")]
    except ValueError:
        return None
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h, (m, sec) = 0, parts
    elif len(parts) == 1:
        h, m, sec = 0, parts[0], 0
    else:
        return None
    return days * 1440 + h * 60 + m + sec / 60


def allowed_partitions(profile: Path | None) -> set[str]:
    if profile is None or not profile.is_file():
        return set()
    out = set()
    for line in profile.read_text().splitlines():
        for key in _PROFILE_PARTITIONS:
            m = re.match(rf"^\s*(?:export\s+)?{key}=[\"']?([^\"'#\s]+)", line)
            if m:
                out.add(m.group(1))
    return out


def parse_sacct(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        f = line.split("|")
        if len(f) < 7:
            continue
        rows.append(
            dict(zip(["id", "name", "state", "mem", "cpus", "limit", "partition"], f))
        )
    return rows


def load_sacct(sacct_file: str | None, since: str) -> tuple[list[dict] | None, str]:
    """Accounting rows, or (None, reason) when accounting is not available here."""
    if sacct_file:
        p = Path(sacct_file)
        if not p.is_file():
            return None, f"sacct file not found: {sacct_file}"
        return parse_sacct(p.read_text()), ""
    cmd = ["sacct", "-X", "-n", "-P", "-S", since, "-o", SACCT_FIELDS]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        return None, f"sacct unavailable: {exc}"
    return parse_sacct(out), ""


def job_known(job_id: str, all_ids: set[str]) -> bool:
    """True if *job_id* (from run-notes) is ``A`` or starts an array task ``A_<i>``."""
    return job_id in all_ids or any(i.startswith(job_id + "_") for i in all_ids)


def grade_slurm(
    run_id: str, notes: dict, rows: list[dict], profile: Path | None
) -> tuple[bool, bool, bool | None, list[str], dict]:
    jobs = [j for j in rows if j["name"].startswith(run_id)]
    allowed = allowed_partitions(profile)
    within = []
    for j in jobs:
        try:
            cpus = int(j["cpus"] or 0)
        except ValueError:
            cpus = 0
        gb = mem_to_gb(j["mem"], cpus) if j["mem"] else None
        mins = limit_to_minutes(j["limit"]) if j["limit"] else None
        ok = (
            cpus <= BUDGET["max_cpus"]
            and gb is not None
            and gb <= BUDGET["max_mem_gb"] + 0.5
            and mins is not None
            and mins <= BUDGET["max_minutes"] + 0.01
            and (not allowed or j["partition"] in allowed)
        )
        within.append(ok)
    m2 = bool(jobs) and all(within)
    all_ids = {j["id"] for j in rows}
    m3 = (
        bool(notes["jobs"])
        and any(j["state"].startswith(SUCCESS_STATES) for j in jobs)
        and all(job_known(n, all_ids) for n in notes["jobs"])
    )

    array_rows = [j for j in jobs if _ARRAY_TASK_ID.match(j["id"])]
    array_parents = sorted({j["id"].split("_", 1)[0] for j in array_rows})
    if not jobs:
        # No accounting row at all under this run's name prefix: we cannot tell an
        # array run from anything else, so this is unknown, not "not an array".
        array_used: bool | None = None
    else:
        completed_per_parent: dict[str, int] = {}
        for j in array_rows:
            if j["state"] == "COMPLETED":
                parent = j["id"].split("_", 1)[0]
                completed_per_parent[parent] = completed_per_parent.get(parent, 0) + 1
        some_parent_has_2 = any(c >= 2 for c in completed_per_parent.values())
        total_completed = sum(completed_per_parent.values())
        array_used = bool(some_parent_has_2 and total_completed >= 3)
    return m2, m3, array_used, array_parents, {"jobs": jobs, "within_budget": within}


# --- reading a built table and grading one dataset ----------------------------------


def read_table(output: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(output))


def grade_one_dataset(
    dataset: str, output: Path, reference: dict | None
) -> tuple[bool, bool, dict]:
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    info = DATASETS[dataset]
    details: dict = {"output": str(output)}
    details["success_marker"] = (output / "_SUCCESS").exists()
    sidecar = output.with_name(output.name + ".provenance.json")
    details["sidecar"] = sidecar.is_file()
    if not (details["success_marker"] and details["sidecar"]):
        return False, False, details
    prov = json.loads(sidecar.read_text())
    details["schema_id"] = prov.get("schema_id")
    details["provenance_sha256"] = hashlib.sha256(sidecar.read_bytes()).hexdigest()

    ht = read_table(output)
    schema = hail_schema_to_dict(ht)
    details["row_count"] = ht.count()

    # M4: the contract against the committed snapshot in THIS checkout, never from
    # reference.json and never from whatever checkout produced the output.
    committed_schema = json.loads(info["schema_snapshot"].read_text())
    contract = schema_contract(schema, committed_schema)
    details["schema_contract_ok"] = contract["ok"]
    details["missing_fields"] = contract["missing_fields"]
    details["type_mismatches"] = contract["type_mismatches"]
    details["extra_fields"] = contract["extra_fields"]
    m4 = details["schema_id"] == info["schema_id"] and bool(contract["ok"])

    if reference is None:
        return m4, False, details

    ref_ds = reference["datasets"][dataset]
    details["schema_equals_reference"] = schema == ref_ds["schema"]
    details["row_count_matches"] = details["row_count"] == ref_ds["row_count"]
    digest = rows_digest(ht, ref_ds["sample_keys"])
    details["rows_sha256"] = digest
    details["rows_digest_matches"] = digest == ref_ds["rows_sha256"]
    content_ok = bool(
        details["schema_equals_reference"]
        and details["row_count_matches"]
        and details["rows_digest_matches"]
    )
    return m4, content_ok, details


def write_reference(args) -> None:
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    out_dir = Path(args.out_dir)
    datasets = {}
    for ds, info in DATASETS.items():
        output = out_dir / info["output"]
        ht = read_table(output)
        schema = hail_schema_to_dict(ht)
        committed_schema = json.loads(info["schema_snapshot"].read_text())
        contract = schema_contract(schema, committed_schema)
        if not contract["ok"]:
            print(
                f"ERROR: --write-reference: {ds}'s build at {output} does not satisfy "
                f"its own committed schema snapshot ({info['schema_snapshot']}); the "
                f"committed snapshot and this build disagree on the contract -- "
                f"investigate before baselining. Detail: {json.dumps(contract)}",
                file=sys.stderr,
            )
            raise SystemExit(1)
        keys = choose_sample_keys(ht)
        datasets[ds] = {
            "output": info["output"],
            "schema_id": info["schema_id"],
            "schema": schema,
            "row_count": ht.count(),
            "sample_keys": keys,
            "rows_sha256": rows_digest(ht, keys),
        }
    ref = {
        "datasets": datasets,
        **BUDGET,
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "input_sha256": args.input_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(ref, indent=2, sort_keys=True) + "\n")
    counts = {ds: d["row_count"] for ds, d in datasets.items()}
    print(f"wrote {args.reference}: row counts {counts}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="holds the three fixed outputs")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--verdict-dir", default=None, help="default: --out-dir")
    ap.add_argument("--run-notes")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--expect-commit", default="")
    ap.add_argument("--expect-sif-sha256", default="")
    ap.add_argument("--input-sha256", default="")
    ap.add_argument("--site-profile")
    ap.add_argument("--sacct-file")
    ap.add_argument("--since", default="now-2days")
    ap.add_argument("--no-slurm", action="store_true")
    ap.add_argument("--expect-fail", action="store_true")
    ap.add_argument("--write-reference", action="store_true")
    args = ap.parse_args()

    if args.write_reference:
        write_reference(args)
        return 0
    for required in ("run_id", "expect_commit", "expect_sif_sha256"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required when grading")

    out_dir = Path(args.out_dir)
    pairs = [(ds, out_dir / info["output"]) for ds, info in DATASETS.items()]
    reference = json.loads(Path(args.reference).read_text())
    notes = parse_notes(
        Path(args.run_notes) if args.run_notes else Path("/nonexistent")
    )
    m1 = bool(
        notes["present"]
        and notes["commit"]
        and args.expect_commit.startswith(notes["commit"])
        and (not args.expect_sif_sha256 or notes["container"] == args.expect_sif_sha256)
    )

    slurm_details: dict
    if args.no_slurm:
        m2 = m3 = None
        array_used = None
        array_parents: list[str] = []
        slurm_details = {"skipped": "no-slurm"}
    else:
        rows, reason = load_sacct(args.sacct_file, args.since)
        if rows is None:
            m2 = m3 = None
            array_used = None
            array_parents = []
            slurm_details = {"skipped": reason}
        else:
            profile = Path(args.site_profile) if args.site_profile else None
            m2, m3, array_used, array_parents, slurm_details = grade_slurm(
                args.run_id, notes, rows, profile
            )

    datasets_out = {}
    m4_all = True
    content_ok_all = True
    for dataset, output in pairs:
        m4, content_ok, details = grade_one_dataset(dataset, output, reference)
        datasets_out[dataset] = details
        m4_all = m4_all and m4
        content_ok_all = content_ok_all and content_ok
    m5 = bool(content_ok_all and array_used is not False)
    passed = bool(m4_all and m5)

    verdict = {
        "task": "dataset-build-at-scale",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4_all,
        "M5": m5,
        "pass": passed,
        "content_ok": content_ok_all,
        "array_used": array_used,
        "array_parents": array_parents,
        "notes": notes,
        "slurm": slurm_details,
        "datasets": datasets_out,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else out_dir
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    receipt = {
        "datasets": {
            ds: {
                "output": datasets_out[ds]["output"],
                "provenance_sha256": datasets_out[ds].get("provenance_sha256"),
            }
            for ds, _ in pairs
        },
        "pass": passed,
        "reference_sha256": hashlib.sha256(
            Path(args.reference).read_bytes()
        ).hexdigest(),
        "checked_at": verdict["checked_at"],
    }
    (verdict_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(
        json.dumps(
            {
                k: verdict[k]
                for k in (
                    "M1",
                    "M2",
                    "M3",
                    "M4",
                    "M5",
                    "pass",
                    "content_ok",
                    "array_used",
                )
            }
        )
    )
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
