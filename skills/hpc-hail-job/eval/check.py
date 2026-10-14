"""Deterministic check for the hpc-hail-job eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json`` next to the output, and exits 0 on pass.
Runs inside the container on a compute node (it needs Hail). Slurm accounting is read
from a file captured on the login node (``--sacct-file``), because the container has no
Slurm client; without a file it tries ``sacct`` and records M2 and M3 as null when that
is unavailable or fails.

    python check.py --output OUT.ht --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA [--expect-sif-sha256 HEX] \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--expect-fail]
    python check.py --output OUT.ht --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

The accounting file is the output of
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
SCHEMA_SNAPSHOT = REPO_ROOT / "hvantk/skills/hgnc/tests/snapshots/schema.json"
SCHEMA_ID = "hgnc-lookup-v1"
BUDGET = {"max_cpus": 8, "max_mem_gb": 32, "max_minutes": 30}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


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


def grade_slurm(
    run_id: str, notes: dict, rows: list[dict], profile: Path | None
) -> tuple[bool, bool, dict]:
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
    m3 = any(j["state"].startswith(SUCCESS_STATES) for j in jobs) and all(
        n in all_ids for n in notes["jobs"]
    )
    return m2, m3, {"jobs": jobs, "within_budget": within}


def read_table(output: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(output))


def grade_output(output: Path, reference: dict | None) -> tuple[bool, bool, dict]:
    from hvantk.tests._snapshot_utils import collect_sample_rows, hail_schema_to_dict

    details: dict = {}
    sidecar = output.with_name(output.name + ".provenance.json")
    details["success_marker"] = (output / "_SUCCESS").exists()
    details["sidecar"] = sidecar.is_file()
    if not (details["success_marker"] and details["sidecar"]):
        return False, False, details
    prov = json.loads(sidecar.read_text())
    details["schema_id"] = prov.get("schema_id")
    ht = read_table(output)
    schema = hail_schema_to_dict(ht)
    expected = json.loads(SCHEMA_SNAPSHOT.read_text())
    details["schema_matches_snapshot"] = schema == expected
    m4 = prov.get("schema_id") == SCHEMA_ID and schema == expected
    details["row_count"] = ht.count()
    if reference is None:
        return m4, False, details
    try:
        collect_sample_rows(ht, keys=reference["sample_keys"])
        details["sample_keys_present"] = True
    except KeyError:
        details["sample_keys_present"] = False
    m5 = (
        details["row_count"] == reference["row_count"]
        and details["sample_keys_present"]
    )
    return m4, m5, details


def write_reference(args) -> None:
    from hvantk.skills.hgnc.tests.test_builder import SAMPLE_KEYS

    ht = read_table(Path(args.output))
    ref = {
        "schema_id": SCHEMA_ID,
        "row_count": ht.count(),
        "sample_keys": SAMPLE_KEYS,
        **BUDGET,
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "input_sha256": args.input_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(ref, indent=2) + "\n")
    print(f"wrote {args.reference}: row_count={ref['row_count']}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--reference", required=True)
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

    output = Path(args.output)
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
        slurm_details = {"skipped": "no-slurm"}
    else:
        rows, reason = load_sacct(args.sacct_file, args.since)
        if rows is None:
            m2 = m3 = None
            slurm_details = {"skipped": reason}
        else:
            profile = Path(args.site_profile) if args.site_profile else None
            m2, m3, slurm_details = grade_slurm(args.run_id, notes, rows, profile)
    m4, m5, out_details = grade_output(output, reference)
    passed = bool(m4 and m5)

    verdict = {
        "task": "hpc-hail-job",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "notes": notes,
        "slurm": slurm_details,
        "output": out_details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    out_dir = output.parent if output.parent.exists() else Path(".")
    (out_dir / "verdict.json").write_text(json.dumps(verdict, indent=2, default=str))
    sidecar = output.with_name(output.name + ".provenance.json")
    receipt = {
        "output": str(output),
        "pass": passed,
        "provenance_sha256": (
            hashlib.sha256(sidecar.read_bytes()).hexdigest()
            if sidecar.is_file()
            else None
        ),
        "reference_sha256": hashlib.sha256(
            Path(args.reference).read_bytes()
        ).hexdigest(),
        "checked_at": verdict["checked_at"],
    }
    (out_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(json.dumps({k: verdict[k] for k in ("M1", "M2", "M3", "M4", "M5", "pass")}))
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
