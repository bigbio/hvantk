"""Deterministic check for the joint-genotyping eval task.

M1 from run-notes.md; M2 and M3 from an accounting file captured on the login node;
M4 = the callset's artifacts exist (VDS with both _SUCCESS markers, MatrixTable with
_SUCCESS, sample QC CSV, VCF); M5 = the sample set, the MatrixTable row count and the
VCF record count equal the reference's. Runs inside the image on a compute node (Hail).

    python check.py --output OUT/cohort --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --no-slurm] [--site-profile FILE] [--prefix cohort] [--expect-fail]
    python check.py --output OUT/cohort --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

BUDGET = {"max_cpus": 8, "max_mem_gb": 32, "max_minutes": 30}
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
    m = _MEM.match(s.strip())
    if not m:
        return None
    unit = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1.0, "T": 1024.0}[m.group(2)]
    n = float(m.group(1)) * unit
    return n * cpus if m.group(3) == "c" else n


def limit_to_minutes(s: str) -> float | None:
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
        if len(f) >= 7:
            rows.append(
                dict(
                    zip(["id", "name", "state", "mem", "cpus", "limit", "partition"], f)
                )
            )
    return rows


def grade_slurm(run_id: str, notes: dict, rows: list[dict], profile: Path | None):
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
        within.append(
            cpus <= BUDGET["max_cpus"]
            and gb is not None
            and gb <= BUDGET["max_mem_gb"] + 0.5
            and mins is not None
            and mins <= BUDGET["max_minutes"] + 0.01
            and (not allowed or j["partition"] in allowed)
        )
    m2 = bool(jobs) and all(within)
    all_ids = {j["id"] for j in rows}
    m3 = (
        bool(notes["jobs"])
        and any(j["state"].startswith(SUCCESS_STATES) for j in jobs)
        and all(n in all_ids for n in notes["jobs"])
    )
    return m2, m3, {"jobs": jobs, "within_budget": within}


def vcf_records(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


def describe(out: Path, prefix: str) -> dict:
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    mt = hl.read_matrix_table(str(out / f"{prefix}.mt"))
    return {
        "samples": sorted(mt.s.collect()),
        "mt_rows": mt.count_rows(),
        "vcf_records": vcf_records(out / f"{prefix}.vcf.bgz"),
    }


def grade_output(out: Path, prefix: str, reference: dict | None):
    details: dict = {
        "vds_success": all(
            (out / f"{prefix}.vds" / part / "_SUCCESS").exists()
            for part in ("reference_data", "variant_data")
        ),
        "mt_success": (out / f"{prefix}.mt" / "_SUCCESS").exists(),
        "sample_qc_csv": (out / "qc" / f"{prefix}_sample_qc.csv").is_file(),
        "vcf": (out / f"{prefix}.vcf.bgz").is_file(),
    }
    m4 = all(details.values())
    if not m4:
        return False, False, details
    details.update(describe(out, prefix))
    if reference is None:
        return m4, False, details
    m5 = (
        details["samples"] == reference["samples"]
        and details["mt_rows"] == reference["mt_rows"]
        and details["vcf_records"] == reference["vcf_records"]
    )
    return m4, m5, details


def write_reference(args) -> None:
    d = describe(Path(args.output), args.prefix)
    ref = {
        "samples": d["samples"],
        "mt_rows": d["mt_rows"],
        "vcf_records": d["vcf_records"],
        **BUDGET,
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "input_sha256": args.input_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(ref, indent=2) + "\n")
    print(
        f"wrote {args.reference}: samples={len(d['samples'])} mt_rows={d['mt_rows']} vcf_records={d['vcf_records']}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--prefix", default="cohort")
    ap.add_argument("--run-notes")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--expect-commit", default="")
    ap.add_argument("--expect-sif-sha256", default="")
    ap.add_argument("--input-sha256", default="")
    ap.add_argument("--site-profile")
    ap.add_argument("--sacct-file")
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

    out = Path(args.output)
    reference = json.loads(Path(args.reference).read_text())
    notes = parse_notes(
        Path(args.run_notes) if args.run_notes else Path("/nonexistent")
    )
    m1 = bool(
        notes["present"]
        and notes["commit"]
        and args.expect_commit.startswith(notes["commit"])
        and notes["container"] == args.expect_sif_sha256
    )
    if args.no_slurm:
        m2 = m3 = None
        slurm_details: dict = {"skipped": "no-slurm"}
    elif args.sacct_file and Path(args.sacct_file).is_file():
        profile = Path(args.site_profile) if args.site_profile else None
        m2, m3, slurm_details = grade_slurm(
            args.run_id, notes, parse_sacct(Path(args.sacct_file).read_text()), profile
        )
    else:
        m2 = m3 = None
        slurm_details = {"skipped": "no accounting file"}
    m4, m5, out_details = grade_output(out, args.prefix, reference)
    passed = bool(m4 and m5)

    verdict = {
        "task": "joint-genotyping",
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
    out_dir = out.parent if out.parent.exists() else Path(".")
    (out_dir / "verdict.json").write_text(json.dumps(verdict, indent=2, default=str))
    state = json.dumps(
        {k: out_details.get(k) for k in ("samples", "mt_rows", "vcf_records")},
        sort_keys=True,
    )
    receipt = {
        "output": str(out),
        "pass": passed,
        "counts_sha256": hashlib.sha256(state.encode()).hexdigest(),
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
