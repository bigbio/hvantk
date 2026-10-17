"""Deterministic check for the callset-subset-export eval task.

M1 from run-notes.md; M2 and M3 from an accounting file; M4 = the delivery is complete
(VCF present, SHA256SUMS validates, README present); M5 = the VCF's header samples
equal the reference's and its record count equals the reference's. No Hail: the
block-gzipped VCF is read as gzip, so this runs through the image on a login node.

    python check.py --output OUT/subset --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --no-slurm] [--site-profile FILE] [--prefix subset] [--expect-fail]
    python check.py --output OUT/subset --reference reference.json --write-reference \
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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sums_validate(out: Path) -> bool:
    sums = out / "SHA256SUMS"
    if not sums.is_file():
        return False
    ok = False
    for line in sums.read_text().splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            return False
        digest, name = parts
        p = out / name.strip().lstrip("*")
        if not p.is_file() or sha256_file(p) != digest:
            return False
        ok = True
    return ok


def inspect_vcf(path: Path) -> dict:
    samples: list[str] = []
    records = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                records += 1
    return {"samples": sorted(samples), "records": records}


def grade_output(out: Path, prefix: str, reference: dict | None):
    vcf = out / f"{prefix}.vcf.bgz"
    details: dict = {
        "vcf": vcf.is_file(),
        "sums_validate": sums_validate(out),
        "readme": (out / "README.md").is_file()
        and (out / "README.md").stat().st_size > 0,
    }
    m4 = all(details.values())
    if not details["vcf"]:
        return False, False, details
    details.update(inspect_vcf(vcf))
    if reference is None:
        return m4, False, details
    m5 = (
        details["samples"] == reference["samples"]
        and details["records"] == reference["records"]
    )
    return m4, m5, details


def write_reference(args) -> None:
    info = inspect_vcf(Path(args.output) / f"{args.prefix}.vcf.bgz")
    ref = {
        "samples": info["samples"],
        "records": info["records"],
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
        f"wrote {args.reference}: samples={info['samples']} records={info['records']}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--prefix", default="subset")
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
        "task": "callset-subset-export",
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
    vcf = out / f"{args.prefix}.vcf.bgz"
    receipt = {
        "output": str(out),
        "pass": passed,
        "vcf_sha256": sha256_file(vcf) if vcf.is_file() else None,
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
