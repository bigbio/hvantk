"""Deterministic check for the ml-features eval task.

Grades one run in five stages (M1 to M5, see skills/workflow-conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json``, and exits 0 on pass. Runs inside the
container on a compute node (it needs Hail). Slurm accounting is read from a file
captured on the login node (``--sacct-file``), because the container has no Slurm
client; without a file it tries ``sacct`` and records M2 and M3 as null when that is
unavailable or fails.

The four outputs are read by the skill's own receipt template,
``templates/describe-training-set.py`` (loaded by path, so the check and the skill's
Verification step cannot disagree on what the outputs say).

M4 is the output contract plus proof the inputs were built: the four files
``clinvar_chr22.{ht,tsv.gz,dictionary.tsv,manifest.json}`` exist, the ``.ht`` is keyed
``[locus, alleles]``, the template finds the four mutually consistent (same rows, same
columns in the same order, a tier on every feature, the manifest's tier counts equal to
the dictionary's), no training feature is typed ``str``, and ``--out-dir`` holds, up to
3 directories deep, built tables whose ``.provenance.json`` sidecars name the schema ids
of ``clinvar:variants``, ``dbnsfp:variants``, ``insider:variants`` and
``uniprot-ptm:sites``. Without that last check an agent could pass by copying a finished
dataset from elsewhere.

M5 is ``content_ok`` (the TSV's SHA-256 -- the export is byte-identical across builds --
and the row count, label counts and tier counts all equal the reference) AND
``reported_ok`` (the run notes' ``rows:`` and ``positives:`` equal the manifest's,
whatever the reference says). The negative control is a real, self-consistent build
from a spec with ``min_review_stars: 0``: ``reported_ok`` true, ``content_ok`` false.

Grading never crashes the run: reading the outputs, checking the inputs and grading the
notes happen inside one try/except in ``main()``, so a corrupt output records ``error``
and fails M4/M5 alone; ``verdict.json`` and ``receipt.json`` are still written.

    python check.py --out-dir DIR --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--verdict-dir DIR] [--expect-fail]
    python check.py --out-dir DIR --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

The accounting file is the output of
``sacct -X -n -P -S <since> -o JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 45}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)
PREFIX_NAME = "clinvar_chr22"
EXPECTED_KEY = ["locus", "alleles"]
#: Schema ids of the variant-level input datasets the agent must build, from their
#: plugin manifests (hvantk/skills/{clinvar,dbnsfp,insider,uniprot_ptm}/plugin.yaml).
REQUIRED_SCHEMAS = (
    "clinvar-variants-v1",
    "dbnsfp-v1",
    "insider-variants-v1",
    "uniprot-ptm-sites-v1",
)
#: How many directories below --out-dir the input builds may sit.
INPUTS_MAX_DEPTH = 3
TEMPLATE = (
    Path(__file__).resolve().parent.parent / "templates" / "describe-training-set.py"
)

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_ROWS = re.compile(r"^rows:\s*(\d+)\s*$", re.M)
_POSITIVES = re.compile(r"^positives:\s*(\d+)\s*$", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


def load_template():
    """The skill's receipt template, imported by path (its name has a hyphen)."""
    spec = importlib.util.spec_from_file_location("describe_training_set", TEMPLATE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- were the inputs built? -----------------------------------------------------------


def find_input_sidecars(out_dir: Path, max_depth: int = INPUTS_MAX_DEPTH) -> list:
    """Every ``*.provenance.json`` sidecar under *out_dir*, at most *max_depth*
    directories below it (depth 0: directly inside ``out_dir``)."""
    out_dir = out_dir.resolve()
    found = []
    for p in out_dir.rglob("*.provenance.json"):
        rel = p.relative_to(out_dir)
        if len(rel.parts) - 1 <= max_depth:
            found.append(p)
    return found


def check_inputs_built(out_dir: Path) -> dict:
    """Did ``--out-dir`` receive a real build of each variant-level input, identified by
    its sidecar's own ``schema_id`` -- never by filename, which the task leaves to the
    agent ("under $HVANTK_EVAL_OUT")."""
    by_schema: dict = {}
    for sidecar in find_input_sidecars(out_dir):
        try:
            payload = json.loads(sidecar.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        schema_id = payload.get("schema_id")
        if schema_id in REQUIRED_SCHEMAS and schema_id not in by_schema:
            by_schema[schema_id] = {
                "sidecar": str(sidecar),
                "sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            }
    missing = [s for s in REQUIRED_SCHEMAS if s not in by_schema]
    return {"ok": not missing, "found": by_schema, "missing": missing}


# --- run notes and Slurm accounting (M1-M3) -----------------------------------------


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {
            "present": False,
            "commit": None,
            "container": None,
            "jobs": [],
            "rows": None,
            "positives": None,
        }
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    r = _ROWS.search(text)
    p = _POSITIVES.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
        "rows": int(r.group(1)) if r else None,
        "positives": int(p.group(1)) if p else None,
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


def allowed_partitions(profile: Path | None) -> set:
    if profile is None or not profile.is_file():
        return set()
    out = set()
    for line in profile.read_text().splitlines():
        for key in _PROFILE_PARTITIONS:
            m = re.match(rf"^\s*(?:export\s+)?{key}=[\"']?([^\"'#\s]+)", line)
            if m:
                out.add(m.group(1))
    return out


def parse_sacct(text: str) -> list:
    rows = []
    for line in text.splitlines():
        f = line.split("|")
        if len(f) < 7:
            continue
        rows.append(
            dict(zip(["id", "name", "state", "mem", "cpus", "limit", "partition"], f))
        )
    return rows


def load_sacct(sacct_file: str | None, since: str):
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


def job_known(job_id: str, all_ids: set) -> bool:
    """True if *job_id* (from run-notes) is ``A`` or starts an array task ``A_<i>``
    (an array submission writes one ``job: A`` line; a finished task is only ever
    accounted under ``A_<i>``, never bare ``A`` -- see dataset-build-at-scale's
    check.py, which this is copied from; an agent may build the inputs as one array)."""
    return job_id in all_ids or any(i.startswith(job_id + "_") for i in all_ids)


def grade_slurm(run_id: str, notes: dict, rows: list, profile: Path | None):
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
        and all(job_known(n, all_ids) for n in notes["jobs"])
    )
    return m2, m3, {"jobs": jobs, "within_budget": within}


# --- reading the training set and grading it against the reference ------------------


def describe_output(out_dir: Path) -> tuple[dict, list]:
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    tpl = load_template()
    desc = tpl.describe(str(out_dir / PREFIX_NAME))
    return desc, tpl.consistency_problems(desc)


def grade_training_set(out_dir: Path, reference: dict):
    details: dict = {"prefix": str(out_dir / PREFIX_NAME)}
    files = {
        s: (out_dir / f"{PREFIX_NAME}{s}").exists()
        for s in (".ht", ".tsv.gz", ".dictionary.tsv", ".manifest.json")
    }
    details["files"] = files
    details["success_marker"] = (out_dir / f"{PREFIX_NAME}.ht" / "_SUCCESS").exists()
    details["inputs_built"] = check_inputs_built(out_dir)
    if not (all(files.values()) and details["success_marker"]):
        return False, False, details

    desc, problems = describe_output(out_dir)
    details["describe"] = {
        k: v for k, v in desc.items() if k not in ("ht_columns", "tsv_header")
    }
    details["problems"] = problems
    m4 = bool(
        desc["key"] == EXPECTED_KEY
        and not problems
        and not desc["string_features"]
        and details["inputs_built"]["ok"]
    )
    for field in ("tsv_sha256", "rows", "label_counts", "tier_counts"):
        if field not in reference:
            raise ValueError(
                f"reference.json has no {field!r}; rerun the control with WRITE_REFERENCE=1"
            )
    content_ok = bool(
        desc["tsv_sha256"] == reference["tsv_sha256"]
        and desc["rows"] == reference["rows"]
        and desc["label_counts"] == reference["label_counts"]
        and desc["tier_counts"] == reference["tier_counts"]
    )
    details["content_ok"] = content_ok
    return m4, content_ok, details


def grade_reported(notes: dict, desc: dict | None) -> bool:
    """Self-consistency between the run notes' ``rows:``/``positives:`` lines and the
    dataset actually graded -- independent of ``content_ok``. Never raises."""
    if not desc or notes.get("rows") is None or notes.get("positives") is None:
        return False
    return bool(
        notes["rows"] == desc.get("rows")
        and notes["positives"] == desc.get("label_counts", {}).get("1")
    )


def write_reference(args) -> None:
    out_dir = Path(args.out_dir)
    desc, problems = describe_output(out_dir)
    inputs = check_inputs_built(out_dir)
    fatal = list(problems)
    if desc["key"] != EXPECTED_KEY:
        fatal.append(f"key {desc['key']}, expected {EXPECTED_KEY}")
    if desc["string_features"]:
        fatal.append(f"str-typed features: {desc['string_features']}")
    if not inputs["ok"]:
        fatal.append(f"input builds missing: {inputs['missing']}")
    if fatal:
        print("ERROR: --write-reference: " + "; ".join(fatal), file=sys.stderr)
        raise SystemExit(1)
    ref = {
        "tsv_sha256": desc["tsv_sha256"],
        "rows": desc["rows"],
        "label_counts": desc["label_counts"],
        "tier_counts": desc["tier_counts"],
        "excluded": desc["excluded"],
        "funnel": desc["funnel"],
        "columns": desc["ht_columns"],
        **BUDGET,
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "input_sha256": args.input_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(ref, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.reference}: rows={desc['rows']} labels={desc['label_counts']} "
        f"tiers={desc['tier_counts']}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out-dir", required=True, help="holds clinvar_chr22.* and the input builds"
    )
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

    # Loading reference.json, reading the outputs, checking the inputs and grading the
    # run notes all happen in one try/except: a missing or corrupt reference.json, or a
    # crash anywhere in here, must still leave verdict.json and receipt.json written.
    error: str | None = None
    try:
        reference = json.loads(Path(args.reference).read_text())
        m4, content_ok, ts_details = grade_training_set(out_dir, reference)
        reported_ok = grade_reported(notes, ts_details.get("describe"))
    except Exception as exc:  # noqa: BLE001 -- see the comment above
        error = f"{type(exc).__name__}: {exc}"
        m4, content_ok, reported_ok = False, False, False
        ts_details = {"prefix": str(out_dir / PREFIX_NAME), "error": error}

    m5 = bool(content_ok and reported_ok)
    passed = bool(m4 and m5)

    verdict = {
        "task": "ml-features",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "content_ok": content_ok,
        "reported_ok": reported_ok,
        "error": error,
        "notes": notes,
        "slurm": slurm_details,
        "training_set": ts_details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else out_dir
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    described = ts_details.get("describe") or {}
    receipt = {
        "prefix": ts_details.get("prefix"),
        "tsv_sha256": described.get("tsv_sha256"),
        "sidecar_sha256": described.get("sidecar_sha256"),
        "rows": described.get("rows"),
        "label_counts": described.get("label_counts"),
        "tier_counts": described.get("tier_counts"),
        "inputs_built": ts_details.get("inputs_built"),
        "pass": passed,
        "reference_sha256": (
            hashlib.sha256(Path(args.reference).read_bytes()).hexdigest()
            if Path(args.reference).is_file()
            else None
        ),
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
                    "reported_ok",
                )
            }
        )
    )
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
