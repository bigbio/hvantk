"""Deterministic check for the dataset-build-at-scale eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json`` next to the outputs, and exits 0 on pass.
Runs inside the container on a compute node (it needs Hail). Slurm accounting is read
from a file captured on the login node (``--sacct-file``), because the container has no
Slurm client; without a file it tries ``sacct`` and records M2 and M3 as null when that
is unavailable or fails.

Array-aware, unlike a single-job eval's check.py: the task asks for one array job
covering all three datasets, so M3's "every listed job exists in accounting" allows a
notes id ``A`` to be satisfied by an accounting row ``A`` *or* ``A_<task>`` (an array
submission writes one ``job: <A>`` line; every finished task is only ever accounted
under ``<A>_<task>``, never bare ``A``). ``array_used`` and ``array_parents`` record
whether the prefixed accounting rows actually show array-task-shaped ids, so the
negative control (the same three builds as three separate single-job submissions) can
be told apart from the positive one on more than content alone.

    python check.py --list TASKS.tsv --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--expect-fail]
    python check.py --list TASKS.tsv --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

TASKS.tsv is the same task list the skill's templates read:
``<provider:dataset>\\t<raw dir>\\t<output path>[\\t<extra options>]`` (only the first
and third fields are used here). The accounting file is the output of
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

#: One entry per dataset the eval builds. ``schema_snapshot`` and the sample-keys
#: import are the plugin's own committed test contract, never duplicated by value here.
DATASETS = {
    "hgnc:lookup": {
        "schema_id": "hgnc-lookup-v1",
        "schema_snapshot": REPO_ROOT / "hvantk/skills/hgnc/tests/snapshots/schema.json",
    },
    "clingen:gene-disease": {
        "schema_id": "clingen-gene-disease-v1",
        "schema_snapshot": REPO_ROOT
        / "hvantk/skills/clingen/tests/snapshots/schema.json",
    },
    "gencc:submissions": {
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


def sample_keys_for(dataset: str) -> list[dict]:
    """The plugin's own committed ``SAMPLE_KEYS``, imported rather than copied."""
    if dataset == "hgnc:lookup":
        from hvantk.skills.hgnc.tests.test_builder import SAMPLE_KEYS
    elif dataset == "clingen:gene-disease":
        from hvantk.skills.clingen.tests.test_builder import SAMPLE_KEYS
    elif dataset == "gencc:submissions":
        from hvantk.skills.gencc.tests.test_builder import SAMPLE_KEYS
    else:
        raise KeyError(dataset)
    return SAMPLE_KEYS


def sample_rows_sha256(rows: list[dict]) -> str:
    blob = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def load_list(path: str) -> list[tuple[str, Path]]:
    pairs = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError(f"{path}: line is missing fields: {line!r}")
        pairs.append((fields[0], Path(fields[2])))
    return pairs


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
    """True if *job_id* (from run-notes) is ``A`` or starts an array task ``A_<i>``.

    An array submission's run-notes line is the array's own id (``job: <A>``), but
    every finished task is accounted only as ``<A>_<i>``, never bare ``<A>`` -- the
    plain membership test a single-job eval uses would never match an array run.
    """
    return job_id in all_ids or any(i.startswith(job_id + "_") for i in all_ids)


def grade_slurm(
    run_id: str, notes: dict, rows: list[dict], profile: Path | None
) -> tuple[bool, bool, bool, list[str], dict]:
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
    completed_array_rows = [j for j in array_rows if j["state"] == "COMPLETED"]
    array_used = len(completed_array_rows) >= 3
    array_parents = sorted({j["id"].split("_", 1)[0] for j in array_rows})
    return (
        m2,
        m3,
        array_used,
        array_parents,
        {"jobs": jobs, "within_budget": within},
    )


def read_table(output: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(output))


def grade_one_dataset(dataset: str, output: Path, reference: dict | None):
    """(m4, content_ok, details) for one dataset. ``content_ok`` is False until a
    reference is available (write-reference mode never calls this)."""
    from hvantk.tests._snapshot_utils import collect_sample_rows, hail_schema_to_dict

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
    ref_ds = reference["datasets"][dataset] if reference else None
    expected_schema = ref_ds["schema"] if ref_ds else None
    details["schema_matches"] = (
        expected_schema is not None and schema == expected_schema
    )
    m4 = details["schema_id"] == DATASETS[dataset]["schema_id"] and bool(
        details["schema_matches"]
    )
    if ref_ds is None:
        return m4, False, details

    details["row_count_matches"] = details["row_count"] == ref_ds["row_count"]
    try:
        rows = collect_sample_rows(ht, keys=sample_keys_for(dataset))
        details["sample_keys_present"] = True
        details["sample_rows_sha256"] = sample_rows_sha256(rows)
    except KeyError:
        details["sample_keys_present"] = False
        details["sample_rows_sha256"] = None
    details["sample_rows_match"] = (
        details["sample_rows_sha256"] == ref_ds["sample_rows_sha256"]
    )
    content_ok = bool(
        details["row_count_matches"]
        and details["sample_keys_present"]
        and details["sample_rows_match"]
    )
    return m4, content_ok, details


def write_reference(args, pairs: list[tuple[str, Path]]) -> None:
    from hvantk.tests._snapshot_utils import collect_sample_rows

    datasets = {}
    for dataset, output in pairs:
        info = DATASETS[dataset]
        schema = json.loads(info["schema_snapshot"].read_text())
        ht = read_table(output)
        keys = sample_keys_for(dataset)
        rows = collect_sample_rows(ht, keys=keys)
        datasets[dataset] = {
            "output": output.name,
            "schema_id": info["schema_id"],
            "schema": schema,
            "row_count": ht.count(),
            "sample_keys": keys,
            "sample_rows_sha256": sample_rows_sha256(rows),
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
    ap.add_argument("--list", dest="list_path", required=True)
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

    pairs = load_list(args.list_path)
    unknown = [ds for ds, _ in pairs if ds not in DATASETS]
    if unknown:
        ap.error(f"unknown dataset(s) in {args.list_path}: {unknown}")

    if args.write_reference:
        write_reference(args, pairs)
        return 0
    for required in ("run_id", "expect_commit", "expect_sif_sha256"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required when grading")

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
    first_output = pairs[0][1] if pairs else None
    out_dir = (
        first_output.parent
        if first_output is not None and first_output.parent.exists()
        else Path(".")
    )
    (out_dir / "verdict.json").write_text(json.dumps(verdict, indent=2, default=str))
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
    (out_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
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
