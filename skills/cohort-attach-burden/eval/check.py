"""Deterministic check for the cohort-attach-burden eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json`` under ``--out-dir`` (or ``--verdict-dir``),
and exits 0 on pass. Runs inside the container on a compute node (it needs Hail for
``attached.ht``). Slurm accounting is read from a file captured on the login node
(``--sacct-file``), because the container has no Slurm client.

There is no plugin-committed schema snapshot for a synthetic burden table or a
dynamically-composed attach output (unlike dataset-build-at-scale/gene-index), so the
CONTRACT FIELD NAMES for ``attached.ht`` come from this checkout's own code (the gene
spine's ``SPINE_FIELDS``, the burden manifest's ``ARCHITECTURE_AXIS_COLUMNS``, plus
``minp`` and ``cohort_tested``) while their observed TYPES come from ``reference.json``,
checked as a floor (extra fields allowed, same spirit as the rest of this batch).

M4 (the output exists and passes its contract):
  * ``gene_burden.tsv`` has every expected column (``BURDEN_COLUMNS`` below) and at
    least one row.
  * ``cohort.yaml`` loads with THIS checkout's own ``load_cohort``; its declared
    ``table`` resolves (absolute as written, else relative to the manifest's own
    directory -- see ``resolve_table_path``; this is the CHECK's own resolution rule
    for finding the file, not hvantk's own: ``hvantk cohort validate``/``attach`` take
    ``table`` verbatim, relative to the process's working directory, never relative to
    the manifest file -- see ``references/cohort-manifest.md``) to that same TSV;
    ``key`` is ``gene_id``; the prior is ``minp``/``lower_is_better``; some axis holds
    all three architecture columns.
  * ``attached.ht`` has ``_SUCCESS``, key ``["gene_id"]``, every contract field present
    with the reference's own stored type (a floor); its report JSON exists with
    ``n_genes`` and ``n_tested``.

M5: ``burden_ok`` (the TSV's row count and a digest over every row's contract columns
equal the reference's own), ``attach_ok`` (``n_genes``/``n_tested`` equal the
reference's; a digest over the TESTED rows' ``gene_id`` + declared columns equal the
reference's), ``reported_ok`` (the run notes' own ``attached: a/b`` line equals the
report's ``n_tested``/``n_genes``, in that order). ``M5 = burden_ok and attach_ok and
reported_ok``; ``pass = M4 and M5``.

    python check.py --out-dir DIR --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--verdict-dir DIR] [--expect-fail]
    python check.py --out-dir DIR --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

DIR holds ``gene_burden.tsv``, ``cohort.yaml``, ``attached.ht`` (and its default report
``attached.ht.report.json``). The accounting file is the output of
``sacct -X -n -P -S <since> -o JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition``.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

BUDGET = {"max_cpus": 8, "max_mem_gb": 32, "max_minutes": 30}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)

BURDEN_OUTPUT = "gene_burden.tsv"
COHORT_OUTPUT = "cohort.yaml"
ATTACHED_OUTPUT = "attached.ht"

# Every column hvantk's burden op writes for a --mtc run (hvantk/algorithms/burden/
# fet.py: pick_min_p + finalize_reductions + apply_mtc). Fixed by the code, never by
# the reference -- a column renamed upstream must fail M4 for every run, not drift the
# reference silently.
BURDEN_COLUMNS = (
    "gene",
    "route",
    "minp",
    "odds_ratio",
    "n_case_var",
    "conc",
    "driver_af",
    "mean_score_case",
    "frac_case_private",
    "p_adj",
)
BURDEN_NUMERIC_COLUMNS = {
    "minp",
    "odds_ratio",
    "n_case_var",
    "conc",
    "driver_af",
    "mean_score_case",
    "frac_case_private",
    "p_adj",
}

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_ATTACHED = re.compile(r"^attached:\s*(\d+)\s*/\s*(\d+)\s*$", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


# --- canonical digest (COMMON.md's shared convention) -------------------------------


def norm(value):
    """JSON-stable form for a digest: floats become a fixed-precision string (stable
    across platforms, finite for nan/inf), lists/tuples map elementwise, sets become
    sorted lists, dicts/Hail Structs become dicts, Locus/Interval render via str()."""
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, float):
        return format(value, ".10g")
    if isinstance(value, (int, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [norm(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(norm(v) for v in value)
    if isinstance(value, dict):
        return {k: norm(v) for k, v in value.items()}
    try:
        import hail as hl

        if isinstance(value, hl.Struct):
            return {k: norm(v) for k, v in dict(value).items()}
        if isinstance(value, (hl.Locus, hl.Interval)):
            return str(value)
    except ImportError:
        pass
    return str(value)


def rows_digest(rows) -> str:
    strings = sorted(
        json.dumps(norm(row), sort_keys=True, separators=(",", ":")) for row in rows
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


# --- run notes and Slurm accounting (M1-M3), the shared pattern across this batch ---


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {
            "present": False,
            "commit": None,
            "container": None,
            "jobs": [],
            "attached": None,
        }
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    a = _ATTACHED.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
        "attached": [int(a.group(1)), int(a.group(2))] if a else None,
    }


def mem_to_gb(s: str, cpus: int):
    m = _MEM.match(s.strip())
    if not m:
        return None
    unit = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1.0, "T": 1024.0}[m.group(2)]
    n = float(m.group(1)) * unit
    return n * cpus if m.group(3) == "c" else n


def limit_to_minutes(s: str):
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


def allowed_partitions(profile):
    if profile is None or not profile.is_file():
        return set()
    out = set()
    for line in profile.read_text().splitlines():
        for key in _PROFILE_PARTITIONS:
            m = re.match(rf"^\s*(?:export\s+)?{key}=[\"']?([^\"'#\s]+)", line)
            if m:
                out.add(m.group(1))
    return out


def parse_sacct(text: str):
    rows = []
    for line in text.splitlines():
        f = line.split("|")
        if len(f) < 7:
            continue
        rows.append(
            dict(zip(["id", "name", "state", "mem", "cpus", "limit", "partition"], f))
        )
    return rows


def load_sacct(sacct_file, since: str):
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


def grade_slurm(run_id, notes, rows, profile):
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


# --- burden TSV (no Hail) ------------------------------------------------------------


def read_burden(path: Path) -> list:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def normalized_burden_row(row: dict) -> dict:
    out = {}
    for col in BURDEN_COLUMNS:
        v = row.get(col, "")
        if col in BURDEN_NUMERIC_COLUMNS:
            out[col] = norm(float(v)) if v not in ("", None) else None
        else:
            out[col] = v
    return out


def grade_burden(out_dir: Path, reference) -> tuple:
    path = out_dir / BURDEN_OUTPUT
    details: dict = {"path": str(path), "exists": path.is_file()}
    if not details["exists"]:
        return False, details
    rows = read_burden(path)
    details["columns"] = list(rows[0].keys()) if rows else []
    details["n_rows"] = len(rows)
    details["columns_ok"] = all(c in details["columns"] for c in BURDEN_COLUMNS)
    m4 = bool(details["columns_ok"] and details["n_rows"] > 0)
    if not m4 or reference is None:
        return m4, details
    normalized = [normalized_burden_row(r) for r in rows]
    digest = rows_digest(normalized)
    details["rows_sha256"] = digest
    ref_burden = reference.get("burden", {})
    details["row_count_matches"] = details["n_rows"] == ref_burden.get("row_count")
    details["rows_digest_matches"] = digest == ref_burden.get("rows_sha256")
    details["burden_ok"] = bool(
        details["row_count_matches"] and details["rows_digest_matches"]
    )
    return m4, details


# --- cohort manifest (no Hail) -------------------------------------------------------


def resolve_table_path(table: str, manifest_path: Path) -> Path:
    """Where THIS CHECK looks for the manifest's declared table: absolute as written,
    else relative to the manifest's own directory. This is the check's own forgiving
    rule for locating the file to verify -- it is NOT what ``hvantk cohort
    validate``/``attach`` do: those take ``manifest.table`` verbatim and open it
    relative to the process's current working directory (see
    ``hvantk/algorithms/burden/manifest.py``'s ``render_cohort_manifest`` docstring and
    ``hvantk/tools/cohort/cohort_cli.py``'s ``_read_header``), never relative to the
    manifest file. A task that writes ``--output`` as an absolute path (as this one is
    told to) never hits the difference."""
    p = Path(table)
    return p if p.is_absolute() else (manifest_path.parent / p)


def grade_cohort(out_dir: Path) -> tuple:
    import jsonschema

    from hvantk.algorithms.burden.manifest import ARCHITECTURE_AXIS_COLUMNS
    from hvantk.algorithms.cohort.spec import load_cohort

    path = out_dir / COHORT_OUTPUT
    details: dict = {"path": str(path), "exists": path.is_file()}
    if not details["exists"]:
        return False, details
    try:
        manifest = load_cohort(path)
    except (jsonschema.ValidationError, ValueError) as exc:
        details["error"] = f"{type(exc).__name__}: {exc}"
        return False, details

    resolved_table = resolve_table_path(manifest.table, path)
    details["key"] = manifest.key
    details["key_ok"] = manifest.key == "gene_id"
    details["table_ok"] = (
        resolved_table.is_file()
        and resolved_table.resolve() == (out_dir / BURDEN_OUTPUT).resolve()
    )
    details["prior_column"] = manifest.prior.column
    details["prior_direction"] = manifest.prior.direction
    details["prior_ok"] = (
        manifest.prior.column == "minp"
        and manifest.prior.direction == "lower_is_better"
    )
    axis_cols = {c for entry in manifest.cohort_axes for c in entry.columns}
    details["axis_ok"] = set(ARCHITECTURE_AXIS_COLUMNS) <= axis_cols
    m4 = bool(
        details["key_ok"]
        and details["table_ok"]
        and details["prior_ok"]
        and details["axis_ok"]
    )
    return m4, details


# --- attached table (Hail) -----------------------------------------------------------


def read_table(path: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(path))


def contract_field_names() -> list:
    from hvantk.algorithms.annotation.spine import SPINE_FIELDS
    from hvantk.algorithms.burden.manifest import ARCHITECTURE_AXIS_COLUMNS

    return sorted(
        set(SPINE_FIELDS)
        | {"minp"}
        | set(ARCHITECTURE_AXIS_COLUMNS)
        | {"cohort_tested"}
    )


def declared_columns() -> list:
    from hvantk.algorithms.burden.manifest import ARCHITECTURE_AXIS_COLUMNS

    return ["minp"] + list(ARCHITECTURE_AXIS_COLUMNS)


def grade_attached(out_dir: Path, reference) -> tuple:
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    path = out_dir / ATTACHED_OUTPUT
    report_path = Path(f"{path}.report.json")
    details: dict = {
        "path": str(path),
        "success_marker": (path / "_SUCCESS").exists(),
        "report_exists": report_path.is_file(),
    }
    if not (details["success_marker"] and details["report_exists"]):
        return False, details

    report = json.loads(report_path.read_text())
    details["report"] = {k: report.get(k) for k in ("n_genes", "n_tested")}
    details["report_ok"] = "n_genes" in report and "n_tested" in report

    ht = read_table(path)
    schema = hail_schema_to_dict(ht)
    row_schema = schema.get("row", {})
    details["key"] = schema.get("key")
    details["key_ok"] = schema.get("key") == ["gene_id"]

    names = contract_field_names()
    details["missing_fields"] = [n for n in names if n not in row_schema]
    m4 = bool(
        details["key_ok"] and details["report_ok"] and not details["missing_fields"]
    )

    if reference is not None:
        ref_types = reference.get("attach", {}).get("contract_fields", {})
        type_mismatches = {
            n: {"expected": ref_types[n], "got": row_schema.get(n)}
            for n in ref_types
            if n in row_schema and row_schema.get(n) != ref_types[n]
        }
        details["type_mismatches"] = type_mismatches
        m4 = bool(m4 and not type_mismatches)

    details["n_genes"] = ht.count()
    if not m4 or reference is None:
        return m4, details

    ref_attach = reference.get("attach", {})
    details["n_genes_matches"] = details["n_genes"] == ref_attach.get("n_genes")
    tested = ht.filter(ht.cohort_tested)
    details["n_tested"] = tested.count()
    details["n_tested_matches"] = details["n_tested"] == ref_attach.get("n_tested")

    rows = [dict(r) for r in tested.select(*declared_columns()).collect()]
    digest = rows_digest(rows)
    details["rows_sha256"] = digest
    details["rows_digest_matches"] = digest == ref_attach.get("rows_sha256")
    details["attach_ok"] = bool(
        details["n_genes_matches"]
        and details["n_tested_matches"]
        and details["rows_digest_matches"]
    )
    return m4, details


# --- writing the reference -----------------------------------------------------------


def write_reference(args) -> None:
    out_dir = Path(args.out_dir)

    burden_m4, burden_details = grade_burden(out_dir, None)
    if not burden_m4:
        print(
            f"ERROR: --write-reference: {BURDEN_OUTPUT} at {out_dir} does not satisfy "
            f"its own contract (expected columns and >=1 row): {burden_details}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    rows = read_burden(out_dir / BURDEN_OUTPUT)
    normalized = [normalized_burden_row(r) for r in rows]
    burden_ref = {"row_count": len(rows), "rows_sha256": rows_digest(normalized)}

    cohort_m4, cohort_details = grade_cohort(out_dir)
    if not cohort_m4:
        print(
            f"ERROR: --write-reference: {COHORT_OUTPUT} at {out_dir} does not satisfy "
            f"the manifest contract: {cohort_details}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    path = out_dir / ATTACHED_OUTPUT
    report_path = Path(f"{path}.report.json")
    if not (path / "_SUCCESS").exists() or not report_path.is_file():
        print(
            f"ERROR: --write-reference: {ATTACHED_OUTPUT} at {out_dir} is missing "
            "_SUCCESS or its report json",
            file=sys.stderr,
        )
        raise SystemExit(1)
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    ht = read_table(path)
    schema = hail_schema_to_dict(ht)
    if schema.get("key") != ["gene_id"]:
        print(
            f"ERROR: --write-reference: {ATTACHED_OUTPUT}'s key is "
            f"{schema.get('key')}, expected ['gene_id']",
            file=sys.stderr,
        )
        raise SystemExit(1)
    names = contract_field_names()
    row_schema = schema.get("row", {})
    missing = [n for n in names if n not in row_schema]
    if missing:
        print(
            f"ERROR: --write-reference: {ATTACHED_OUTPUT} is missing contract "
            f"field(s): {missing}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    contract_types = {n: row_schema[n] for n in names}

    tested = ht.filter(ht.cohort_tested)
    rows = [dict(r) for r in tested.select(*declared_columns()).collect()]
    attach_ref = {
        "n_genes": ht.count(),
        "n_tested": tested.count(),
        "contract_fields": contract_types,
        "rows_sha256": rows_digest(rows),
    }

    ref = {
        "burden": burden_ref,
        "attach": attach_ref,
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
        f"wrote {args.reference}: burden rows={burden_ref['row_count']} "
        f"attach n_tested/n_genes={attach_ref['n_tested']}/{attach_ref['n_genes']}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--verdict-dir", default=None)
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
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else out_dir
    verdict_dir.mkdir(parents=True, exist_ok=True)

    reference = None
    ref_path = Path(args.reference)
    if ref_path.is_file():
        try:
            loaded = json.loads(ref_path.read_text())
        except ValueError:
            loaded = None
        # The committed placeholder ({"produced_by": null}) has no "burden"/"attach"
        # keys; treat it the same as "no reference yet" rather than crashing below.
        if loaded is not None and "burden" in loaded and "attach" in loaded:
            reference = loaded

    notes = parse_notes(
        Path(args.run_notes) if args.run_notes else Path("/nonexistent")
    )
    m1 = bool(
        notes["present"]
        and notes["commit"]
        and args.expect_commit.startswith(notes["commit"])
        and (not args.expect_sif_sha256 or notes["container"] == args.expect_sif_sha256)
    )

    if args.no_slurm:
        m2 = m3 = None
        slurm_details: dict = {"skipped": "no-slurm"}
    else:
        rows, reason = load_sacct(args.sacct_file, args.since)
        if rows is None:
            m2 = m3 = None
            slurm_details = {"skipped": reason}
        else:
            profile = Path(args.site_profile) if args.site_profile else None
            m2, m3, slurm_details = grade_slurm(args.run_id, notes, rows, profile)

    burden_details: dict = {}
    cohort_details: dict = {}
    attached_details: dict = {}
    error = None
    try:
        burden_m4, burden_details = grade_burden(out_dir, reference)
        cohort_m4, cohort_details = grade_cohort(out_dir)
        attached_m4, attached_details = grade_attached(out_dir, reference)
        m4 = bool(burden_m4 and cohort_m4 and attached_m4)
    except Exception as exc:  # noqa: BLE001 -- grading must not crash before the
        # verdict is written; record the failure and fail M4.
        m4 = False
        error = f"{type(exc).__name__}: {exc}"

    burden_ok = bool(reference is not None and burden_details.get("burden_ok"))
    attach_ok = bool(reference is not None and attached_details.get("attach_ok"))
    report = attached_details.get("report") or {}
    reported_ok = bool(
        reference is not None
        and notes.get("attached") is not None
        and report.get("n_tested") is not None
        and report.get("n_genes") is not None
        and notes["attached"] == [report["n_tested"], report["n_genes"]]
    )
    m5 = bool(burden_ok and attach_ok and reported_ok)
    passed = bool(m4 and m5)

    verdict = {
        "task": "cohort-attach-burden",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "burden_ok": burden_ok,
        "attach_ok": attach_ok,
        "reported_ok": reported_ok,
        "error": error,
        "notes": notes,
        "slurm": slurm_details,
        "burden": burden_details,
        "cohort": cohort_details,
        "attached": attached_details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )

    def _sha(p: Path):
        return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None

    receipt = {
        "burden_sha256": _sha(out_dir / BURDEN_OUTPUT),
        "cohort_sha256": _sha(out_dir / COHORT_OUTPUT),
        "pass": passed,
        "reference_sha256": hashlib.sha256(ref_path.read_bytes()).hexdigest(),
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
                    "burden_ok",
                    "attach_ok",
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
