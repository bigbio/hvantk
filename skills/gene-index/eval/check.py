"""Deterministic check for the gene-index eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json``, and exits 0 on pass. Runs inside the
container on a compute node (it needs Hail). Slurm accounting is read from a file
captured on the login node (``--sacct-file``), because the container has no Slurm
client; without a file it tries ``sacct`` and records M2 and M3 as null when that is
unavailable or fails.

M4 is the schema CONTRACT plus proof the two inputs were actually built, not just the
gene index: the gene index's key must be ``["gene_id"]``, every name in THIS
checkout's own ``SPINE_FIELDS`` (``hvantk.algorithms.annotation.spine``) must be
present with the type ``reference.json``'s ``schema`` recorded for it (extra fields
beyond ``SPINE_FIELDS`` are allowed, and the key field itself is never counted as
"extra" even though Hail's own row schema always includes it), AND ``--out-dir``
must contain, up to 3 directories deep, a built table whose ``.provenance.json``
sidecar names ``schema_id: ensembl-gene-structure-v1`` and one whose sidecar names
``schema_id: hgnc-lookup-v1`` (the two schema ids the ``ensembl-gene`` and ``hgnc``
plugin manifests declare). Without that last check, an agent could pass the task by
copying a finished gene index from another skill's staged data, with no build of its
own. There is no fixture-derived schema snapshot for the gene index itself (unlike a
plugin's own committed contract): the reference IS the only schema-type source there
is, because `annotate spine` stamps no sidecar of its own.

M5 is ``content_ok`` (this run reproduces the reference build: row count, HGNC mapping
rate within a tolerance, and a digest over every row restricted to the contract fields
-- key plus ``SPINE_FIELDS`` -- all equal to the reference) AND ``reported_ok`` (the
run notes' own ``genes:`` and ``mapping_rate:`` lines are self-consistent with the
table actually graded, independent of whether that table matches the reference). Both
tolerances come from ``reference.json``'s own ``tolerances`` object (written by
``--write-reference``), falling back to this module's defaults for an older reference
that predates it. The negative control is a build that is perfectly self-consistent
(``reported_ok`` true) but does not reproduce the reference (``content_ok`` false) --
see eval/control.sh.

Grading never crashes the run: reading the table, checking the inputs, and grading the
run notes are all called inside one try/except in ``main()``, so a corrupt or
wrong-type output, or a malformed reported value, records ``error`` in the details and
fails M4/M5 alone -- ``verdict.json`` and ``receipt.json`` are still written.

    python check.py --out-dir DIR --run-notes NOTES.md --reference reference.json \
        --run-id PREFIX --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--verdict-dir DIR] [--expect-fail]
    python check.py --out-dir DIR --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

DIR holds the fixed-name output ``gene-index.ht``, plus (somewhere within 3 directories
of it) the two input builds. The accounting file is the output of
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

BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 30}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)
GENE_INDEX_NAME = "gene-index.ht"
EXPECTED_KEY = ["gene_id"]
MAPPING_RATE_TOLERANCE = 1e-9
REPORTED_RATE_TOLERANCE = 6e-5
#: The two input datasets' schema ids, from their plugin manifests
#: (hvantk/skills/ensembl_gene/plugin.yaml, hvantk/skills/hgnc/plugin.yaml).
ENSEMBL_SCHEMA_ID = "ensembl-gene-structure-v1"
HGNC_SCHEMA_ID = "hgnc-lookup-v1"
#: How many directories below --out-dir the two input builds may sit.
INPUTS_MAX_DEPTH = 3

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_GENES = re.compile(r"^genes:\s*(\d+)\s*$", re.M)
_MAPPING_RATE = re.compile(r"^mapping_rate:\s*(\d+(?:\.\d+)?%?)\s*$", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


# --- the digest: order-insensitive, contract fields only (key + SPINE_FIELDS) -------


def norm(value):
    """Floats to a stable string (so ``nan``/``inf`` digest the same way every run);
    everything else recurses or passes through. See skills/_conventions' shared digest
    convention -- each eval/check.py keeps its own copy, self-contained."""
    if isinstance(value, float):
        return format(value, ".10g")
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [norm(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(norm(v) for v in value)
    if isinstance(value, dict):
        return {k: norm(v) for k, v in value.items()}
    return value


def rows_digest(ht, fields: list) -> str:
    """SHA-256 over every row, projected to *fields* (the CONTRACT: key fields stay on
    the row regardless of ``select()``, so only the non-key fields need naming).
    Order-insensitive. The gene index is small enough (tens of thousands of rows) that
    collecting every row, rather than a fixed-size sample, is cheap."""
    from hvantk.tests._snapshot_utils import _to_jsonable

    key_field_names = list(ht.key)
    non_key_fields = [f for f in fields if f not in key_field_names]
    collected = ht.select(*non_key_fields).collect()
    strings = sorted(
        json.dumps(norm(_to_jsonable(dict(row))), sort_keys=True, separators=(",", ":"))
        for row in collected
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


# --- schema contract (a floor: SPINE_FIELDS names present with reference's stored
# --- type; extra fields allowed; the key field is never "extra" even though Hail's
# --- own row schema always carries it alongside the row-value fields) ----------------


def schema_contract(
    built_row: dict, ref_row: dict, spine_fields: list, key_fields: list
) -> dict:
    missing = [f for f in spine_fields if f not in built_row]
    mismatches = {}
    for f in spine_fields:
        if f in built_row and f in ref_row and built_row[f] != ref_row[f]:
            mismatches[f] = {"expected": ref_row[f], "got": built_row[f]}
    extra = sorted(set(built_row) - set(spine_fields) - set(key_fields))
    return {
        "ok": not missing and not mismatches,
        "missing_fields": missing,
        "type_mismatches": mismatches,
        "extra_fields": extra,
    }


# --- proof the two inputs were actually built (M4), not merely assumed -------------


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
    """Did ``--out-dir`` actually receive a real ``ensembl-gene:structure`` build and a
    real ``hgnc:lookup`` build, each identified by its sidecar's own ``schema_id`` --
    never by filename, which the task leaves to the agent ("anywhere under
    $HVANTK_EVAL_OUT")."""
    by_schema: dict = {}
    for sidecar in find_input_sidecars(out_dir):
        try:
            payload = json.loads(sidecar.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        schema_id = payload.get("schema_id")
        if (
            schema_id in (ENSEMBL_SCHEMA_ID, HGNC_SCHEMA_ID)
            and schema_id not in by_schema
        ):
            suffix = ".provenance.json"
            table = sidecar.with_name(sidecar.name[: -len(suffix)])
            by_schema[schema_id] = {
                "table": str(table),
                "sidecar": str(sidecar),
                "sidecar_sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            }
    ensembl = by_schema.get(ENSEMBL_SCHEMA_ID)
    hgnc = by_schema.get(HGNC_SCHEMA_ID)
    return {
        "ok": bool(ensembl and hgnc),
        "ensembl_gene_structure": ensembl,
        "hgnc_lookup": hgnc,
    }


# --- run notes and Slurm accounting (M1-M3) -----------------------------------------


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {
            "present": False,
            "commit": None,
            "container": None,
            "jobs": [],
            "genes": None,
            "mapping_rate_raw": None,
        }
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    g = _GENES.search(text)
    r = _MAPPING_RATE.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
        "genes": int(g.group(1)) if g else None,
        "mapping_rate_raw": r.group(1) if r else None,
    }


def parse_reported_rate(raw: str | None) -> float | None:
    """``mapping_rate:``'s value as a fraction, or ``None`` if it cannot be parsed
    (never raises -- a malformed value must fail ``reported_ok``, not crash the run).

    A trailing ``%`` always means percentage. With no ``%``, a value greater than 1
    is still read as a percentage written without the sign (e.g. an agent writing
    ``96.66`` instead of ``96.66%`` or ``0.9666``); a value at or below 1 is already a
    fraction.
    """
    if raw is None:
        return None
    raw = raw.strip()
    try:
        if raw.endswith("%"):
            return float(raw[:-1]) / 100.0
        value = float(raw)
    except ValueError:
        return None
    return value / 100.0 if value > 1 else value


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
    check.py, which this is copied from; the gene-index procedure recommends an array
    for the two input builds)."""
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


# --- reading the gene index and grading it against the reference -------------------


def read_table(output: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(output))


def grade_gene_index(output: Path, reference: dict | None):
    from hvantk.algorithms.annotation.spine import SPINE_FIELDS, spine_mapping_rate
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    out_dir = output.parent
    details: dict = {"output": str(output)}
    details["success_marker"] = (output / "_SUCCESS").exists()
    details["inputs_built"] = check_inputs_built(out_dir)

    if not details["success_marker"]:
        return False, False, details

    ht = read_table(output)
    key_ok = list(ht.key) == EXPECTED_KEY
    details["key"] = list(ht.key)
    details["key_ok"] = key_ok

    schema = hail_schema_to_dict(ht)
    details["row_count"] = ht.count()
    details["mapping_rate"] = spine_mapping_rate(ht)

    if reference is None:
        # --write-reference path never reaches here (write_reference() returns before
        # main() calls this); a None reference during grading means M4 cannot check the
        # contract's stored types, so it is reported as such rather than guessed at.
        details["contract"] = {"error": "no reference to check types against"}
        return False, False, details

    if "schema" not in reference or "contract_fields" not in reference:
        raise ValueError(
            "reference.json predates the gene-index contract fields; rerun the "
            "control with WRITE_REFERENCE=1"
        )
    ref_row = reference["schema"].get("row", {})
    missing_from_reference = [f for f in SPINE_FIELDS if f not in ref_row]
    if missing_from_reference:
        raise ValueError(
            "reference.json's stored schema is missing SPINE_FIELDS "
            f"{missing_from_reference}; rerun the control with WRITE_REFERENCE=1"
        )
    contract = schema_contract(
        schema.get("row", {}), ref_row, SPINE_FIELDS, EXPECTED_KEY
    )
    details["schema_contract_ok"] = contract["ok"]
    details["missing_fields"] = contract["missing_fields"]
    details["type_mismatches"] = contract["type_mismatches"]
    details["extra_fields"] = contract["extra_fields"]
    m4 = bool(
        details["success_marker"]
        and key_ok
        and contract["ok"]
        and details["inputs_built"]["ok"]
    )

    tolerances = reference.get("tolerances", {})
    mapping_rate_tol = tolerances.get("mapping_rate", MAPPING_RATE_TOLERANCE)

    digest = rows_digest(ht, reference["contract_fields"])
    details["content_sha256"] = digest
    content_ok = bool(
        details["row_count"] == reference["row_count"]
        and abs(details["mapping_rate"] - reference["mapping_rate"]) <= mapping_rate_tol
        and digest == reference["content_sha256"]
    )
    details["content_ok"] = content_ok
    return m4, content_ok, details


def grade_reported(
    notes: dict, row_count, mapping_rate, tol: float = REPORTED_RATE_TOLERANCE
) -> bool:
    """Self-consistency between the run notes' own ``genes:``/``mapping_rate:`` lines
    and the table actually graded -- independent of ``content_ok``. Never raises: a
    missing or malformed value (on either side) is simply ``False``."""
    if row_count is None or mapping_rate is None:
        return False
    if notes.get("genes") is None or notes.get("mapping_rate_raw") is None:
        return False
    reported_rate = parse_reported_rate(notes["mapping_rate_raw"])
    if reported_rate is None:
        return False
    return bool(
        notes["genes"] == row_count and abs(reported_rate - mapping_rate) <= tol
    )


def write_reference(args) -> None:
    from hvantk.algorithms.annotation.spine import SPINE_FIELDS, spine_mapping_rate
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    output = Path(args.out_dir) / GENE_INDEX_NAME
    ht = read_table(output)
    if list(ht.key) != EXPECTED_KEY:
        print(
            f"ERROR: --write-reference: {output} has key {list(ht.key)}, expected "
            f"{EXPECTED_KEY}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    missing = [f for f in SPINE_FIELDS if f not in ht.row]
    if missing:
        print(
            f"ERROR: --write-reference: {output} is missing SPINE_FIELDS: {missing}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    inputs = check_inputs_built(Path(args.out_dir))
    if not inputs["ok"]:
        print(
            f"ERROR: --write-reference: {args.out_dir} does not contain both input "
            f"builds (ensembl-gene:structure schema_id={ENSEMBL_SCHEMA_ID!r}, "
            f"hgnc:lookup schema_id={HGNC_SCHEMA_ID!r}); found: {inputs}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    schema = hail_schema_to_dict(ht)
    row_count = ht.count()
    mapping_rate = spine_mapping_rate(ht)
    # The contract digest is restricted to key + SPINE_FIELDS, sorted, so a later build
    # that legitimately carries extra fields beyond SPINE_FIELDS never costs it a
    # digest match (M4's contract check already covers the schema separately).
    contract_fields = sorted(set(EXPECTED_KEY) | set(SPINE_FIELDS))
    digest = rows_digest(ht, contract_fields)

    ref = {
        "schema": schema,
        "contract_fields": contract_fields,
        "row_count": row_count,
        "mapping_rate": mapping_rate,
        "content_sha256": digest,
        "tolerances": {
            "mapping_rate": MAPPING_RATE_TOLERANCE,
            "reported_mapping_rate": REPORTED_RATE_TOLERANCE,
        },
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
        f"wrote {args.reference}: row_count={row_count} mapping_rate={mapping_rate:.6f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="holds gene-index.ht")
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
    output = out_dir / GENE_INDEX_NAME
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

    # Reading the table, checking the two inputs, and grading the run notes against it
    # all happen in one try/except: a crash anywhere in here must still leave
    # verdict.json and receipt.json written (with the error recorded), not abort the
    # script before either file exists.
    try:
        m4, content_ok, gi_details = grade_gene_index(output, reference)
        reported_tol = reference.get("tolerances", {}).get(
            "reported_mapping_rate", REPORTED_RATE_TOLERANCE
        )
        reported_ok = grade_reported(
            notes,
            gi_details.get("row_count"),
            gi_details.get("mapping_rate"),
            reported_tol,
        )
    except Exception as exc:  # noqa: BLE001 -- a corrupt/wrong-type output or a
        # reference-schema problem must not crash the whole verdict; record it and
        # fail M4/M5 alone.
        m4, content_ok, reported_ok = False, False, False
        gi_details = {"output": str(output), "error": f"{type(exc).__name__}: {exc}"}

    m5 = bool(content_ok and reported_ok)
    passed = bool(m4 and m5)

    verdict = {
        "task": "gene-index",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "content_ok": content_ok,
        "reported_ok": reported_ok,
        "notes": notes,
        "slurm": slurm_details,
        "gene_index": gi_details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else out_dir
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    receipt = {
        "output": str(output),
        "content_sha256": gi_details.get("content_sha256"),
        "row_count": gi_details.get("row_count"),
        "mapping_rate": gi_details.get("mapping_rate"),
        "inputs_built": gi_details.get("inputs_built"),
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
