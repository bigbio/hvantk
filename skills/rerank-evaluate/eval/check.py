"""Deterministic check for the rerank-evaluate eval task.

Grades one run in five stages (M1 to M5, see skills/_conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json``, and exits 0 on pass. No Hail: it reads
``config.yaml`` (PyYAML) and ``ranked.tsv`` (pandas + scikit-learn) only, so it runs
through the image on a login node (a plain container exec, no Slurm job needed for the
check itself) exactly like ``ancestry-inference``'s own check.py. Slurm accounting is
read from a file captured on the login node (``--sacct-file``), because the container
has no Slurm client; without a file it tries ``sacct`` and records M2 and M3 as null
when that is unavailable or fails.

M4 is config-and-shape only (``config.yaml`` parses with exactly the allowed keys, one
feature axis declared; ``ranked.tsv`` carries the nine columns ``hvantk rerank``
always emits, with at least one label-positive and one label-negative row). M5 is
five independent sub-checks, ANDed together:

* ``universe_ok`` -- the digest of the ranked table's own ``gene`` column (sorted)
  equals the reference's. This is what catches the negative control: a universe
  built from the whole gene index instead of the cohort's tested genes changes this
  digest while changing nothing else about how the run was configured.
* ``labels_ok`` -- the digest of the genes with ``y == 1`` equals the reference's.
* ``audit_ok`` -- the digest of ``(gene, flag, flag_reason)`` triples equals the
  reference's.
* ``metrics_ok`` -- ``roc_auc_score``/``average_precision_score`` recomputed directly
  from ``ranked.tsv`` (never trusted from ``run-notes.md`` or a prior receipt) are
  within 0.002 of the reference's.
* ``reported_ok`` -- the ``auc:`` line in ``run-notes.md`` is within 0.0005 of the
  AUC this check just recomputed, i.e. the run notes were not hand-typed or copied
  from a different run.

One dataset's digest never depends on row order (every digest below sorts first), and
every digest is restricted to exactly the fields it names -- never the live table's
full column set -- so extra columns a future config might add cannot flip a digest
that has nothing to do with them.

    python check.py --config config.yaml --ranked ranked.tsv --run-notes NOTES.md \\
        --reference reference.json --run-id PREFIX --expect-commit SHA \\
        --expect-sif-sha256 HEX [--sacct-file FILE | --since now-2days | --no-slurm] \\
        [--site-profile FILE] [--verdict-dir DIR] [--expect-fail]
    python check.py --config config.yaml --ranked ranked.tsv --reference reference.json \\
        --write-reference --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

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

BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 45}
SUCCESS_STATES = ("COMPLETED",)
AUC_TOLERANCE = 0.002
AP_TOLERANCE = 0.002
REPORTED_TOLERANCE = 0.0005
TOP_DECILE_PCT = 90.0
REQUIRED_COLUMNS = (
    "gene",
    "prior_stat",
    "score",
    "score_percentile",
    "tier",
    "verdict",
    "flag",
    "flag_reason",
    "y",
)
ALLOWED_CONFIG_KEYS = {"name", "cohort", "features", "labels", "min_label_coverage"}
REQUIRED_CONFIG_KEYS = {"cohort", "features", "labels"}

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_AUC_LINE = re.compile(r"^auc:\s*([0-9.]+)\s*$", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


# --- run notes and Slurm accounting (M1-M3); copied, not imported, per
# --- skills/_conventions/SKILL.md section 6: every eval/check.py is self-contained ---


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {
            "present": False,
            "commit": None,
            "container": None,
            "jobs": [],
            "auc": None,
        }
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    a = _AUC_LINE.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
        "auc": float(a.group(1)) if a else None,
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


SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"


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
    m3 = (
        bool(notes["jobs"])
        and any(j["state"].startswith(SUCCESS_STATES) for j in jobs)
        and all(n in all_ids for n in notes["jobs"])
    )
    return m2, m3, {"jobs": jobs, "within_budget": within}


# --- order-insensitive digests, following skills/_conventions/SKILL.md section 5's
# --- rows_digest convention: sha256 of the sorted, newline-joined json.dumps of each
# --- normalised row; floats as format(x, ".10g") so nan/inf hash as stable strings ---


def _norm(x):
    if isinstance(x, float):
        return format(x, ".10g")
    if isinstance(x, dict):
        return {k: _norm(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_norm(v) for v in x]
    return x


def _digest(rows) -> str:
    strings = sorted(
        json.dumps(_norm(r), sort_keys=True, separators=(",", ":")) for r in rows
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


def _to_bool(series):
    """A bool Series, whether pandas inferred ``flag`` as bool or left it as text."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1"})


# --- reading the ranked table and the config -----------------------------------------


def load_ranked(path: Path):
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    missing = sorted(set(REQUIRED_COLUMNS) - set(df.columns))
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}")
    df["flag"] = _to_bool(df["flag"])
    # Pandas' default NA tokens turn an empty flag_reason ("" for a non-flagged gene,
    # exactly what hvantk rerank writes) into NaN on re-read; restore it so audit_ok's
    # digest sees the same value the run itself produced, not a round-trip artifact.
    df["flag_reason"] = df["flag_reason"].fillna("")
    df["gene"] = df["gene"].astype(str)
    return df


def grade_config(path: Path) -> tuple[bool, dict]:
    import yaml

    details: dict = {"path": str(path)}
    if not path.is_file():
        details["error"] = "config not found"
        return False, details
    try:
        doc = yaml.safe_load(path.read_text())
    except yaml.YAMLError as exc:
        details["error"] = f"invalid YAML: {exc}"
        return False, details
    if not isinstance(doc, dict):
        details["error"] = f"expected a YAML mapping, got {type(doc).__name__}"
        return False, details
    keys = set(doc)
    details["keys"] = sorted(keys)
    unknown = sorted(keys - ALLOWED_CONFIG_KEYS)
    missing_required = sorted(REQUIRED_CONFIG_KEYS - keys)
    details["unknown_keys"] = unknown
    details["missing_required_keys"] = missing_required
    features = doc.get("features")
    details["n_features"] = len(features) if isinstance(features, list) else None
    ok = (
        not unknown
        and not missing_required
        and isinstance(features, list)
        and len(features) == 1
    )
    return ok, details


def grade_ranked_shape(df) -> tuple[bool, dict]:
    details = {
        "columns": list(df.columns),
        "n_rows": int(len(df)),
        "n_positive": int((df["y"] == 1).sum()),
        "n_negative": int((df["y"] == 0).sum()),
    }
    ok = (
        set(df.columns) == set(REQUIRED_COLUMNS)
        and details["n_positive"] > 0
        and details["n_negative"] > 0
    )
    return ok, details


def compute_metrics(df) -> dict:
    from sklearn.metrics import average_precision_score, roc_auc_score

    scored = df[df["score"].notna()]
    y = scored["y"].astype(int)
    score = scored["score"].astype(float)
    auc = float(roc_auc_score(y, score))
    ap = float(average_precision_score(y, score))
    top_decile = scored[scored["score_percentile"].astype(float) >= TOP_DECILE_PCT]
    tier_counts = {str(k): int(v) for k, v in scored["tier"].value_counts().items()}
    verdict_counts = {
        str(k): int(v) for k, v in scored["verdict"].value_counts().items()
    }
    flag_counts = {
        str(k): int(v)
        for k, v in scored.loc[scored["flag"], "flag_reason"].value_counts().items()
    }
    return {
        "n_genes": int(len(scored)),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "auc": auc,
        "average_precision": ap,
        "n_top_decile_positive": int((top_decile["y"] == 1).sum()),
        "tier_counts": tier_counts,
        "verdict_counts": verdict_counts,
        "flag_counts": flag_counts,
    }


def grade_content(df, reference: dict | None) -> tuple[dict, dict]:
    """(sub_checks, details) for universe/labels/audit/metrics -- M5 minus reported_ok."""
    universe_sha256 = _digest(sorted(df["gene"]))
    labels_sha256 = _digest(sorted(df.loc[df["y"] == 1, "gene"]))
    audit_sha256 = _digest(
        {"gene": g, "flag": bool(f), "flag_reason": r}
        for g, f, r in zip(df["gene"], df["flag"], df["flag_reason"])
    )
    metrics = compute_metrics(df)
    details = {
        "universe_sha256": universe_sha256,
        "labels_sha256": labels_sha256,
        "audit_sha256": audit_sha256,
        "metrics": metrics,
    }
    if reference is None:
        return {}, details
    checks = {
        "universe_ok": universe_sha256 == reference["universe_sha256"],
        "labels_ok": labels_sha256 == reference["labels_sha256"],
        "audit_ok": audit_sha256 == reference["audit_sha256"],
        "metrics_ok": (
            abs(metrics["auc"] - reference["auc"]) <= AUC_TOLERANCE
            and abs(metrics["average_precision"] - reference["average_precision"])
            <= AP_TOLERANCE
        ),
    }
    return checks, details


def write_reference(args) -> None:
    config_ok, config_details = grade_config(Path(args.config))
    if not config_ok:
        print(
            f"ERROR: --write-reference: {args.config} does not satisfy the config "
            f"contract: {json.dumps(config_details)}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    df = load_ranked(Path(args.ranked))
    shape_ok, shape_details = grade_ranked_shape(df)
    if not shape_ok:
        print(
            f"ERROR: --write-reference: {args.ranked} does not satisfy the ranked-table "
            f"contract: {json.dumps(shape_details)}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    _, details = grade_content(df, reference=None)
    metrics = details["metrics"]
    ref = {
        "universe_sha256": details["universe_sha256"],
        "labels_sha256": details["labels_sha256"],
        "audit_sha256": details["audit_sha256"],
        "auc": metrics["auc"],
        "average_precision": metrics["average_precision"],
        "n_genes": metrics["n_genes"],
        "n_positive": metrics["n_positive"],
        "n_negative": metrics["n_negative"],
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
        f"wrote {args.reference}: n_genes={metrics['n_genes']} "
        f"n_positive={metrics['n_positive']} auc={metrics['auc']:.4f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--ranked", required=True)
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
    ap.add_argument("--verdict-dir", default=None)
    ap.add_argument("--expect-fail", action="store_true")
    ap.add_argument("--write-reference", action="store_true")
    args = ap.parse_args()

    if args.write_reference:
        write_reference(args)
        return 0
    for required in ("run_id", "expect_commit", "expect_sif_sha256"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required when grading")

    ranked_path = Path(args.ranked)
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

    config_ok, config_details = grade_config(Path(args.config))
    try:
        df = load_ranked(ranked_path)
        shape_ok, shape_details = grade_ranked_shape(df)
        content_checks, content_details = grade_content(df, reference)
        reported_ok = bool(
            notes["auc"] is not None
            and abs(notes["auc"] - content_details["metrics"]["auc"])
            <= REPORTED_TOLERANCE
        )
        content_checks["reported_ok"] = reported_ok
        m5 = bool(content_checks) and all(content_checks.values())
    except Exception as exc:  # noqa: BLE001 -- a corrupt or missing ranked table
        # must not crash the whole verdict; record it and fail M5 (M4 still reflects
        # the config/shape checks that did succeed).
        shape_ok, shape_details = False, {"error": f"{type(exc).__name__}: {exc}"}
        content_checks, content_details = {}, {"error": f"{type(exc).__name__}: {exc}"}
        reported_ok = False
        m5 = False
    m4 = bool(config_ok and shape_ok)
    passed = bool(m4 and m5)

    verdict = {
        "task": "rerank-evaluate",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "config": config_details,
        "ranked": shape_details,
        "content_checks": content_checks,
        "content": content_details,
        "notes": notes,
        "slurm": slurm_details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else ranked_path.parent
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    receipt = {
        "ranked": str(ranked_path),
        "pass": passed,
        "ranked_sha256": (
            hashlib.sha256(ranked_path.read_bytes()).hexdigest()
            if ranked_path.is_file()
            else None
        ),
        "reference_sha256": hashlib.sha256(
            Path(args.reference).read_bytes()
        ).hexdigest(),
        "checked_at": verdict["checked_at"],
    }
    (verdict_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(
        json.dumps(
            {
                **{k: verdict[k] for k in ("M1", "M2", "M3", "M4", "M5", "pass")},
                "content_checks": content_checks,
            }
        )
    )
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
