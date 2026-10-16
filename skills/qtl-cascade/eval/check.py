"""Deterministic check for the qtl-cascade eval task.

M1 from run-notes.md; M2 and M3 from an accounting file captured on the login node
(no array: task.md asks for at least four plain jobs — the three table builds and
the cascade run; eval/control.sh's own exemplar run also submits
summarize-cascade.py as a fifth, but that is not required of the agent under
evaluation). M4 = the per-tissue `cascade.ht` exists with `_SUCCESS`, is
keyed `(locus, alleles, gene_id)` and satisfies its field contract (read fresh from
THIS checkout's own schema below, never from `reference.json`); `gene_summary.tsv`
and `coloc_results.tsv` carry their expected columns; the report file exists and
contains the effective tissue label and the coloc section's exact heading from
`report.py`. The per-tissue directory and that label are resolved by
`resolve_tissue_dir`: `qtlcascade run --tissues <tissue>` writes under
`per_tissue/<_sanitize(tissue)>`, but a run given no `--tissues` filter at all
writes under `per_tissue/all` instead (identical numbers, for this eval's
single-tissue tables) — both are accepted. M5 has three parts:

- ``cascade_ok``: the cascade class counts equal the reference's, and a sha256 over
  `gene_summary.tsv`'s cascade-only contract columns (deliberately excluding the
  coloc overlay columns `coloc_max_h4`/`coloc_n_tissues` — see the negative control
  below) equals the reference's.
- ``coloc_ok``: the coloc gene set (a sha256 of the sorted gene_id list) equals the
  reference's, and for every one of those genes `n_variants` is exactly equal and
  every `H0`-`H4` is within ``H_TOLERANCE`` of the reference's, and the colocalised
  count (`H4 > 0.8`) is exactly equal.
- ``reported_ok``: the run notes' own `cascade_genes:`/`colocalised:` lines equal
  what this run's own outputs actually contain (a self-consistency check, not a
  comparison against the reference) — the gene_summary row count and the coloc
  colocalised count.

`cascade_ok` and `coloc_ok` are graded from different inputs on purpose: the negative
control rebuilds the pQTL table at a significant-pairs-only threshold and passes it
as BOTH the cascade and the coloc pQTL input. Cascade classification applies that
same threshold again internally (a no-op on an already-filtered table), so the
cascade output is unaffected and `cascade_ok` stays true; but coloc needs every
pair in the regional window, so feeding it significant-pairs-only input changes the
coloc results and `coloc_ok` goes false. Runs inside the image on a compute node
(Hail, for `cascade.ht`; pandas for the TSVs).

    python check.py --out-dir OUT --tissue Liver --run-notes NOTES.md \
        --reference reference.json --run-id PREFIX --expect-commit SHA \
        --expect-sif-sha256 HEX [--sacct-file FILE | --since now-2days | --no-slurm] \
        [--site-profile FILE] [--verdict-dir DIR] [--expect-fail]
    python check.py --out-dir OUT --tissue Liver --reference reference.json \
        --write-reference --expect-commit SHA --expect-sif-sha256 HEX --input-sha256 HEX

OUT is a `qtlcascade run -o` directory (holds `per_tissue/<tissue>/...` and
`qtlcascade_report.html`). The accounting file is the output of
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

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 45}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)
H4_THRESHOLD = 0.8
H_TOLERANCE = 1e-6

#: cascade.ht's field contract: a floor, read fresh from this checkout every time.
CASCADE_KEY = ["locus", "alleles", "gene_id"]
CASCADE_ROW = {
    "locus": "locus<GRCh38>",
    "alleles": "array<str>",
    "gene_id": "str",
    "eqtl_beta": "float64",
    "eqtl_se": "float64",
    "eqtl_pvalue": "float64",
    "pqtl_beta": "float64",
    "pqtl_se": "float64",
    "pqtl_pvalue": "float64",
    "cascade_class": "str",
    "attenuation_ratio": "float64",
    "tissue": "str",
}
#: gene_summary.tsv columns the task must produce (coloc ran, so the overlay columns
#: are expected too); CASCADE_CONTRACT_FIELDS is the subset `cascade_ok`'s digest
#: covers — deliberately excluding the coloc overlay (see module docstring).
GENE_SUMMARY_COLUMNS = [
    "gene_id",
    "n_eqtl_variants",
    "n_pqtl_variants",
    "n_concordant",
    "n_discordant",
    "best_eqtl_pvalue",
    "best_pqtl_pvalue",
    "has_complete_cascade",
    "tissues",
    "coloc_max_h4",
    "coloc_n_tissues",
]
CASCADE_CONTRACT_FIELDS = [
    "gene_id",
    "n_eqtl_variants",
    "n_pqtl_variants",
    "n_concordant",
    "n_discordant",
    "best_eqtl_pvalue",
    "best_pqtl_pvalue",
    "has_complete_cascade",
]
COLOC_COLUMNS = ["gene_id", "tissue", "H0", "H1", "H2", "H3", "H4", "n_variants"]
COLOC_H_COLUMNS = ["H0", "H1", "H2", "H3", "H4"]

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_CASCADE_GENES_NOTE = re.compile(r"^cascade_genes:\s*(\d+)\s*$", re.M)
_COLOCALISED_NOTE = re.compile(r"^colocalised:\s*(\d+)\s*$", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")


# --- run notes and Slurm accounting (M1-M3) -----------------------------------------


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {
            "present": False,
            "commit": None,
            "container": None,
            "jobs": [],
            "cascade_genes": None,
            "colocalised": None,
        }
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    g = _CASCADE_GENES_NOTE.search(text)
    m = _COLOCALISED_NOTE.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
        "jobs": _JOB.findall(text),
        "cascade_genes": int(g.group(1)) if g else None,
        "colocalised": int(m.group(1)) if m else None,
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


# --- the sha256 digest convention shared across this project's evals ---------------


def _norm(value):
    """Normalise one TSV-derived cell for a stable digest: pandas'
    `to_dict(orient="records")` boxes numpy scalars as native Python types in
    practice, but this stays correct even if a future/older pandas does not."""
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".10g")
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return None
    return value


def rows_digest(rows: list[dict]) -> str:
    """sha256 over the newline-joined, sorted list of normalised, sorted-key JSON
    rows. Order-insensitive. `rows` is the FULL table here (one eval window's worth
    of genes is small enough that sampling would add complexity for no benefit)."""
    strings = sorted(
        json.dumps(
            {k: _norm(v) for k, v in row.items()}, sort_keys=True, separators=(",", ":")
        )
        for row in rows
    )
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


def ids_digest(ids: list[str]) -> str:
    return hashlib.sha256("\n".join(sorted(set(ids))).encode()).hexdigest()


# --- reading the outputs -------------------------------------------------------------


def _fields_contract(built: dict, expected: dict) -> dict:
    missing = [name for name in expected if name not in built]
    mismatches = {
        name: {"expected": exp, "got": built.get(name)}
        for name, exp in expected.items()
        if name in built and built[name] != exp
    }
    extra = sorted(set(built) - set(expected))
    return {
        "ok": not missing and not mismatches,
        "missing_fields": missing,
        "type_mismatches": mismatches,
        "extra_fields": extra,
    }


def cascade_schema_and_counts(cascade_ht_path: Path) -> dict:
    from hvantk.core.utils.hail_context import init_hail
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    init_hail(quiet=True)
    import hail as hl

    ht = hl.read_table(str(cascade_ht_path))
    schema = hail_schema_to_dict(ht)
    key_ok = schema["key"] == CASCADE_KEY
    row_contract = _fields_contract(schema["row"], CASCADE_ROW)
    counts = dict(
        ht.group_by(ht.cascade_class)
        .aggregate(n=hl.agg.count())
        .to_pandas()
        .set_index("cascade_class")["n"]
    )
    class_counts = {str(k): int(v) for k, v in counts.items()}
    return {
        "key_ok": key_ok,
        "row_contract": row_contract,
        "schema_ok": key_ok and row_contract["ok"],
        "class_counts": class_counts,
        "n_rows": sum(class_counts.values()),
    }


def read_gene_summary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def read_coloc(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def resolve_tissue_dir(out_dir: Path, tissue: str) -> tuple[Path, str]:
    """The per-tissue directory ``qtlcascade run`` actually wrote to, and the
    label that run used for it.

    ``run --tissues <tissue>`` writes under ``per_tissue/<_sanitize(tissue)>``
    (``pipeline.py``'s own sanitiser, imported rather than re-derived, so this can
    never disagree with it); ``run`` given no ``--tissues`` filter at all writes
    under ``per_tissue/all`` instead, with the report's own tissue label also
    ``"all"`` — numerically identical to a single-tissue run for this eval's
    single-tissue tables, so both are accepted and the label that was actually
    used is returned for the report-content check.
    """
    from hvantk.algorithms.qtlcascade.pipeline import _sanitize

    named = out_dir / "per_tissue" / _sanitize(tissue)
    if named.exists():
        return named, tissue
    return out_dir / "per_tissue" / "all", "all"


def grade_output(out_dir: Path, tissue: str, reference: dict | None) -> tuple:
    """Returns (m4, m5, details)."""
    details: dict = {}
    tissue_dir, tissue_label = resolve_tissue_dir(out_dir, tissue)
    details["tissue_dir"] = str(tissue_dir)
    details["tissue_label"] = tissue_label
    cascade_ht = tissue_dir / "cascade.ht"
    gene_summary_tsv = tissue_dir / "gene_summary.tsv"
    coloc_tsv = tissue_dir / "coloc_results.tsv"
    report_html = out_dir / "qtlcascade_report.html"

    details["cascade_success"] = (cascade_ht / "_SUCCESS").exists()
    details["gene_summary_present"] = gene_summary_tsv.is_file()
    details["coloc_present"] = coloc_tsv.is_file()
    details["report_present"] = report_html.is_file()
    if not (
        details["cascade_success"]
        and details["gene_summary_present"]
        and details["coloc_present"]
        and details["report_present"]
    ):
        return False, False, details

    cascade_info = cascade_schema_and_counts(cascade_ht)
    details["cascade_schema_ok"] = cascade_info["schema_ok"]
    details["cascade_row_contract"] = cascade_info["row_contract"]
    details["class_counts"] = cascade_info["class_counts"]

    gene_df = read_gene_summary(gene_summary_tsv)
    missing_gene_cols = [c for c in GENE_SUMMARY_COLUMNS if c not in gene_df.columns]
    details["gene_summary_columns_ok"] = not missing_gene_cols
    details["gene_summary_missing_columns"] = missing_gene_cols
    details["n_cascade_genes"] = len(gene_df)

    coloc_df = read_coloc(coloc_tsv)
    missing_coloc_cols = [c for c in COLOC_COLUMNS if c not in coloc_df.columns]
    details["coloc_columns_ok"] = not missing_coloc_cols
    details["coloc_missing_columns"] = missing_coloc_cols
    details["coloc_n_tested"] = len(coloc_df)
    details["coloc_colocalised"] = int((coloc_df["H4"] > H4_THRESHOLD).sum())

    report_text = report_html.read_text(encoding="utf-8", errors="replace")
    details["report_has_tissue"] = tissue_label in report_text
    details["report_has_coloc_heading"] = "Colocalization Results" in report_text

    m4 = bool(
        cascade_info["schema_ok"]
        and details["gene_summary_columns_ok"]
        and details["coloc_columns_ok"]
        and details["report_has_tissue"]
        and details["report_has_coloc_heading"]
    )
    if not m4:
        return False, False, details

    # cascade_ok: class counts + a digest over the cascade-only contract columns
    # (never the coloc overlay columns — see module docstring).
    gene_rows = gene_df[CASCADE_CONTRACT_FIELDS].to_dict(orient="records")
    details["gene_summary_sha256"] = rows_digest(gene_rows)

    # coloc_ok inputs, independent of cascade_ok's.
    coloc_by_gene = {
        str(r["gene_id"]): {
            "n_variants": int(r["n_variants"]),
            **{h: float(r[h]) for h in COLOC_H_COLUMNS},
        }
        for _, r in coloc_df.iterrows()
    }
    details["coloc_gene_ids_sha256"] = ids_digest(list(coloc_by_gene))

    if reference is None:
        return m4, False, details

    details["class_counts_match"] = details["class_counts"] == reference["class_counts"]
    details["gene_summary_sha256_match"] = (
        details["gene_summary_sha256"] == reference["gene_summary_sha256"]
    )
    cascade_ok = bool(
        details["class_counts_match"] and details["gene_summary_sha256_match"]
    )
    details["cascade_ok"] = cascade_ok

    ref_coloc = reference["coloc"]
    gene_ids_match = details["coloc_gene_ids_sha256"] == ref_coloc["gene_ids_sha256"]
    details["coloc_gene_ids_match"] = gene_ids_match
    per_gene_ok = gene_ids_match
    if gene_ids_match:
        for gene_id, ref_vals in ref_coloc["genes"].items():
            got = coloc_by_gene.get(gene_id)
            if got is None or got["n_variants"] != ref_vals["n_variants"]:
                per_gene_ok = False
                break
            if any(abs(got[h] - ref_vals[h]) > H_TOLERANCE for h in COLOC_H_COLUMNS):
                per_gene_ok = False
                break
    details["coloc_per_gene_ok"] = per_gene_ok
    details["coloc_colocalised_match"] = (
        details["coloc_colocalised"] == ref_coloc["colocalised_count"]
    )
    coloc_ok = bool(per_gene_ok and details["coloc_colocalised_match"])
    details["coloc_ok"] = coloc_ok

    m5 = bool(cascade_ok and coloc_ok)
    return m4, m5, details


def write_reference(args) -> None:
    out_dir = Path(args.out_dir)
    m4, _, details = grade_output(out_dir, args.tissue, reference=None)
    if not m4:
        print(
            "ERROR: --write-reference refuses because M4 is not true for this "
            f"output; details: {json.dumps(details, indent=2, default=str)}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    tissue_dir, _ = resolve_tissue_dir(out_dir, args.tissue)
    coloc_df = read_coloc(tissue_dir / "coloc_results.tsv")
    genes = {
        str(r["gene_id"]): {
            "n_variants": int(r["n_variants"]),
            **{h: float(r[h]) for h in COLOC_H_COLUMNS},
        }
        for _, r in coloc_df.iterrows()
    }
    ref = {
        "tissue": args.tissue,
        "class_counts": details["class_counts"],
        "n_cascade_genes": details["n_cascade_genes"],
        "gene_summary_sha256": details["gene_summary_sha256"],
        "coloc": {
            "gene_ids_sha256": details["coloc_gene_ids_sha256"],
            "colocalised_count": details["coloc_colocalised"],
            "genes": genes,
        },
        "tolerance": {"H": H_TOLERANCE},
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
        f"wrote {args.reference}: class_counts={details['class_counts']} "
        f"n_cascade_genes={details['n_cascade_genes']} "
        f"coloc genes={len(genes)} colocalised={details['coloc_colocalised']}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="a qtlcascade run's -o directory")
    ap.add_argument("--tissue", required=True)
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

    try:
        # Loaded inside the try, alongside grade_output: a missing or corrupt
        # reference.json must still leave a verdict.json on disk (M4 can still be
        # computed; M5 cannot), never crash before anything is written.
        reference = json.loads(Path(args.reference).read_text())
        m4, m5, details = grade_output(out_dir, args.tissue, reference)
    except Exception as exc:  # noqa: BLE001 -- a corrupt/missing output or
        # reference must not crash the verdict; record it and fail M4/M5, same
        # pattern as this project's other evals.
        m4, m5 = False, False
        details = {"error": f"{type(exc).__name__}: {exc}"}

    reported_ok = False
    if m4:
        reported_ok = bool(
            notes["cascade_genes"] == details.get("n_cascade_genes")
            and notes["colocalised"] == details.get("coloc_colocalised")
        )
    details["reported_ok"] = reported_ok
    m5 = bool(m5 and reported_ok)
    passed = bool(m4 and m5)

    verdict = {
        "task": "qtl-cascade",
        "tissue": args.tissue,
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "notes": notes,
        "slurm": slurm_details,
        "output": details,
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = Path(args.verdict_dir) if args.verdict_dir else out_dir
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    try:
        reference_sha256 = hashlib.sha256(Path(args.reference).read_bytes()).hexdigest()
    except OSError:
        reference_sha256 = None
    receipt = {
        "tissue": args.tissue,
        "pass": passed,
        "gene_summary_sha256": details.get("gene_summary_sha256"),
        "coloc_gene_ids_sha256": details.get("coloc_gene_ids_sha256"),
        "reference_sha256": reference_sha256,
        "checked_at": verdict["checked_at"],
    }
    (verdict_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(
        json.dumps(
            {k: verdict[k] for k in ("M1", "M2", "M3", "M4", "M5", "pass")}
            | {
                "cascade_ok": details.get("cascade_ok"),
                "coloc_ok": details.get("coloc_ok"),
                "reported_ok": details.get("reported_ok"),
            }
        )
    )
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
