"""Deterministic check for the drift-triage eval task.

Grades one run in five stages (M1 to M5, see skills/workflow-conventions/SKILL.md section 5),
writes ``verdict.json`` and ``receipt.json``, and exits 0 on pass. Runs inside the
container on a compute node (it needs Hail to read the built table and, separately,
to run a small subprocess build -- see `grade_repair_effective`). Slurm accounting is
read from a file captured on the login node (``--sacct-file``), because the container
has no Slurm client; without a file it tries ``sacct`` and records M2 and M3 as null
when that is unavailable or fails.

This grades THREE things that can vary independently: the built table's content (M4,
and the ``content_ok`` half of M5; ``grade_content`` is called inside a try/except in
``main()``, so a corrupt or wrong-type output records ``error`` in its own details and
fails M4/``content_ok`` alone), whether the GRADED checkout's own code can actually
reproduce the repair from scratch (``repair_effective``, part of M5, proven by two
small isolated subprocess builds -- a big ``--output`` table could in principle have
been produced some other way; one build is the header-mutated fixture, the other is
the fixture UNCHANGED, because a correct repair must keep the old source name working
too, not just add the new one), and the state of the checkout that produced it -- its
fingerprint, its ledger row, and whether its committed snapshots were left alone (the
rest of M5). ``--repo`` names the checkout to read the latter two from; it may be a
copy of the real checkout, which is exactly how the positive and negative controls
exercise it.

The schema check is a CONTRACT, not an equality, against the committed
``hvantk/skills/clingen/tests/snapshots/schema.json`` -- read from THIS script's own
checkout (never ``--repo``, so a tampered copy's own snapshot cannot loosen the
contract): the key must match exactly, and every field the snapshot declares must be
present with the same type, but a full build may legitimately carry MORE fields than
a small fixture-derived snapshot does (the gencc fixture, for one, ships 14 of
upstream's 31 columns). ``content_ok`` respects the same floor rather than tightening
it back into an equality: the rows digest is computed only over ``contract_fields``
(the committed snapshot's key fields plus its row field names, stored in
``reference.json`` at ``--write-reference`` time) -- row count and digest equal is the
whole of ``content_ok``; schema equality against the reference is recorded for
information only, since M4's contract above already covers the schema.

    python check.py --output OUT.ht --repo REPO --reference reference.json \
        [--verdict-dir DIR] --run-notes NOTES.md --run-id PREFIX \
        --expect-commit SHA --expect-sif-sha256 HEX \
        [--sacct-file FILE | --since now-2days | --no-slurm] [--site-profile FILE] \
        [--expect-fail]
    python check.py --output OUT.ht --repo REPO --reference reference.json \
        --write-reference --staged-file MUTATED \
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
import os
import re
import subprocess
import sys
from pathlib import Path

#: This script's own checkout root (skills/drift-triage/eval/check.py -> eval ->
#: drift-triage -> skills -> repo root), deliberately NEVER ``--repo``: the schema
#: CONTRACT and the fixture `grade_repair_effective` builds from must come from a
#: checkout this script trusts, not from whatever copy is under grading.
OWN_REPO = Path(__file__).resolve().parents[3]

DATASET = "clingen:gene-disease"
SCHEMA_ID = "clingen-gene-disease-v1"
FILENAME = "Clingen-Gene-Disease-Summary.csv"
LEDGER = "hvantk/resources/drift_ledger.json"
FINGERPRINT = "hvantk/skills/clingen/tests/drift_fingerprint.json"
SNAPSHOTS = "hvantk/skills/clingen/tests/snapshots"
FIXTURE = "hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv"
BUDGET = {"max_cpus": 4, "max_mem_gb": 16, "max_minutes": 30}
SACCT_FIELDS = "JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
SUCCESS_STATES = ("COMPLETED",)
N_SAMPLE_KEYS = 12
REPAIR_TIMEOUT_S = 600

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)
_JOB = re.compile(r"^job:\s*(\d+(?:_\d+)?)", re.M)
_MEM = re.compile(r"^([0-9.]+)([KMGT])([nc]?)$")
_PROFILE_PARTITIONS = ("HVANTK_CPU_PARTITION", "HVANTK_FAT_PARTITION")

# Run in its own subprocess (never in-process, and always before this process starts
# Hail), with --repo on PYTHONPATH: builds BOTH the committed ClinGen fixture
# header-mutated the same way the eval's staged file is, AND the fixture unmodified --
# proving the GRADED checkout's own code (not whatever produced --output) can build
# the repaired layout AND still builds the OLD source name, which a correct repair
# must keep (the plugin's own tests build from the unmutated fixture and would break
# if the old rename were replaced instead of extended). Asserts the anchor was
# replaced exactly once before building either layout. Reads the fixture path and
# nothing else from the environment, so quoting the path into this source is never a
# concern.
_REPAIR_SNIPPET = """
import json
import os
import sys
import tempfile
from pathlib import Path

fixture = Path(os.environ["DRIFT_TRIAGE_FIXTURE"])
anchor = '"GENE ID (HGNC)"'
replacement = '"HGNC ID"'

try:
    original_text = fixture.read_text()
    if original_text.count(anchor) != 1:
        print(json.dumps({"ok": False, "reason": "mutation did not apply"}))
        sys.exit(1)
    mutated_text = original_text.replace(anchor, replacement, 1)

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)

    from hvantk.core.plugin.loader import get_registry
    from hvantk.tests._snapshot_utils import phase_b_snapshot_adapter

    spec = get_registry().get_dataset("clingen:gene-disease")
    adapter = phase_b_snapshot_adapter(spec.builder, "clingen:gene-disease")

    def build(contents):
        tmp_dir = tempfile.mkdtemp(prefix="drift-triage-repair-")
        raw = Path(tmp_dir) / fixture.name
        raw.write_text(contents)
        out_path = str(Path(tmp_dir) / "out.ht")
        ht = adapter(input_path=tmp_dir, output_path=out_path)
        return list(ht.key), ht.count()

    # BOTH layouts, built by the GRADED checkout's own code: the header-mutated
    # fixture (the new source name) and the unmodified committed fixture (the old
    # source name) -- a repair that replaced the old rename instead of adding to it
    # would build the first and fail the second.
    key_new, n_rows_new = build(mutated_text)
    key_old, n_rows_old = build(original_text)
    print(
        json.dumps(
            {
                "ok": True,
                "key_new": key_new,
                "n_rows_new": n_rows_new,
                "key_old": key_old,
                "n_rows_old": n_rows_old,
            }
        )
    )
except Exception as exc:
    print(json.dumps({"ok": False, "error": "%s: %s" % (type(exc).__name__, exc)}))
    sys.exit(1)
"""

# --- the common pattern, shared (near-verbatim) with skills/hpc-hail-job/eval/check.py ---


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
    m3 = (
        bool(notes["jobs"])
        and any(j["state"].startswith(SUCCESS_STATES) for j in jobs)
        and all(n in all_ids for n in notes["jobs"])
    )
    return m2, m3, {"jobs": jobs, "within_budget": within}


# --- robustness helpers (fix D): naive timestamps read as UTC; a fingerprint or
#     ledger that is not a JSON object -- at any nesting level this code reads --
#     counts as absent rather than raising. ---


def _parse_iso(s: str | None) -> dt.datetime | None:
    if not s:
        return None
    try:
        parsed = dt.datetime.fromisoformat(str(s).replace("Z", "+00:00"))
    except (ValueError, AttributeError, TypeError):
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed


def _read_json_object(path: Path) -> dict | None:
    """``path`` read as a JSON object, or None: missing, unreadable, unparseable, or
    parsing to something other than a JSON object (a list, a string, ...) are all
    "absent", never a crash -- a fingerprint or ledger a repair never touched is
    exactly this shape in a freshly-copied, never-regenerated checkout."""
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _nested_get(d, *keys):
    """``d.get(keys[0]).get(keys[1])...``, returning None the moment anything along
    the way is not a dict, instead of raising ``AttributeError``."""
    cur = d
    for k in keys:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(k)
    return cur


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


# --- dataset/content specific (fix A: no SAMPLE_KEYS; fix B: schema contract) ---


def schema_contract(schema: dict, committed: dict) -> tuple[bool, dict]:
    """True if ``schema`` satisfies ``committed``'s CONTRACT: the same key, and every
    row field ``committed`` declares present with the same type. Extra fields in
    ``schema`` beyond what ``committed`` declares are fine -- a committed snapshot is
    a FLOOR, not an exact description of a full build (the gencc fixture carries 14
    of upstream's 31 columns, and every builder imports and renames every column
    present, so a full build legitimately carries more).
    """
    details: dict = {}
    key_ok = schema.get("key") == committed.get("key")
    details["key_ok"] = key_ok
    row = schema.get("row") if isinstance(schema.get("row"), dict) else {}
    committed_row = (
        committed.get("row") if isinstance(committed.get("row"), dict) else {}
    )
    missing = [f for f in committed_row if f not in row]
    wrong_type = [
        {"field": f, "expected": committed_row[f], "got": row.get(f)}
        for f in committed_row
        if f in row and row[f] != committed_row[f]
    ]
    details["missing_fields"] = missing
    details["wrong_type_fields"] = wrong_type
    ok = bool(key_ok and not missing and not wrong_type)
    return ok, details


def select_sample_keys(ht) -> list[dict]:
    """``N_SAMPLE_KEYS`` keys chosen deterministically from the table's row order,
    independent of content -- a dataset's hand-authored SAMPLE_KEYS is not guaranteed
    to exist in a real build (confirmed on the cluster: ``collect_sample_rows`` raised
    ``KeyError`` for a committed ClinGen key against a real build), so this selects
    from whatever table is actually in front of it.
    """
    import hail as hl
    from hvantk.tests._snapshot_utils import _to_jsonable

    n = ht.count()
    if n == 0:
        return []
    positions = sorted({(i * n) // N_SAMPLE_KEYS for i in range(N_SAMPLE_KEYS)})
    indexed = ht.add_index("_row_index")
    chosen = indexed.filter(
        hl.literal(set(positions), dtype=hl.tset(hl.tint64)).contains(
            indexed["_row_index"]
        )
    )
    chosen = chosen.select()  # Table.select() with NO args: keeps exactly the key
    # fields. select(*key_field_names) would raise ExpressionException ("cannot
    # overwrite key field ... with annotate, select or drop") since the key fields
    # are already there implicitly -- naming them again is an overwrite, not a select.
    rows = chosen.collect()
    return [{k: _to_jsonable(v) for k, v in dict(row).items()} for row in rows]


def rows_digest(ht, keys: list[dict], fields: list[str]) -> str:
    """sha256 over every row of ``ht`` whose key is in ``keys``, projected down to
    *fields* first (the dataset's CONTRACT: its committed snapshot's key fields plus
    its row field names -- never the live schema's full field list, so a build that
    legitimately carries more columns than the fixture-derived snapshot never fails
    the digest over that alone): order-insensitive (the per-row JSON strings are
    collected then SORTED before joining) and duplicate-key-safe (a semi-join against
    ``ht``'s own unique keys contributes each matching row once, however many times
    its key appears in ``keys``).
    """
    import hail as hl
    from hvantk.tests._snapshot_utils import _jsonable_to_hail_python, _to_jsonable

    if not keys:
        return hashlib.sha256(b"").hexdigest()
    key_field_names = list(ht.key)
    key_dtype = ht.key.dtype
    key_rows = [
        {
            name: _jsonable_to_hail_python(k[name], key_dtype[name])
            for name in key_field_names
        }
        for k in keys
    ]
    requested = hl.Table.parallelize(key_rows, schema=key_dtype, key=key_field_names)
    matched = ht.semi_join(requested)
    # Table.select() rejects key-field names; key fields stay on the row regardless,
    # so only the non-key contract fields are named here.
    non_key_fields = [f for f in fields if f not in key_field_names]
    collected = matched.select(*non_key_fields).collect()
    strings = sorted(
        json.dumps(
            {k: _to_jsonable(v) for k, v in dict(row).items()},
            sort_keys=True,
            separators=(",", ":"),
        )
        for row in collected
    )
    return hashlib.sha256("\n".join(strings).encode("utf-8")).hexdigest()


def read_table(output: Path):
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    return hl.read_table(str(output))


def grade_content(
    output: Path, reference: dict, committed_schema: dict
) -> tuple[bool, bool, dict]:
    """M4 (the schema CONTRACT) and the content half of M5 (exact reproducibility
    against the reference's own build), from --output and the reference alone."""
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

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

    contract_ok, contract_details = schema_contract(schema, committed_schema)
    details["schema_contract"] = contract_details
    m4 = bool(prov.get("schema_id") == SCHEMA_ID and contract_ok)

    if "contract_fields" not in reference:
        raise ValueError(
            "reference predates contract_fields; rerun the control with "
            "WRITE_REFERENCE=1"
        )
    # Informational only: M4's contract above already covers the schema, and a full
    # build may legitimately carry more fields than the fixture-derived snapshot the
    # reference's digest is restricted to, so schema equality is not part of
    # content_ok below.
    details["schema_matches_reference_exactly"] = schema == reference["schema"]
    details["row_count"] = ht.count()
    row_count_ok = details["row_count"] == reference["row_count"]
    rows_sha = rows_digest(ht, reference["sample_keys"], reference["contract_fields"])
    details["rows_sha256"] = rows_sha
    rows_ok = rows_sha == reference["rows_sha256"]
    details["row_count_ok"] = row_count_ok
    details["rows_ok"] = rows_ok
    content_ok = bool(row_count_ok and rows_ok)
    return m4, content_ok, details


def grade_repair_effective(repo: Path, contract_key) -> tuple[bool, dict]:
    """Prove ``repo``'s OWN code builds BOTH the mutated and the unmutated layout,
    independent of whatever produced ``--output``: two fresh, tiny subprocess builds
    of the committed ClinGen fixture (header-mutated, and as committed), with ``repo``
    on PYTHONPATH. A repair that renamed the column instead of adding the new name
    alongside the old one would build the first and fail the second -- exactly the
    regression the plugin's own tests (run from the unmutated fixture) exist to catch.
    Always run before this process starts its own Hail/JVM, so only one JVM is ever up
    at a time.

    ``cwd=repo`` matters as much as ``PYTHONPATH=repo`` here: ``python -c`` puts the
    *current working directory* at the front of ``sys.path`` (exactly like
    ``python -m``), ahead of ``PYTHONPATH``. Without pinning it, the subprocess
    inherits whatever directory this check happened to be run from, and if that
    directory also has an ``hvantk/`` package next to it (the real checkout, say),
    ``import hvantk`` resolves there instead of into ``repo`` -- silently grading the
    wrong code. Pointing ``cwd`` at ``repo`` makes the two agree instead of race.
    """
    fixture = OWN_REPO / FIXTURE
    env = dict(os.environ)
    env["PYTHONPATH"] = str(repo)
    env["PYSPARK_SUBMIT_ARGS"] = "--master local[1] --driver-memory 1500m pyspark-shell"
    env["DRIFT_TRIAGE_FIXTURE"] = str(fixture)
    details: dict = {}
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _REPAIR_SNIPPET],
            env=env,
            cwd=str(repo),
            capture_output=True,
            text=True,
            timeout=REPAIR_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired as exc:
        details["error"] = f"timed out after {REPAIR_TIMEOUT_S}s"
        details["stderr_tail"] = (exc.stderr or "")[-2000:]
        return False, details

    details["returncode"] = proc.returncode
    details["stderr_tail"] = proc.stderr[-2000:] if proc.stderr else ""
    last_line = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
    try:
        payload = json.loads(last_line)
    except ValueError:
        details["parse_error"] = True
        return False, details
    details["payload"] = payload
    ok = bool(
        payload.get("ok") is True
        and payload.get("key_new") == contract_key
        and isinstance(payload.get("n_rows_new"), int)
        and payload.get("n_rows_new") >= 1
        and payload.get("key_old") == contract_key
        and isinstance(payload.get("n_rows_old"), int)
        and payload.get("n_rows_old") >= 1
    )
    return ok, details


def grade_bookkeeping(repo: Path, reference: dict) -> tuple[bool, bool, bool, dict]:
    """The rest of M5: snapshots untouched, fingerprint regenerated, ledger row
    recorded. Everything read from ``repo`` -- unlike the schema contract -- because
    this IS a question about what ``repo`` itself recorded.
    """
    details: dict = {}

    schema_sha = sha256_file(repo / SNAPSHOTS / "schema.json")
    rows_sha_snap = sha256_file(repo / SNAPSHOTS / "sample_rows.json")
    details["schema_snapshot_sha256"] = schema_sha
    details["row_snapshot_sha256"] = rows_sha_snap
    snapshots_unchanged = bool(
        schema_sha == reference["snapshots"]["schema_sha256"]
        and rows_sha_snap == reference["snapshots"]["sample_rows_sha256"]
    )
    details["snapshots_unchanged"] = snapshots_unchanged

    fp = _read_json_object(repo / FINGERPRINT)
    details["fingerprint_parses"] = fp is not None
    expected_fp = reference["expected_fingerprint"]
    fingerprint_regenerated = bool(
        fp is not None
        and _nested_get(fp, "headers", FILENAME) == expected_fp["headers"]
        and _nested_get(fp, "checksums", FILENAME) == expected_fp["checksum"]
    )
    details["fingerprint_regenerated"] = fingerprint_regenerated

    ledger = _read_json_object(repo / LEDGER) or {}
    row = ledger.get(DATASET)
    details["ledger_row"] = row
    ledger_row_ok = False
    if isinstance(row, dict):
        last_change = _parse_iso(row.get("last_upstream_change"))
        rebuilt_at = _parse_iso(row.get("rebuilt_at"))
        committed = reference.get("committed_ledger_row") or {}
        committed_change = _parse_iso(
            committed.get("last_upstream_change")
            if isinstance(committed, dict)
            else None
        )
        accepted_in = row.get("accepted_in")
        ledger_row_ok = bool(
            row.get("signal") == "schema"
            and isinstance(accepted_in, str)
            and accepted_in.strip()
            and last_change is not None
            and (committed_change is None or last_change > committed_change)
            and rebuilt_at is not None
            and rebuilt_at >= last_change
        )
    details["ledger_row_ok"] = ledger_row_ok

    return snapshots_unchanged, fingerprint_regenerated, ledger_row_ok, details


def _load_probe_offline():
    """Import templates/probe-offline.py by path, relative to this file's location."""
    path = Path(__file__).resolve().parents[1] / "templates" / "probe-offline.py"
    spec = importlib.util.spec_from_file_location("probe_offline", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_reference(args) -> None:
    from hvantk.tests._snapshot_utils import hail_schema_to_dict

    repo = Path(args.repo)
    ht = read_table(Path(args.output))
    schema = hail_schema_to_dict(ht)

    committed_schema_path = OWN_REPO / SNAPSHOTS / "schema.json"
    committed_schema = json.loads(committed_schema_path.read_text())
    contract_ok, contract_details = schema_contract(schema, committed_schema)
    if not contract_ok:
        raise SystemExit(
            "the live build's schema does not satisfy the committed contract at "
            f"{committed_schema_path}: {contract_details}; investigate before "
            "writing a reference from it"
        )

    # The CONTRACT, not the live schema: the committed snapshot's key fields plus its
    # own row field names, sorted -- what the digest below is restricted to, so a
    # later full build's extra columns can never cost it a digest match.
    contract_fields = sorted(
        set(committed_schema.get("key", []))
        | set(committed_schema.get("row", {}).keys())
    )
    keys = select_sample_keys(ht)
    rows_sha = rows_digest(ht, keys, contract_fields)

    probe_offline = _load_probe_offline()
    probe = probe_offline._dataset_probe(DATASET)
    with probe_offline._answer_from_file(args.staged_file):
        observed = probe()

    ledger = _read_json_object(repo / LEDGER) or {}
    committed_row = ledger.get(DATASET)

    ref = {
        "schema_id": SCHEMA_ID,
        "schema": schema,
        "row_count": ht.count(),
        "sample_keys": keys,
        "contract_fields": contract_fields,
        "rows_sha256": rows_sha,
        "snapshots": {
            "schema_sha256": sha256_file(committed_schema_path),
            "sample_rows_sha256": sha256_file(
                OWN_REPO / SNAPSHOTS / "sample_rows.json"
            ),
        },
        "expected_fingerprint": {
            "file": FILENAME,
            "headers": observed["headers"][FILENAME],
            "checksum": observed["checksums"][FILENAME],
        },
        "committed_ledger_row": committed_row,
        **BUDGET,
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "input_sha256": args.input_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(ref, indent=2) + "\n")
    print(f"wrote {args.reference}: row_count={ref['row_count']}, keys={len(keys)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--repo", required=True, help="checkout to read bookkeeping from")
    ap.add_argument("--reference", required=True)
    ap.add_argument(
        "--verdict-dir", default=None, help="default: --output's parent directory"
    )
    ap.add_argument("--run-notes")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--expect-commit", default="")
    ap.add_argument("--expect-sif-sha256", default="")
    ap.add_argument("--input-sha256", default="")
    ap.add_argument("--staged-file", help="the mutated raw file; required to write")
    ap.add_argument("--site-profile")
    ap.add_argument("--sacct-file")
    ap.add_argument("--since", default="now-2days")
    ap.add_argument("--no-slurm", action="store_true")
    ap.add_argument("--expect-fail", action="store_true")
    ap.add_argument("--write-reference", action="store_true")
    args = ap.parse_args()

    if args.write_reference:
        if not args.staged_file:
            ap.error("--staged-file is required with --write-reference")
        write_reference(args)
        return 0
    for required in ("run_id", "expect_commit", "expect_sif_sha256"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required when grading")

    output = Path(args.output)
    repo = Path(args.repo)
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

    # repair_effective FIRST, in its own subprocess, before this process touches Hail.
    committed_schema = json.loads((OWN_REPO / SNAPSHOTS / "schema.json").read_text())
    repair_effective, repair_details = grade_repair_effective(
        repo, committed_schema.get("key")
    )

    try:
        m4, content_ok, content_details = grade_content(
            output, reference, committed_schema
        )
    except Exception as exc:  # noqa: BLE001 -- a corrupt or wrong-type output must
        # not crash the whole verdict; record it and fail M4/content_ok for it alone.
        m4, content_ok = False, False
        content_details = {"error": f"{type(exc).__name__}: {exc}"}
    snapshots_unchanged, fingerprint_regenerated, ledger_row_ok, bk_details = (
        grade_bookkeeping(repo, reference)
    )
    m5 = bool(
        content_ok
        and snapshots_unchanged
        and fingerprint_regenerated
        and ledger_row_ok
        and repair_effective
    )
    passed = bool(m4 and m5)

    verdict = {
        "task": "drift-triage",
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "notes": notes,
        "slurm": slurm_details,
        "output": {
            **content_details,
            "content_ok": content_ok,
            "repair_effective": repair_effective,
            "repair_effective_details": repair_details,
            **bk_details,
        },
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    verdict_dir = (
        Path(args.verdict_dir)
        if args.verdict_dir
        else (output.parent if output.parent.exists() else Path("."))
    )
    verdict_dir.mkdir(parents=True, exist_ok=True)
    (verdict_dir / "verdict.json").write_text(
        json.dumps(verdict, indent=2, default=str)
    )
    sidecar = output.with_name(output.name + ".provenance.json")
    receipt = {
        "output": str(output),
        "repo": str(repo),
        "pass": passed,
        "provenance_sha256": sha256_file(sidecar),
        "reference_sha256": sha256_file(Path(args.reference)),
        "fingerprint_sha256": sha256_file(repo / FINGERPRINT),
        "ledger_sha256": sha256_file(repo / LEDGER),
        "repair_effective": repair_effective,
        "checked_at": verdict["checked_at"],
    }
    (verdict_dir / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(json.dumps({k: verdict[k] for k in ("M1", "M2", "M3", "M4", "M5", "pass")}))
    if args.expect_fail:
        return 0 if not passed else 1
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
