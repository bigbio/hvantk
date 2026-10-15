"""Deterministic check for the hpc-sync eval task.

Grades one run against the staged scenario: M1 from ``run-notes.md``; M2 and M3 are
null because the task needs no job; M4 is "the cluster remote's main and dev equal the
GitHub remote's"; M5 is "the cluster-only branch is untouched and the GitHub remote was
not changed", both against ``reference.json``. No git, Hail or Slurm: refs are read from the
bare repositories' files, so it runs inside the image on a login node.

    python check.py --scenario DIR --run-notes NOTES.md --reference reference.json \
        --expect-commit SHA --expect-sif-sha256 HEX [--expect-fail]
    python check.py --scenario DIR --reference reference.json --write-reference \
        --expect-commit SHA --expect-sif-sha256 HEX
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

BRANCHES = ("main", "dev")
CLUSTER_ONLY = ("cluster-work",)

_COMMIT = re.compile(r"^commit:\s*([0-9a-f]{7,40})\s*$", re.M)
_CONTAINER = re.compile(r"^container:\s*([0-9a-f]{64})\s*$", re.M)


def parse_notes(path: Path) -> dict:
    if not path.is_file():
        return {"present": False, "commit": None, "container": None}
    text = path.read_text()
    c = _COMMIT.search(text)
    k = _CONTAINER.search(text)
    return {
        "present": True,
        "commit": c.group(1) if c else None,
        "container": k.group(1) if k else None,
    }


def ref(repo: Path, name: str) -> str | None:
    """The commit a branch of a bare repository points at, read from the repository's
    own files: no git binary is needed, and the image has none."""
    loose = repo / "refs" / "heads" / name
    if loose.is_file():
        return loose.read_text().strip() or None
    packed = repo / "packed-refs"
    if packed.is_file():
        for line in packed.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1] == f"refs/heads/{name}":
                return parts[0]
    return None


def refs(repo: Path, names: tuple[str, ...]) -> dict[str, str | None]:
    return {n: ref(repo, n) for n in names}


def write_reference(args) -> None:
    scenario = Path(args.scenario)
    reference = {
        "github": refs(scenario / "github.git", BRANCHES),
        "cluster_only": refs(scenario / "cluster.git", CLUSTER_ONLY),
        "produced_by": {
            "commit": args.expect_commit,
            "sif_sha256": args.expect_sif_sha256,
            "at": dt.datetime.now(dt.timezone.utc).isoformat(),
        },
    }
    Path(args.reference).write_text(json.dumps(reference, indent=2) + "\n")
    print(f"wrote {args.reference}: {reference['github']}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--run-notes")
    ap.add_argument("--expect-commit", default="")
    ap.add_argument("--expect-sif-sha256", default="")
    ap.add_argument("--expect-fail", action="store_true")
    ap.add_argument("--write-reference", action="store_true")
    args = ap.parse_args()

    if args.write_reference:
        write_reference(args)
        return 0
    for required in ("expect_commit", "expect_sif_sha256"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required when grading")

    scenario = Path(args.scenario)
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

    github = refs(scenario / "github.git", BRANCHES)
    cluster = refs(scenario / "cluster.git", BRANCHES)
    cluster_only = refs(scenario / "cluster.git", CLUSTER_ONLY)
    m4 = all(cluster[b] is not None and cluster[b] == github[b] for b in BRANCHES)
    github_untouched = github == reference["github"]
    cluster_only_kept = cluster_only == reference["cluster_only"]
    m5 = github_untouched and cluster_only_kept
    passed = bool(m4 and m5)

    verdict = {
        "task": "hpc-sync",
        "M1": m1,
        "M2": None,
        "M3": None,
        "M4": m4,
        "M5": m5,
        "pass": passed,
        "notes": notes,
        "slurm": {"skipped": "task needs no job"},
        "output": {
            "github": github,
            "cluster": cluster,
            "cluster_only": cluster_only,
            "github_untouched": github_untouched,
            "cluster_only_kept": cluster_only_kept,
        },
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    out_dir = scenario.parent
    (out_dir / "verdict.json").write_text(json.dumps(verdict, indent=2))
    state = json.dumps(
        {"github": github, "cluster": cluster, "cluster_only": cluster_only},
        sort_keys=True,
    )
    receipt = {
        "scenario": str(scenario),
        "pass": passed,
        "refs_sha256": hashlib.sha256(state.encode()).hexdigest(),
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
