"""The shipped examples/rerank/ example still runs and still matches its committed output.

The Availability statement points reviewers at this example, so it has to keep working:
without a test, a change to the engine, the Config surface, or the CLI would break it
silently and nobody would notice until a reviewer tried it.

Unmarked on purpose, so the default ``pytest`` selection collects it -- pytest.ini's
testpaths are hvantk/tests and hvantk/skills, so nothing under examples/ is ever collected
on its own. Needs no Hail; the rerank engine is pandas + scikit-learn only.

Scores are compared with a tolerance rather than byte-for-byte. They are byte-identical on
a fixed dependency set (and verified stable across OMP thread counts), but scikit-learn
patch releases can move the last decimal places, and asserting equality would turn a
harmless upgrade into a spurious failure. Gene set, column set and ranking are asserted
exactly, because those must not move.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("sklearn", reason="rerank needs the 'ml' extra (scikit-learn)")

REPO_ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_DIR = REPO_ROOT / "examples" / "rerank"
CONFIG = EXAMPLE_DIR / "config.yaml"
REFERENCE = EXAMPLE_DIR / "results" / "rerank_scores.tsv"

EXPECTED_COLUMNS = [
    "gene",
    "prior_stat",
    "score",
    "score_percentile",
    "tier",
    "verdict",
    "flag",
    "flag_reason",
    "y",
]


def test_example_files_are_present():
    """Guards the Availability statement: the example a reviewer is pointed at exists."""
    missing = [
        p.relative_to(REPO_ROOT)
        for p in (
            CONFIG,
            REFERENCE,
            EXAMPLE_DIR / "README.md",
            EXAMPLE_DIR / "make_fixtures.py",
            EXAMPLE_DIR / "data" / "cohort.yaml",
            EXAMPLE_DIR / "data" / "constraint.tsv",
            EXAMPLE_DIR / "data" / "expression.tsv",
            EXAMPLE_DIR / "data" / "labels.txt",
            EXAMPLE_DIR / "data" / "cohort.tsv",
        )
        if not p.exists()
    ]
    assert not missing, f"examples/rerank is incomplete: {missing}"


def test_reference_output_has_the_documented_shape():
    df = pd.read_csv(REFERENCE, sep="\t")
    assert list(df.columns) == EXPECTED_COLUMNS
    assert len(df) == 150
    assert df["gene"].is_unique
    # The audit is advisory, so flagged genes must still carry a score and a tier.
    flagged = df[df["flag"]]
    assert not flagged.empty, "fixtures should exercise the audit"
    assert flagged["tier"].notna().all()
    assert set(df.loc[df["flag"], "flag_reason"]) <= {
        "insufficient_data",
        "recurrent_variant",
        "common_driver",
    }


@pytest.mark.slow
def test_example_run_reproduces_the_committed_output(tmp_path):
    """Run the documented command and compare against results/rerank_scores.tsv."""
    out = tmp_path / "rerank_scores.tsv"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "hvantk.hvantk",
            "rerank",
            "-c",
            str(CONFIG),
            "-o",
            str(out),
        ],
        cwd=REPO_ROOT,  # config.yaml paths resolve against the working directory
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 0, f"rerank failed:\n{proc.stdout}\n{proc.stderr}"

    got = pd.read_csv(out, sep="\t")
    ref = pd.read_csv(REFERENCE, sep="\t")

    assert list(got.columns) == list(ref.columns)
    assert got["gene"].tolist() == ref["gene"].tolist()
    # Ranking is the contract; the raw float is allowed to wobble in its last places.
    assert got["score"].rank().tolist() == ref["score"].rank().tolist()
    assert got["score"].to_numpy() == pytest.approx(ref["score"].to_numpy(), abs=1e-6)
    assert got["flag"].tolist() == ref["flag"].tolist()
    assert (
        got["flag_reason"].fillna("").tolist() == ref["flag_reason"].fillna("").tolist()
    )
