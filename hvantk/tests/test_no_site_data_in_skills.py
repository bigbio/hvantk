"""Reject site data in anything under ``skills/`` and ``.claude-plugin/``.

Workflow skills run on a shared cluster, so the material they are distilled from is
full of usernames, home paths, hostnames, partition names and e-mail addresses. None of
that belongs in a public repository. This test is the generic, always-on line: it knows
nothing about any particular site and rejects the *shapes* such data takes. The
site-specific list lives in the author's local publish guard and never enters the repo.

Deliberately generic patterns only. Do not add a real username, host or partition name
here to "catch it": that would publish the thing being protected.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCANNED_DIRS = ("skills", ".claude-plugin")
TEXT_SUFFIXES = {
    ".md",
    ".sh",
    ".sbatch",
    ".py",
    ".json",
    ".yaml",
    ".yml",
    ".txt",
    ".toml",
}

#: Public hosts a skill may legitimately name. Add to this list on purpose, never a
#: cluster, institute or personal host.
ALLOWED_HOSTS = {
    "github.com",
    "raw.githubusercontent.com",
    "storage.googleapis.com",
    "www.genenames.org",
    "genenames.org",
    "hail.is",
    "docs.anthropic.com",
    "anthropic.com",
    "claude.com",
    "agentskills.io",
    "example.org",
    "example.com",
    "ensembl.org",
    "ftp.ensembl.org",
    "ebi.ac.uk",
    "ftp.ebi.ac.uk",
    "proot.gitlab.io",
    "apptainer.org",
    "sylabs.io",
    "slurm.schedmd.com",
    "schedmd.com",
    "python.org",
    "readthedocs.io",
}

_HOME_PATH = re.compile(r"(?<![\w/])/(?:Users|home|user)/[A-Za-z0-9._-]+")
_EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@([A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+)")
_SSH_USER_AT = re.compile(r"\bssh\s+(?:-\S+\s+)*[A-Za-z0-9._-]+@")
_SBATCH_SITE_DIRECTIVE = re.compile(
    r"#SBATCH\s+--(?:partition|account|mail-user|qos)\b"
)
_LITERAL_PARTITION = re.compile(r"--(?:partition|account|qos)=(?!['\"]?[$<{])\S+")
_MAIL_USER = re.compile(r"--mail-user[= ](?!['\"]?[$<])\S+")
_HOST = re.compile(
    r"\b(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+"
    r"(?:ac\.uk|de|edu|org|com|net|eu|uk|io|fr|ch|nl|it|es)\b"
)


def findings_in(text: str) -> list[str]:
    """Every site-data shape found in ``text``, as ``<rule>: <match>`` strings."""
    out: list[str] = []
    for m in _HOME_PATH.finditer(text):
        out.append(f"home-path: {m.group(0)}")
    for m in _EMAIL.finditer(text):
        if m.group(1).lower() not in ALLOWED_HOSTS:
            out.append(f"email: {m.group(0)}")
    for m in _SSH_USER_AT.finditer(text):
        out.append(f"ssh-user-at-host: {m.group(0)}")
    for m in _SBATCH_SITE_DIRECTIVE.finditer(text):
        out.append(f"sbatch-site-directive: {m.group(0)}")
    for m in _LITERAL_PARTITION.finditer(text):
        out.append(f"literal-site-flag: {m.group(0)}")
    for m in _MAIL_USER.finditer(text):
        out.append(f"mail-user: {m.group(0)}")
    for m in _HOST.finditer(text):
        host = m.group(0).lower()
        if host not in ALLOWED_HOSTS and not any(
            host.endswith("." + a) for a in ALLOWED_HOSTS
        ):
            out.append(f"host: {host}")
    return out


def scanned_files() -> list[Path]:
    files: list[Path] = []
    for d in SCANNED_DIRS:
        root = REPO_ROOT / d
        if not root.exists():
            continue
        for p in sorted(root.rglob("*")):
            if not p.is_file() or p.suffix not in TEXT_SUFFIXES:
                continue
            if ".local." in p.name:  # the filled site profile, gitignored
                continue
            files.append(p)
    return files


# --- the detector's own tests: each rule must fire on its shape and stay quiet on the
# --- placeholder form a skill is supposed to use instead.


@pytest.mark.parametrize(
    "text, rule",
    [
        ("see /Users/someone/projects/x", "home-path"),
        ("cd /home/someone/work", "home-path"),
        ("ls /user/someone/git", "home-path"),
        ("mail someone@institute.example.net", "email"),
        ("ssh someone@login1", "ssh-user-at-host"),
        ("#SBATCH --partition=batch", "sbatch-site-directive"),
        ("#SBATCH --mail-user=x", "sbatch-site-directive"),
        ("sbatch --partition=batch job.sbatch", "literal-site-flag"),
        ("sbatch --account=proj123 job.sbatch", "literal-site-flag"),
        ("--mail-user=someone@institute.example.net", "mail-user"),
        ("copy to login.cluster.institute.de", "host"),
    ],
)
def test_each_rule_fires_on_its_shape(text, rule):
    assert any(f.startswith(rule + ":") for f in findings_in(text)), findings_in(text)


@pytest.mark.parametrize(
    "text",
    [
        "write under $HVANTK_WORK/out",
        'sbatch --partition="$HVANTK_CPU_PARTITION" job.sbatch',
        "sbatch --partition=$HVANTK_CPU_PARTITION job.sbatch",
        "sbatch --partition=<your partition> job.sbatch",
        'ssh "$HPC_ALIAS" true',
        "export SLURM_EMAIL_NOTIFICATION=you@example.org",
        "https://github.com/bigbio/hvantk",
        "https://storage.googleapis.com/public-download-files/x",
        "hgnc.ht and hgnc.ht.provenance.json",
        "site-profile.local.sh",
        '--mail-user="$SLURM_EMAIL_NOTIFICATION"',
    ],
)
def test_placeholder_forms_are_clean(text):
    assert findings_in(text) == []


def test_scanner_sees_the_skills_tree():
    """Guard the guard: an empty file list would make the tree test vacuous."""
    assert any("skills" in p.parts for p in scanned_files()), (
        "no files found under skills/"
    )


def test_no_site_data_in_skills_tree():
    problems: list[str] = []
    for p in scanned_files():
        for f in findings_in(p.read_text(encoding="utf-8", errors="replace")):
            problems.append(f"{p.relative_to(REPO_ROOT)}: {f}")
    assert not problems, (
        "site data found under skills/ or .claude-plugin/ (use $HVANTK_* variables, "
        "<placeholders> and the site profile instead):\n  " + "\n  ".join(problems)
    )
