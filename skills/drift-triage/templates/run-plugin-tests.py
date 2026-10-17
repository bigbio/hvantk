"""Run one dataset's own test command from the checkout root.

A plugin's round-trip test resolves its fixture path relative to the current working
directory, not to the test file -- see `hvantk/skills/_conventions/SKILL.md` section 8
and the plugin's own test module. Running pytest from anywhere else makes every
fixture read fail with a path that looks right but is not, so this script chdirs into
the checkout before running anything, rather than relying on the caller's own working
directory.

The command itself is not hardcoded: it is read from the dataset's own manifest
(`plugin.yaml`'s `tests.command`, surfaced as `DatasetSpec.test_paths.command`), so a
provider whose test command differs from the usual `pytest .../tests -m hail` --
multi-dataset providers sometimes narrow it with `-k`, for example -- still runs
correctly.

    python run-plugin-tests.py --repo /path/to/checkout --dataset clingen:gene-disease [-- -k gene_disease]
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--repo", required=True, help="hvantk checkout root; becomes the cwd"
    )
    ap.add_argument(
        "--dataset", required=True, metavar="PROVIDER:DATASET", dest="dataset"
    )
    ap.add_argument(
        "extra_args",
        nargs=argparse.REMAINDER,
        metavar="-- [extra args]",
        help="extra arguments appended to the resolved command, after a literal --",
    )
    args = ap.parse_args()

    extra = args.extra_args
    if extra and extra[0] == "--":
        extra = extra[1:]

    os.chdir(args.repo)

    from hvantk.core.plugin.loader import get_registry

    command = get_registry().get_dataset(args.dataset).test_paths.command
    cmd = shlex.split(command)
    if not cmd:
        print(f"error: empty tests.command for {args.dataset!r}", file=sys.stderr)
        return 1

    if cmd[0] == "pytest":
        import pytest

        return pytest.main([*cmd[1:], "-p", "no:cacheprovider", *extra])
    return subprocess.run([*cmd, *extra]).returncode


if __name__ == "__main__":
    sys.exit(main())
