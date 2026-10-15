"""Run one plugin's own test suite from the checkout root.

A plugin's round-trip test resolves its fixture path (for example
``hvantk/skills/clingen/tests/testdata/raw/clingen/...``) relative to the current
working directory, not to the test file -- see `hvantk/skills/_conventions/SKILL.md`
section 8 ("fixture paths are relative to the checkout root") and the plugin's own test
module. Running pytest from anywhere else makes every fixture read fail with a path
that looks right but is not, so this script chdirs into the checkout before invoking
pytest rather than relying on the caller's own working directory.

    python run-plugin-tests.py --repo /path/to/checkout --provider clingen [-- -k gene_disease]
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--repo", required=True, help="hvantk checkout root; becomes the cwd"
    )
    ap.add_argument(
        "--provider", required=True, help="plugin directory name under hvantk/skills/"
    )
    ap.add_argument(
        "pytest_args",
        nargs=argparse.REMAINDER,
        metavar="-- [pytest args]",
        help="extra pytest arguments, after a literal --",
    )
    args = ap.parse_args()

    extra = args.pytest_args
    if extra and extra[0] == "--":
        extra = extra[1:]

    os.chdir(args.repo)
    import pytest

    return pytest.main(
        [
            f"hvantk/skills/{args.provider}/tests",
            "-m",
            "hail",
            "-p",
            "no:cacheprovider",
            *extra,
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
