"""Assert the built wheel ships what it should and nothing else.

Run by the `packaging-smoke` CI job after `poetry build -f wheel`. Nothing else in the
pipeline builds a wheel, so before this existed the include/exclude globs in pyproject.toml
were never checked and had silently rotted: the wheel carried 339 files under a tests/
directory and 43.3 MB of fixtures (a 14.7 MB VDS zip, an 11.4 MB expression-atlas file) in a
46.2 MB package, while the excludes that were supposed to stop it either named the wrong
directory (`tests/data` where the skills use `tests/testdata`) or were overridden by a
broader `include`.

The size ceiling is deliberately generous. It is a tripwire for "someone added a fixture to a
shipped path", not a byte budget -- a real regression here is an order of magnitude, not a few
percent.
"""
from __future__ import annotations

import glob
import re
import sys
import zipfile
from pathlib import Path

MAX_UNCOMPRESSED_MB = 10.0
REPO = Path(__file__).resolve().parents[2]


def main() -> int:
    wheels = glob.glob(str(REPO / "dist" / "*.whl"))
    if len(wheels) != 1:
        print(f"FAIL expected exactly one wheel in dist/, found {len(wheels)}: {wheels}")
        return 1
    names = zipfile.ZipFile(wheels[0]).infolist()
    paths = [i.filename for i in names]
    failures: list[str] = []

    def count(pred) -> int:
        return sum(1 for p in paths if pred(p))

    # --- things that must NOT ship
    top_level_tests = count(lambda p: p.startswith("hvantk/tests/"))
    if top_level_tests:
        failures.append(f"{top_level_tests} files under hvantk/tests/ (expected 0)")

    test_modules = count(lambda p: re.search(r"/test_[^/]*\.py$", p))
    if test_modules:
        failures.append(f"{test_modules} test_*.py modules (expected 0)")

    # Drift fingerprints are the deliberate exception: `hvantk drift` compares a live probe
    # against them at runtime, so they belong in the wheel. Multi-dataset providers suffix the
    # name per dataset (drift_fingerprint_samples.json), hence the prefix match rather than an
    # exact one.
    def _is_fingerprint(p: str) -> bool:
        return re.search(r"/drift_fingerprint[^/]*\.json$", p) is not None

    fixtures = count(
        lambda p: "/tests/" in p and not (_is_fingerprint(p) or p.endswith(".gitkeep"))
    )
    if fixtures:
        failures.append(f"{fixtures} test fixtures under a tests/ dir (only fingerprints ship)")

    mb = sum(i.file_size for i in names) / 1e6
    if mb > MAX_UNCOMPRESSED_MB:
        failures.append(f"wheel is {mb:.1f} MB uncompressed, ceiling is {MAX_UNCOMPRESSED_MB}")

    # --- things that MUST ship. Counted against the tree, so adding a provider needs no edit
    # here; a packaging glob that silently drops one is what this catches.
    # Manifests alone are not enough. `hvantk plugins list` -- the CI job's other check --
    # only reads YAML (loader Pass 1, no callables imported), so a packaging glob that dropped
    # every builder.py would still print the right provider count and leave the job green.
    # The skills-module count is what actually covers that, and the fingerprints are the
    # runtime data `hvantk drift` compares against.
    # `keep_disk` filters the TREE side so both sides count the same population. It matters for
    # skills modules (test_*.py are excluded from the wheel on purpose, so counting them on
    # disk would fail every run) and must NOT be applied to fingerprints, which legitimately
    # live under tests/.
    not_a_test = lambda p: "/tests/" not in Path(p).as_posix()  # noqa: E731
    everything = lambda p: True                                  # noqa: E731

    for label, pattern, wheel_pred, keep_disk in (
        ("plugin.yaml manifests", "hvantk/skills/**/plugin.yaml",
         lambda p: p.endswith("plugin.yaml"), everything),
        ("SKILL.md files", "hvantk/skills/**/SKILL.md",
         lambda p: p.endswith("SKILL.md"), everything),
        ("catalog datasets.json", "hvantk/skills/**/catalog/*.json",
         lambda p: "/catalog/" in p and p.endswith(".json"), everything),
        ("drift fingerprints", "hvantk/skills/**/tests/drift_fingerprint*.json",
         _is_fingerprint, everything),
        ("skills modules", "hvantk/skills/**/*.py",
         lambda p: p.startswith("hvantk/skills/") and p.endswith(".py"), not_a_test),
    ):
        on_disk = len([p for p in glob.glob(str(REPO / pattern), recursive=True)
                       if keep_disk(p)])
        in_wheel = count(wheel_pred)
        if on_disk != in_wheel:
            failures.append(f"{label}: {on_disk} in tree but {in_wheel} in wheel")

    if not count(lambda p: p.endswith("entry_points.txt")):
        failures.append("no entry_points.txt -- console script and provider plugins are unregistered")

    for f in failures:
        print(f"FAIL {f}")
    if failures:
        return 1
    print(f"OK wheel is {mb:.2f} MB uncompressed, {len(paths)} files, no test payload")
    return 0


if __name__ == "__main__":
    sys.exit(main())
