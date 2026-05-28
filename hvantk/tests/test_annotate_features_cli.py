"""Regression: --output_ht must be required (issue #123)."""
import subprocess
import sys


def test_output_ht_is_required():
    """The CLI must fail with a non-zero exit when --output_ht is omitted.

    Prior to #123 the default was the literal string ``"None/data/features"``,
    silently producing an oddly-named directory in the working tree.
    """
    result = subprocess.run(
        [sys.executable, "-m", "hvantk.tools.annotation.annotate_features",
         "--variant_ht", "/tmp/does-not-matter.ht"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0, (
        "annotate_features should fail without --output_ht; "
        f"got returncode={result.returncode}, stderr={result.stderr!r}"
    )
    assert "output_ht" in result.stderr.lower() or "required" in result.stderr.lower(), (
        f"Expected argparse to complain about --output_ht; stderr was: {result.stderr!r}"
    )
