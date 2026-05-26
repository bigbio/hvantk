"""Regression test for D1 — drift_runner._IGNORED_KEYS is inconsistent with
the run_builder fingerprint-hash ignore set.

`run_builder._coerce_fingerprint` excludes both `fetched_at` and `probe_version`
when canonicalizing a fingerprint hash for an artifact's provenance. The drift
comparison in `drift_runner._compare_fingerprints` excludes only `fetched_at`.

Consequence: bumping a probe's PROBE_VERSION (a routine code change) will mark
every dataset's drift status as `drifted` even though source data is unchanged,
while the stored artifact fingerprint won't reflect the bump. The drift check
fires false positives that operators have to manually distinguish from real
upstream drift.

Fix: unify the two ignore lists — `_IGNORED_KEYS = frozenset({"fetched_at",
"probe_version"})` in both modules, ideally by hoisting the constant into
`hvantk.core.plugin.api` and importing it from both consumers.
"""

from __future__ import annotations

from hvantk.core.plugin.drift_runner import _compare_fingerprints


def test_probe_version_bump_does_not_count_as_drift():
    """A probe-implementation version bump should not flip status to 'drifted'.

    Fails on 9dbb654: probe_version differs → _compare_fingerprints returns
    a non-None diff → drift status would be 'drifted'.
    Passes after fix: probe_version is excluded → returns None.
    """
    expected = {
        "probe_version": 1,
        "source_version": "Sun, 04 May 2026 12:08:00 GMT",
        "headers": {"clinvar.vcf.gz": {"content_length": "97582033"}},
        "fetched_at": "2026-05-16T00:00:00+00:00",
    }
    observed = {
        "probe_version": 2,
        "source_version": "Sun, 04 May 2026 12:08:00 GMT",
        "headers": {"clinvar.vcf.gz": {"content_length": "97582033"}},
        "fetched_at": "2026-05-22T00:00:00+00:00",
    }

    diff = _compare_fingerprints(expected, observed)
    assert diff is None, (
        f"probe_version bump should be ignored in drift comparison "
        f"(consistent with run_builder._coerce_fingerprint), got diff={diff!r}"
    )
