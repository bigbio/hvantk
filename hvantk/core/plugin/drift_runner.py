"""Compare a plugin's live drift-probe fingerprint against its committed expected
fingerprint. Returns a structured DriftResult that the CLI / CI can consume.
"""

from __future__ import annotations

import json
import signal
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable

from .api import (
    DatasetSpec,
    DriftProbeError,
    PROBE_FINGERPRINT_IGNORED_KEYS,
    PROBE_STATUS_STUB,
    placeholder_baseline_reason,
)


@dataclass(frozen=True)
class DriftResult:
    dataset_name: str
    status: str  # clean | drifted | probe_failed | stub
    observed: dict[str, Any] | None = None
    expected: dict[str, Any] | None = None
    diff: dict[str, Any] | None = None
    probe_error: BaseException | None = None
    #: The committed baseline this result was diffed against. Datasets that share one
    #: share a drift signal, which is what lets the CI bot collapse them into a single
    #: PR instead of one per dataset. See `run_drift_checks`.
    fingerprint_path: str | None = None


def run_drift_check(dataset_name: str, *, timeout: int = 60) -> DriftResult:
    """Resolve dataset from the module registry, invoke probe, diff."""
    from . import loader as plugin_loader

    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(dataset_name)
    return _run_drift_check_with_spec(spec, timeout=timeout)


def run_drift_checks(
    specs: "Iterable[DatasetSpec]", *, timeout: int = 60
) -> list[DriftResult]:
    """Drift-check many datasets, probing once per distinct drift signal.

    Datasets that declare the SAME ``drift_fingerprint`` baseline and resolve to the
    SAME probe callable do not have separate drift signals -- they have one, reported
    several times. `ucsc-cellbrowser` is the worked example: `default`, `adult-ctx` and
    `dev-ctx` are distinct *schema* variants (their obs cell-type column is `celltype`,
    `Class` and `Type_v2` respectively, which is why each earns its own snapshot), but
    `fetch_fingerprint()` takes no arguments and fingerprints the provider-wide catalog
    at cells.ucsc.edu/dataset.json. One upstream event therefore produced three
    identical drift reports, three branches writing the same file, and three mutually
    conflicting PRs -- merging any one made the other two conflict.

    Grouping is keyed on ``(fingerprint path, probe callable)`` rather than the path
    alone. Two datasets sharing a baseline but resolving to DIFFERENT probes is a
    manifest error (the two would overwrite each other's baseline), and this key makes
    the runtime conservative about it: they simply do not group, and each is probed and
    reported on its own. `hvantk plugins validate` reports the misconfiguration.

    Every dataset still gets its own entry in the returned list, so the report shape is
    unchanged; members of a group share ``fingerprint_path``, which is what lets the CI
    bot open one PR per signal instead of one per dataset.
    """
    cache: dict[tuple[str, int], DriftResult] = {}
    results: list[DriftResult] = []
    for spec in specs:
        key = (str(spec.test_paths.drift_fingerprint), id(spec.drift_probe))
        cached = cache.get(key)
        if cached is None:
            cached = _run_drift_check_with_spec(spec, timeout=timeout)
            cache[key] = cached
        # `replace` rather than reuse: the shared result carries the FIRST member's
        # dataset_name, and every caller keys off that field.
        results.append(replace(cached, dataset_name=spec.name))
    return results


def _probe_failed(spec: DatasetSpec, exc: DriftProbeError) -> DriftResult:
    """Build a probe_failed result, also flagging a missing baseline fingerprint.

    The probe runs before the baseline is read (so intentional stubs classify
    correctly). Without this, a plugin whose probe fails AND ships no committed
    baseline would surface only the probe error and silently hide the
    missing-fingerprint configuration issue. The underlying probe error is
    preserved either way.
    """
    fp_path = Path(spec.test_paths.drift_fingerprint)
    if not fp_path.exists():
        exc = DriftProbeError(
            f"{exc}; additionally, expected fingerprint is missing at {fp_path}"
        )
    return DriftResult(
        dataset_name=spec.name,
        status="probe_failed",
        probe_error=exc,
        fingerprint_path=str(fp_path),
    )


def _run_drift_check_with_spec(
    spec: DatasetSpec, *, timeout: int = 60
) -> DriftResult:
    # Resolved up front because every return below reports it, including the stub
    # branch, which returns before the baseline is read. Reading the PATH is not
    # reading the FILE, so this does not disturb the probe-before-baseline ordering
    # described below.
    fp_path = Path(spec.test_paths.drift_fingerprint)

    # Invoke the probe before loading the baseline so an intentional stub
    # (documentation-only source with no probeable URL) is reported as
    # status="stub" — these plugins ship no committed baseline, so a
    # baseline-first ordering would mislabel them as "probe_failed".
    try:
        observed = _invoke_with_timeout(spec.drift_probe, timeout=timeout)
    except DriftProbeError as exc:
        return _probe_failed(spec, exc)
    except Exception as exc:  # noqa: BLE001
        return _probe_failed(
            spec, DriftProbeError(f"probe raised {type(exc).__name__}: {exc}")
        )

    if observed.get("probe_status") == PROBE_STATUS_STUB:
        return DriftResult(
            dataset_name=spec.name,
            status="stub",
            observed=observed,
            fingerprint_path=str(fp_path),
        )

    try:
        expected = json.loads(fp_path.read_text())
    except FileNotFoundError:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            observed=observed,
            probe_error=DriftProbeError(
                f"missing expected fingerprint at {fp_path}"
            ),
            fingerprint_path=str(fp_path),
        )

    # A hand-seeded baseline cannot equal a live observation, so diffing it would
    # report drift forever. That is a missing baseline wearing a committed file's
    # clothes, and it classifies as probe_failed rather than drifted.
    seeded = placeholder_baseline_reason(expected)
    if seeded is not None:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            observed=observed,
            expected=expected,
            probe_error=DriftProbeError(
                f"committed baseline at {fp_path} was never captured from a live "
                f"probe ({seeded}); run `hvantk drift --regenerate {spec.name}`"
            ),
            fingerprint_path=str(fp_path),
        )

    diff = _compare_fingerprints(expected, observed)
    if diff is None:
        return DriftResult(
            dataset_name=spec.name,
            status="clean",
            observed=observed,
            expected=expected,
            fingerprint_path=str(fp_path),
        )
    return DriftResult(
        dataset_name=spec.name,
        status="drifted",
        observed=observed,
        expected=expected,
        diff=diff,
        fingerprint_path=str(fp_path),
    )


def _invoke_with_timeout(fn, *, timeout: int) -> dict:
    """Run fn() under a signal-based timeout (POSIX). Falls back to no
    timeout on platforms where SIGALRM is unavailable.

    Caveat: on POSIX this overwrites any pre-existing SIGALRM handler and
    cancels any pending alarm. Safe to use in CLI entrypoints; nesting calls
    or running alongside other alarm-using code is not supported.
    """
    if not hasattr(signal, "SIGALRM"):
        return _coerce_fingerprint(fn())

    def _handler(signum, frame):
        raise DriftProbeError(f"probe timed out after {timeout}s")

    prev = signal.signal(signal.SIGALRM, _handler)
    try:
        signal.alarm(timeout)
        result = fn()
        return _coerce_fingerprint(result)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


def _coerce_fingerprint(result) -> dict:
    """Convert a probe's return value to a plain dict, raising a clear
    DriftProbeError if the value isn't dict-like."""
    if isinstance(result, dict):
        return result
    try:
        return dict(result)
    except (TypeError, ValueError) as exc:
        raise DriftProbeError(
            f"probe returned non-mapping value: {type(result).__name__} "
            f"({result!r})"
        ) from exc


def _compare_fingerprints(
    expected: dict[str, Any], observed: dict[str, Any]
) -> dict[str, Any] | None:
    """Return None if equal (ignoring PROBE_FINGERPRINT_IGNORED_KEYS), else a structured diff."""
    def strip(d: dict) -> dict:
        return {k: v for k, v in d.items() if k not in PROBE_FINGERPRINT_IGNORED_KEYS}

    e = strip(expected)
    o = strip(observed)
    if e == o:
        return None
    added = {k: o[k] for k in o.keys() - e.keys()}
    removed = {k: e[k] for k in e.keys() - o.keys()}
    changed = {
        k: {"expected": e[k], "observed": o[k]}
        for k in e.keys() & o.keys()
        if e[k] != o[k]
    }
    return {"added": added, "removed": removed, "changed": changed}
