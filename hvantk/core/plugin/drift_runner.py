"""Compare a plugin's live drift-probe fingerprint against its committed expected
fingerprint. Returns a structured DriftResult that the CLI / CI can consume.
"""

from __future__ import annotations

import json
import signal
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .api import (
    DatasetSpec,
    DriftProbeError,
    PROBE_FINGERPRINT_IGNORED_KEYS,
    PROBE_STATUS_STUB,
)


@dataclass(frozen=True)
class DriftResult:
    dataset_name: str
    status: str  # clean | drifted | probe_failed | stub
    observed: dict[str, Any] | None = None
    expected: dict[str, Any] | None = None
    diff: dict[str, Any] | None = None
    probe_error: BaseException | None = None


def run_drift_check(dataset_name: str, *, timeout: int = 60) -> DriftResult:
    """Resolve dataset from the module registry, invoke probe, diff."""
    from . import loader as plugin_loader

    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(dataset_name)
    return _run_drift_check_with_spec(spec, timeout=timeout)


def _run_drift_check_with_spec(
    spec: DatasetSpec, *, timeout: int = 60
) -> DriftResult:
    # Invoke the probe before loading the baseline so an intentional stub
    # (documentation-only source with no probeable URL) is reported as
    # status="stub" — these plugins ship no committed baseline, so a
    # baseline-first ordering would mislabel them as "probe_failed".
    try:
        observed = _invoke_with_timeout(spec.drift_probe, timeout=timeout)
    except DriftProbeError as exc:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            probe_error=exc,
        )
    except Exception as exc:  # noqa: BLE001
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            probe_error=DriftProbeError(f"probe raised {type(exc).__name__}: {exc}"),
        )

    if observed.get("probe_status") == PROBE_STATUS_STUB:
        return DriftResult(
            dataset_name=spec.name,
            status="stub",
            observed=observed,
        )

    fp_path = Path(spec.test_paths.drift_fingerprint)
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
        )

    diff = _compare_fingerprints(expected, observed)
    if diff is None:
        return DriftResult(
            dataset_name=spec.name,
            status="clean",
            observed=observed,
            expected=expected,
        )
    return DriftResult(
        dataset_name=spec.name,
        status="drifted",
        observed=observed,
        expected=expected,
        diff=diff,
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
