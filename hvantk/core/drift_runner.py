"""Compare a plugin's live drift-probe fingerprint against its committed expected
fingerprint. Returns a structured DriftResult that the CLI / CI can consume.
"""

from __future__ import annotations

import json
import signal
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .plugin_api import DatasetSpec, DriftProbeError

# Keys excluded from the fingerprint comparison.
_IGNORED_KEYS = frozenset({"fetched_at"})


@dataclass
class DriftResult:
    dataset_name: str
    status: str  # clean | drifted | probe_failed
    observed: dict[str, Any] | None = None
    expected: dict[str, Any] | None = None
    diff: dict[str, Any] | None = None
    probe_error: BaseException | None = None


def run_drift_check(dataset_name: str, *, timeout: int = 60) -> DriftResult:
    """Resolve dataset from the module registry, invoke probe, diff."""
    from . import plugin_loader

    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(dataset_name)
    return _run_drift_check_with_spec(spec, timeout=timeout)


def _run_drift_check_with_spec(
    spec: DatasetSpec, *, timeout: int = 60
) -> DriftResult:
    fp_path = Path(spec.test_paths.drift_fingerprint)
    try:
        expected = json.loads(fp_path.read_text())
    except FileNotFoundError:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            probe_error=DriftProbeError(
                f"missing expected fingerprint at {fp_path}"
            ),
        )

    try:
        observed = _invoke_with_timeout(spec.drift_probe, timeout=timeout)
    except DriftProbeError as exc:
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            expected=expected,
            probe_error=exc,
        )
    except Exception as exc:  # noqa: BLE001
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            expected=expected,
            probe_error=DriftProbeError(f"probe raised {type(exc).__name__}: {exc}"),
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
    timeout on platforms where SIGALRM is unavailable."""
    if not hasattr(signal, "SIGALRM"):
        return dict(fn())

    def _handler(signum, frame):
        raise DriftProbeError(f"probe timed out after {timeout}s")

    prev = signal.signal(signal.SIGALRM, _handler)
    try:
        signal.alarm(timeout)
        result = fn()
        return dict(result)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)


def _compare_fingerprints(
    expected: dict[str, Any], observed: dict[str, Any]
) -> dict[str, Any] | None:
    """Return None if equal (ignoring _IGNORED_KEYS), else a structured diff."""
    def strip(d: dict) -> dict:
        return {k: v for k, v in d.items() if k not in _IGNORED_KEYS}

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
