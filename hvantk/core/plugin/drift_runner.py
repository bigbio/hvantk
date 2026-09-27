"""Compare a plugin's live drift-probe fingerprint against its committed expected
fingerprint. Returns a structured DriftResult that the CLI / CI can consume.
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Mapping

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
    #: The committed baseline this result was diffed against. Datasets sharing a
    #: baseline share a single drift signal, which is what lets the CI bot collapse them
    #: into one PR instead of one per dataset. See `run_drift_checks`.
    fingerprint_path: str | None = None
    #: ``module:function`` of the probe that produced this result. Pairs with
    #: ``fingerprint_path`` to identify the SIGNAL: the runner refuses to group two
    #: datasets that share a baseline but declare different probes (a manifest error),
    #: and the bot must apply the same rule or it would re-merge what the runner
    #: deliberately kept apart -- opening one PR whose regeneration covers only the first
    #: dataset and silently leaving the second's drift unaddressed.
    probe_ref: str | None = None


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
    # Values keep the probe object alive alongside its result. The key holds only
    # id(probe), and CPython reuses addresses: `specs` is an Iterable, so a caller may
    # pass a generator whose specs become unreachable as it advances, and a freshly
    # allocated probe could land on a freed address. Two distinct probes would then
    # collide on one key and a dataset would receive another dataset's drift result.
    # Holding the reference makes the address un-reusable for the loop's lifetime.
    cache: dict[tuple[str, int], tuple[object, DriftResult]] = {}
    results: list[DriftResult] = []
    for spec in specs:
        probe = spec.drift_probe
        key = (str(spec.test_paths.drift_fingerprint), id(probe))
        entry = cache.get(key)
        if entry is None:
            cached = _run_drift_check_with_spec(spec, timeout=timeout)
            cache[key] = (probe, cached)
        else:
            cached = entry[1]
        # `replace` rather than reuse: the shared result carries the FIRST member's
        # dataset_name, and every caller keys off that field.
        results.append(replace(cached, dataset_name=spec.name))
    return results


def _probe_ref(spec: DatasetSpec) -> str:
    """``module:qualname`` of a spec's probe -- a serialisable stand-in for the callable.

    The runner can compare callables by identity; the JSON report the CI bot consumes
    cannot, so the identity has to survive serialisation for the bot to apply the same
    grouping rule.
    """
    fn = spec.drift_probe
    return f"{getattr(fn, '__module__', '?')}:{getattr(fn, '__qualname__', repr(fn))}"


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
        probe_ref=_probe_ref(spec),
    )


def _run_drift_check_with_spec(spec: DatasetSpec, *, timeout: int = 60) -> DriftResult:
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
            probe_ref=_probe_ref(spec),
        )

    def _bad_baseline(message: str, *, expected: dict[str, Any] | None = None) -> DriftResult:
        """Build this dataset's probe_failed row for a baseline that is present but
        unusable (missing, unreadable, unparseable, wrong shape, or hand-seeded).

        Every call site below differs only in ``message`` (and the seeded-baseline
        case, which also has an ``expected`` to show alongside the observed value);
        collapsing them here keeps `dataset_name`/`observed`/`fingerprint_path`/
        `probe_ref` from drifting out of sync across the four sites.
        """
        return DriftResult(
            dataset_name=spec.name,
            status="probe_failed",
            observed=observed,
            expected=expected,
            probe_error=DriftProbeError(message),
            fingerprint_path=str(fp_path),
            probe_ref=_probe_ref(spec),
        )

    try:
        expected = json.loads(fp_path.read_text())
    except FileNotFoundError:
        return _bad_baseline(f"missing expected fingerprint at {fp_path}")
    except (OSError, ValueError) as exc:
        # A baseline that exists but cannot be read (permissions, a directory) or parsed
        # (truncated, empty -- json.JSONDecodeError is a ValueError). Before #361 this
        # propagated as a traceback, which exits 1 == EXIT_DRIFTED: drift.yml gates on
        # rc > 2, so the step went green, and because the JSON report is echoed only after
        # every check completes, one bad file discarded the report for all 26 datasets.
        # It is this dataset's problem, so it is this dataset's probe_failed row.
        return _bad_baseline(
            f"committed baseline at {fp_path} could not be read as JSON "
            f"({type(exc).__name__}: {exc}); run `hvantk drift --regenerate {spec.name}`"
        )
    if not isinstance(expected, dict):
        # `placeholder_baseline_reason` and `_compare_fingerprints` both call `.get`/`.items`
        # on this, so a JSON list or scalar died there with an AttributeError instead.
        return _bad_baseline(
            f"committed baseline at {fp_path} is not a JSON object (got "
            f"{type(expected).__name__}); run `hvantk drift --regenerate {spec.name}`"
        )

    # A hand-seeded baseline cannot equal a live observation, so diffing it would
    # report drift forever. That is a missing baseline wearing a committed file's
    # clothes, and it classifies as probe_failed rather than drifted.
    seeded = placeholder_baseline_reason(expected)
    if seeded is not None:
        return _bad_baseline(
            f"committed baseline at {fp_path} was never captured from a live "
            f"probe ({seeded}); run `hvantk drift --regenerate {spec.name}`",
            expected=expected,
        )

    diff = _compare_fingerprints(expected, observed)
    if diff is None:
        return DriftResult(
            dataset_name=spec.name,
            status="clean",
            observed=observed,
            expected=expected,
            fingerprint_path=str(fp_path),
            probe_ref=_probe_ref(spec),
        )
    return DriftResult(
        dataset_name=spec.name,
        status="drifted",
        observed=observed,
        expected=expected,
        diff=diff,
        fingerprint_path=str(fp_path),
        probe_ref=_probe_ref(spec),
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


def _current_umask() -> int:
    """Read the process umask without changing it.

    ``os.umask`` is set-and-return-the-old-value -- there is no read-only form -- so
    the only way to read it is to set some value and immediately set it back.
    """
    mask = os.umask(0)
    os.umask(mask)
    return mask


def write_fingerprint(path: Path, fingerprint: Mapping[str, Any]) -> None:
    """Write a fingerprint file atomically, in the byte format `--regenerate` always used.

    Temp file in the SAME directory, then ``os.replace``: rename is atomic only within one
    filesystem, and the drift bot runs ``--regenerate`` unattended, so an interrupted or
    out-of-disk write must leave either the old baseline or the new one -- never a
    truncated file, which is exactly the crash the baseline reader above now reports.
    The serialisation is deliberately identical to what the CLI wrote before (indent=2,
    default=str, no trailing newline) so regenerated files do not churn.

    ``flush`` + ``fsync`` happen before the rename so the new bytes are durable on disk
    before the name is ever repointed at them -- otherwise a crash could let the rename
    outrun the write.

    ``tempfile.mkstemp`` always creates its file mode 0600, which -- unlike the in-place
    ``write_text`` this replaced -- would silently tighten a committed baseline's mode on
    every regenerate. When ``path`` already exists, its current mode is copied onto the
    temp file before the replace; for a brand-new baseline, the mode a normal file create
    would get (0666 masked by the process umask) is applied instead.
    """
    path = Path(path)
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(json.dumps(fingerprint, indent=2, default=str))
            fh.flush()
            os.fsync(fh.fileno())
        if path.exists():
            shutil.copymode(path, tmp_name)
        else:
            os.chmod(tmp_name, 0o666 & ~_current_umask())
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def regenerate_fingerprint(spec: DatasetSpec, *, timeout: int = 60) -> dict:
    """Run ``spec``'s probe and overwrite its committed baseline with the result.

    Goes through the same ``_invoke_with_timeout`` / ``_coerce_fingerprint`` path as a
    drift check, and -- like that path -- any exception the probe raises, not only a
    ``DriftProbeError``, is re-raised as one. So a probe that hangs, raises ANYTHING, or
    returns a non-mapping value fails here with a ``DriftProbeError`` and the file on
    disk is not touched. Before #361 the CLI called ``spec.drift_probe()`` bare and
    ``write_text``'d whatever came back, which is how a corrupt baseline could be
    committed by the unattended bot in the first place.
    """
    try:
        observed = _invoke_with_timeout(spec.drift_probe, timeout=timeout)
    except DriftProbeError:
        raise
    except Exception as exc:  # noqa: BLE001
        raise DriftProbeError(f"probe raised {type(exc).__name__}: {exc}") from exc
    write_fingerprint(Path(spec.test_paths.drift_fingerprint), observed)
    return observed


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
