"""Global Hail initialization context utilities.

Provides a single, idempotent, thread-safe initializer for Hail so multiple
streamers or commands do not attempt concurrent hl.init() calls (which Hail
does not support cleanly).

Public API:
- init_hail(**kwargs): initialize Hail once (idempotent). Subsequent calls are
  no-ops unless conflicting kwargs are passed (then a warning is logged).
- hail_initialized() -> bool: whether init_hail has successfully run.
- get_hail_init_args() -> dict: kwargs captured from the first initialization.
- shutdown_hail(): optional explicit stop; after calling, re-init is allowed.

Environment variables:
- HVANTK_SKIP_HAIL_INIT: if set to a non-empty value, init_hail() will skip
  calling hl.init() (assumes caller/environment already initialized Hail).
"""

from __future__ import annotations

import os
import threading
import logging
from typing import Dict, Any

import numpy as np

# Hail 0.2.x references ``np.bool``, removed in NumPy 1.24+. Restore the alias
# before importing Hail so its internal modules don't AttributeError on import.
if not hasattr(np, "bool"):
    np.bool = np.bool_  # type: ignore[attr-defined]

import hail as hl  # noqa: E402

_logger = logging.getLogger(__name__)

_INIT_LOCK = threading.Lock()
_HAIL_INITIALIZED = False
_HAIL_INIT_ARGS: Dict[str, Any] = {}


def _hail_already_initialized() -> bool:
    """Whether a Hail backend exists, WITHOUT creating one.

    ``hl.current_backend()`` looks like the natural probe but is not a predicate: it
    resolves through ``Env.hc()``, which calls ``hl.init()`` with defaults when no
    backend exists (printing "Initializing Hail with default parameters..."). Used as a
    guard it is self-fulfilling -- it starts Hail with the wrong settings, then reports
    "already running", so the caller's kwargs are dropped on the floor.

    That silently disabled ``init_hail(tmp_dir=..., local_tmpdir=...)`` in every fresh
    process. On a multi-node Spark cluster the consequence is not cosmetic: Hail spills
    to node-local ``/tmp``, and executors on other nodes fail with
    ``FileNotFoundException: file:/tmp/aggregate_intermediates/...``.

    Hail exposes no public "is initialized" predicate, so this reads ``Env._hc``
    directly. It is a plain attribute check with no side effects, which is exactly the
    property the guard needs.
    """
    try:
        from hail.utils.java import Env

        return Env._hc is not None
    except Exception:  # pragma: no cover - hail missing or internals moved
        return False


def hail_initialized() -> bool:
    return _HAIL_INITIALIZED


def get_hail_init_args() -> Dict[str, Any]:
    return dict(_HAIL_INIT_ARGS)


def init_hail(**kwargs) -> None:
    """Initialize Hail once (thread-safe, idempotent).

    Subsequent calls with the same (or no) kwargs are ignored. If different
    kwargs are provided after initialization, a warning is logged and the
    existing session is kept.
    """
    global _HAIL_INITIALIZED, _HAIL_INIT_ARGS

    if os.environ.get("HVANTK_SKIP_HAIL_INIT"):
        if not _HAIL_INITIALIZED:
            _logger.info(
                "Skipping hl.init() due to HVANTK_SKIP_HAIL_INIT environment variable"
            )
            _HAIL_INITIALIZED = True  # Treat as initialized to avoid later attempts
        return

    if _HAIL_INITIALIZED:
        if kwargs and any(_HAIL_INIT_ARGS.get(k) != v for k, v in kwargs.items()):
            _logger.warning(
                "Hail already initialized with %s; ignoring conflicting re-init kwargs %s",
                _HAIL_INIT_ARGS,
                kwargs,
            )
        return

    with _INIT_LOCK:
        # Re-check inside lock
        if _HAIL_INITIALIZED:
            if kwargs and any(_HAIL_INIT_ARGS.get(k) != v for k, v in kwargs.items()):
                _logger.warning(
                    "(post-lock) Hail already initialized with %s; ignoring conflicting re-init kwargs %s",
                    _HAIL_INIT_ARGS,
                    kwargs,
                )
            return
        # A Hail backend may already exist from a path other than init_hail()
        # — e.g. a raw hl.init() in a notebook or a test fixture. Calling
        # hl.init() again raises Hail's "already initialized" error, so adopt
        # the existing session instead of re-initializing.
        # NB: probe with _hail_already_initialized(), never hl.current_backend() —
        # the latter CREATES a default-configured backend rather than reporting one.
        if _hail_already_initialized():
            _HAIL_INITIALIZED = True
            # The kwargs were NOT applied (Hail was initialized elsewhere), so
            # record an empty dict rather than the requested kwargs — keeps
            # get_hail_init_args() honest and lets later calls still surface
            # conflicting-kwargs warnings.
            _HAIL_INIT_ARGS = {}
            if kwargs:
                _logger.warning(
                    "Hail already initialized outside init_hail(); "
                    "cannot apply init kwargs %s",
                    kwargs,
                )
            else:
                _logger.info(
                    "Adopting Hail backend already initialized outside init_hail()"
                )
            return
        hl.init(**kwargs)
        _HAIL_INITIALIZED = True
        _HAIL_INIT_ARGS = dict(kwargs)
        _logger.info("Hail initialized (once) with args: %s", _HAIL_INIT_ARGS)


def shutdown_hail() -> None:
    """Explicitly stop Hail session and reset guard.

    Not called automatically to avoid stopping a shared global context.
    """
    global _HAIL_INITIALIZED, _HAIL_INIT_ARGS
    if not _HAIL_INITIALIZED:
        return
    try:
        hl.stop()
        _logger.info("Hail stopped via shutdown_hail()")
    except Exception as e:  # pragma: no cover - defensive
        _logger.warning("Error while stopping Hail: %s", e)
    finally:
        _HAIL_INITIALIZED = False
        _HAIL_INIT_ARGS = {}
