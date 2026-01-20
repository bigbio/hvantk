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
import hail as hl

_logger = logging.getLogger(__name__)

_INIT_LOCK = threading.Lock()
_HAIL_INITIALIZED = False
_HAIL_INIT_ARGS: Dict[str, Any] = {}


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
