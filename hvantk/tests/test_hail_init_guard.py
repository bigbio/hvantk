"""Guard logic for init_hail's "is Hail already running?" check.

Deliberately NOT in ``test_hail_init.py``: that module is gated behind
``RUN_HAIL_TESTS=1`` because its tests start a real Hail session. These do not --
they exercise the detection branch with monkeypatched state, so they run in the
default fast suite where a regression would actually be noticed.
"""

import pytest


def test_detection_helper_does_not_itself_initialize_hail(monkeypatch):
    """Detecting an existing backend must not CREATE one.

    ``hl.current_backend()`` is not a predicate: it resolves through ``Env.hc()``,
    which prints "Initializing Hail with default parameters..." and calls ``init()``
    when no backend exists. Using it as the guard makes the guard self-fulfilling --
    it starts Hail with defaults, then reports "already running" and silently drops
    the caller's kwargs.

    Observed consequence: every ``init_hail(tmp_dir=..., local_tmpdir=...)`` in a
    fresh process was a no-op, so on a multi-node Spark cluster Hail spilled to
    node-local ``/tmp`` and executors on other nodes died with
    ``FileNotFoundException: file:/tmp/aggregate_intermediates/...``.
    """
    from hail.utils.java import Env

    from hvantk.core.utils.hail_context import _hail_already_initialized

    monkeypatch.setattr(Env, "_hc", None, raising=False)
    assert _hail_already_initialized() is False

    monkeypatch.setattr(Env, "_hc", object(), raising=False)
    assert _hail_already_initialized() is True


def test_init_hail_applies_kwargs_when_no_backend_exists(monkeypatch):
    """A fresh process must forward the caller's kwargs to hl.init, not discard them."""
    import hail as hl
    from hail.utils.java import Env

    import hvantk.core.utils.hail_context as ctx

    monkeypatch.setattr(ctx, "_HAIL_INITIALIZED", False, raising=False)
    monkeypatch.setattr(ctx, "_HAIL_INIT_ARGS", {}, raising=False)
    monkeypatch.setattr(Env, "_hc", None, raising=False)

    seen = {}
    monkeypatch.setattr(hl, "init", lambda **kw: seen.update(kw))

    ctx.init_hail(tmp_dir="/shared/tmp", local_tmpdir="/shared/tmp")
    assert seen == {"tmp_dir": "/shared/tmp", "local_tmpdir": "/shared/tmp"}


def test_adoption_still_works_when_a_backend_really_exists(monkeypatch):
    """The original behaviour must survive: adopt a session started elsewhere."""
    import hail as hl
    from hail.utils.java import Env

    import hvantk.core.utils.hail_context as ctx

    monkeypatch.setattr(ctx, "_HAIL_INITIALIZED", False, raising=False)
    monkeypatch.setattr(ctx, "_HAIL_INIT_ARGS", {}, raising=False)
    monkeypatch.setattr(Env, "_hc", object(), raising=False)

    called = []
    monkeypatch.setattr(hl, "init", lambda **kw: called.append(kw))

    ctx.init_hail(tmp_dir="/shared/tmp")
    assert called == []  # adopted, not re-initialized
    assert ctx.hail_initialized() is True
    # kwargs were NOT applied, so they must not be reported as if they had been
    assert ctx.get_hail_init_args() == {}


@pytest.mark.parametrize("sentinel", [None, object()])
def test_helper_is_side_effect_free(monkeypatch, sentinel):
    """Calling the helper must never mutate Env._hc, whatever its state."""
    from hail.utils.java import Env

    from hvantk.core.utils.hail_context import _hail_already_initialized

    monkeypatch.setattr(Env, "_hc", sentinel, raising=False)
    _hail_already_initialized()
    assert Env._hc is sentinel
