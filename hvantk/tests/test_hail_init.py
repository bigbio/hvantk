import os
import hail as hl
import pytest
from hvantk.core.utils.hail_context import (
    init_hail,
    hail_initialized,
    get_hail_init_args,
    shutdown_hail,
)
from unittest.mock import patch
import logging

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_HAIL_TESTS") != "1",
    reason="Hail tests require RUN_HAIL_TESTS=1",
)


@pytest.mark.hail
def test_init_hail_idempotent():
    shutdown_hail()  # ensure clean state
    call_count = {"n": 0}
    real_init = hl.init

    def fake_init(**kwargs):
        call_count["n"] += 1
        return real_init(**kwargs)

    with patch("hail.init", side_effect=fake_init):
        assert not hail_initialized()
        init_hail(default_reference="GRCh38")
        assert hail_initialized()
        assert call_count["n"] == 1
        init_hail(default_reference="GRCh38")
        assert call_count["n"] == 1  # still 1
        # conflicting args -> warning but no new init
        init_hail(log="/tmp/another.log")
        assert call_count["n"] == 1
    args = get_hail_init_args()
    assert args.get("default_reference") == "GRCh38"
    shutdown_hail()
    assert not hail_initialized()


@pytest.mark.hail
def test_reinit_after_shutdown():
    shutdown_hail()
    call_count = {"n": 0}
    real_init = hl.init

    def fake_init(**kwargs):
        call_count["n"] += 1
        return real_init(**kwargs)

    with patch("hail.init", side_effect=fake_init):
        init_hail()
        init_hail()
        assert call_count["n"] == 1
        shutdown_hail()
        init_hail()
        assert call_count["n"] == 2  # second init after shutdown
    shutdown_hail()


@pytest.mark.hail
def test_hail_matrix_generation():
    """
    Test if Hail is installed and working correctly.

    This test uses the session-scoped hail_session fixture for Hail initialization.
    """
    # Create a MatrixTable using Hail's balding_nichols_model
    mt = hl.balding_nichols_model(n_populations=3, n_samples=10, n_variants=100)
    logging.info(mt.show())
    assert mt.count() == (100, 10)
