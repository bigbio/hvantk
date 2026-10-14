# Global pytest configuration for tests
# - Use a non-interactive matplotlib backend for faster, headless testing

import os
import pathlib
import matplotlib
import pytest

# Force a headless backend for any plotting tests
matplotlib.use("Agg")

# Make pytest output shorter for CLI commands
os.environ.setdefault("PYTHONWARNINGS", "ignore")


# ----------------------------
# Auto-marking by filename
# ----------------------------
_SUFFIX_TO_MARK = {
    "_hail.py": pytest.mark.hail,
    "_network.py": pytest.mark.network,
    "_slow.py": pytest.mark.slow,
    "_integration.py": pytest.mark.integration,
}


def pytest_collection_modifyitems(config, items):
    """Auto-apply markers based on filename suffix."""
    for item in items:
        fname = pathlib.Path(item.fspath).name
        # Apply suffix-based markers
        for suffix, marker in (
            _SUFFIX_TO_MARK.items() if False else _SUFFIX_TO_MARK.items()
        ):
            if fname.endswith(suffix):
                item.add_marker(marker)


# ----------------------------
# Session-scoped Hail init
# ----------------------------
@pytest.fixture(scope="session")
def hail_session(tmp_path_factory, request):
    """
    Initialize Hail once per test session (and per xdist worker) to avoid repeated startup cost
    and spurious re-init errors. Set up for every test marked with @pytest.mark.hail by
    ``_hail_session_for_hail_tests`` below.

    Honors optional env vars to tune Spark:
      - PYSPARK_MASTER (default: local[*])
      - SPARK_DRIVER_MEMORY (default: 2g)
      - SPARK_EXECUTOR_MEMORY (default: 2g)
      - SPARK_SQL_SHUFFLE_PARTITIONS (default: 4)
      - HAIL_DEFAULT_REFERENCE (optional)
    """
    import importlib

    # Determine worker id (xdist) for isolated tmp dirs
    worker_id = "master"
    try:
        # present when xdist plugin runs
        workerinput = getattr(request.config, "workerinput", None)
        if isinstance(workerinput, dict):
            worker_id = workerinput.get("workerid", "master")
    except Exception:
        pass

    # Unique tmp dir per worker
    base_tmp = tmp_path_factory.mktemp(f"hail_tmp_{worker_id}")
    spark_local_dir = base_tmp / "spark_local"
    spark_local_dir.mkdir(exist_ok=True, parents=True)

    # Lazy import hail
    hl = importlib.import_module("hail")

    # Detect if already initialized. Not with `hl.current_backend()`: when no backend
    # exists it starts Hail with default settings, so the init below never ran (see
    # `_hail_already_initialized`).
    from hvantk.core.utils.hail_context import _hail_already_initialized

    started_here = False
    already_initialized = _hail_already_initialized()

    if not already_initialized:
        master = os.getenv("PYSPARK_MASTER", "local[*]")
        driver_mem = os.getenv("SPARK_DRIVER_MEMORY", "2g")
        exec_mem = os.getenv("SPARK_EXECUTOR_MEMORY", "2g")
        shuffle_parts = os.getenv("SPARK_SQL_SHUFFLE_PARTITIONS", "4")
        default_ref = os.getenv("HAIL_DEFAULT_REFERENCE")

        spark_conf = {
            "spark.ui.enabled": "false",
            "spark.ui.showConsoleProgress": "false",
            "spark.sql.shuffle.partitions": str(shuffle_parts),
            "spark.driver.memory": str(driver_mem),
            "spark.executor.memory": str(exec_mem),
            "spark.local.dir": str(spark_local_dir),
            "spark.driver.maxResultSize": "0",
        }

        init_kwargs = {
            "app_name": f"hvantk-tests-{worker_id}",
            "master": master,
            "tmp_dir": str(base_tmp),
            "spark_conf": spark_conf,
        }
        if default_ref:
            init_kwargs["default_reference"] = default_ref

        hl.init(**init_kwargs)
        started_here = True

        # Optional tiny smoke-check to fail fast
        try:
            _ = hl.eval(hl.len(hl.array([1, 2, 3])))
        except Exception:
            # If a transient error happens here, ensure we stop to avoid hanging sessions
            hl.stop()
            raise

    # Yield the hail module for convenience (optional to use in tests)
    yield hl

    # Only stop if we started the session here
    if started_here:
        try:
            hl.stop()
        except Exception:
            pass


@pytest.fixture(autouse=True)
def _hail_session_for_hail_tests(request):
    """Set up ``hail_session`` for each test marked ``hail``.

    Requested at setup because pytest fixes a test's fixtures at collection: a
    ``usefixtures`` mark added in ``pytest_collection_modifyitems`` is ignored.
    Each plugin's ``tests/conftest.py`` re-exports this with ``hail_session``.
    """
    if request.node.get_closest_marker("hail") is not None:
        request.getfixturevalue("hail_session")


# ----------------------------
# Skill snapshot regeneration flag
# ----------------------------
# ``pytest_addoption`` for --regenerate-snapshots lives in the repo-root conftest.py,
# not here: this file is only loaded for paths under hvantk/tests/, so registering it
# here left the option undefined when a plugin's own suite was run on its own (pytest
# conftest discovery walks ancestors, and hvantk/skills/ is a sibling of hvantk/tests/).
# The fixture below stays here, because each plugin's tests/conftest.py re-exports it.
@pytest.fixture
def regenerate_snapshots(request) -> bool:
    return request.config.getoption("--regenerate-snapshots")
