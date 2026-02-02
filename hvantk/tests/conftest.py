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
    """Auto-apply markers based on filename suffix and ensure hail tests use hail_session."""
    for item in items:
        fname = pathlib.Path(item.fspath).name
        # Apply suffix-based markers
        for suffix, marker in (
            _SUFFIX_TO_MARK.items() if False else _SUFFIX_TO_MARK.items()
        ):
            if fname.endswith(suffix):
                item.add_marker(marker)
        # Ensure hail tests use the hail_session fixture (single-shot Hail startup)
        if item.get_closest_marker("hail") is not None:
            item.add_marker(pytest.mark.usefixtures("hail_session"))


# ----------------------------
# Session-scoped Hail init
# ----------------------------
@pytest.fixture(scope="session")
def hail_session(tmp_path_factory, request):
    """
    Initialize Hail once per test session (and per xdist worker) to avoid repeated startup cost
    and spurious re-init errors. Auto-applied to tests marked with @pytest.mark.hail.

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

    # Detect if already initialized
    started_here = False
    already_initialized = False
    try:
        backend = hl.current_backend()
        if backend is not None:
            already_initialized = True
    except Exception:
        already_initialized = False

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
