# hvantk/tests/rerank/conftest.py
# Limit OpenMP / BLAS thread count so GBM cross-validation tests don't
# spawn 100+ threads on many-core HPC nodes and get OOM-killed.
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
