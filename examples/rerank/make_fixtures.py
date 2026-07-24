"""Regenerate the synthetic fixtures under ``examples/rerank/data/``.

The data is invented, not derived from any cohort: this example exists so a reviewer can
run ``hvantk rerank`` end-to-end and get a deterministic table, not to reproduce a
scientific result. Everything is seeded (``default_rng(0)``), so re-running this rewrites
the same bytes.

The generated signal is deliberately weak-but-real. Both feature axes carry information
about the label, and expression carries a little less than constraint, so the drop-one
ablation reports a non-trivial ordering instead of two identical rows. Do not read the
resulting AUCs as anything but a smoke test.

Run from the repository root::

    python examples/rerank/make_fixtures.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "data"

N_GENES = 150
POSITIVE_RATE = 0.30
SEED = 0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    genes = [f"GENE{i:03d}" for i in range(N_GENES)]
    y = (rng.random(N_GENES) < POSITIVE_RATE).astype(int)

    # Prior: a burden-style p-value, folded into the single cohort table below (a cohort
    # manifest can only point at one gene-level table -- see data/cohort.yaml). Passed
    # through to the output as `prior_stat`; it is NOT a model feature (engine.py builds
    # feat_cols from the feature matrix only). Positives get smaller p-values so the
    # prior is a plausible starting ranking.
    minp = np.where(
        y == 1,
        rng.uniform(1e-6, 5e-2, N_GENES),
        rng.uniform(1e-3, 1.0, N_GENES),
    )

    # Axis 1 -- constraint. Listed first in config.yaml, so the Evaluator treats it as the
    # ablation baseline.
    pd.DataFrame(
        {
            "gene": genes,
            "mis_z": y * 1.0 + rng.normal(0, 0.6, N_GENES),
            "pli": np.clip(y * 0.4 + rng.uniform(0, 1, N_GENES), 0, 1),
        }
    ).to_csv(OUT / "constraint.tsv", sep="\t", index=False, float_format="%.6g")

    # Axis 2 -- expression. A second axis is required for the drop-one ablation to emit
    # per-axis rows at all; with a single axis the Evaluator reports only the baseline.
    pd.DataFrame(
        {
            "gene": genes,
            "dev_expr": y * 0.7 + rng.normal(0, 0.8, N_GENES),
        }
    ).to_csv(OUT / "expression.tsv", sep="\t", index=False, float_format="%.6g")

    # Labels: the positive set. Coverage is 1.0 (every labelled gene is in the matrix),
    # so the default min_label_coverage=0.5 gate passes without being overridden.
    (OUT / "labels.txt").write_text(
        "\n".join(g for g, lab in zip(genes, y) if lab == 1) + "\n"
    )

    # Cohort: one gene-level table carrying both the prior statistic (`minp`) and the
    # case/control architecture columns CaseControlArchitectureAudit requires --
    # data/cohort.yaml declares `minp` as its prior and `n_case_var`/`conc`/`driver_af`
    # as its "architecture" axis. Seed each flag branch so the advisory columns in the
    # output are non-empty and a reviewer can see the audit working (audit.py):
    #   n_case_var <= 2                  -> insufficient_data
    #   conc >= 0.6 and driver_af > 5e-5 -> recurrent_variant
    #   driver_af > 1e-3                 -> common_driver   (highest precedence)
    n_case_var = rng.integers(3, 40, N_GENES)
    conc = rng.uniform(0.0, 0.5, N_GENES)
    driver_af = rng.uniform(0, 5e-5, N_GENES)

    n_case_var[:4] = rng.integers(0, 3, 4)  # insufficient_data
    conc[10:14] = rng.uniform(0.6, 0.95, 4)  # recurrent_variant
    driver_af[10:14] = rng.uniform(1e-4, 9e-4, 4)
    driver_af[20:23] = rng.uniform(2e-3, 8e-3, 3)  # common_driver

    pd.DataFrame(
        {
            "gene": genes,
            "minp": minp,
            "n_case_var": n_case_var,
            "conc": conc,
            "driver_af": driver_af,
        }
    ).to_csv(OUT / "cohort.tsv", sep="\t", index=False, float_format="%.6g")

    print(f"wrote 4 fixtures to {OUT}")
    print(
        f"  {N_GENES} genes, {int(y.sum())} positives, {N_GENES - int(y.sum())} negatives"
    )


if __name__ == "__main__":
    main()
