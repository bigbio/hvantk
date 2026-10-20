# Re-rank example

Walks through [`examples/rerank/`](https://github.com/bigbio/hvantk/tree/main/examples/rerank),
a complete, runnable `hvantk rerank` run on small synthetic fixtures. It needs no Hail, no
network access and no real cohort -- the whole thing runs in a few seconds -- so it is the
fastest way to see the re-ranking workflow end to end before pointing it at real data. For
the full option reference, the config schema and the opt-in controls (blocked
cross-validation, the multiplicity correction, seed sweeps), see
[Re-rank](../tools/rerank.md).

## What the fixtures are

`examples/rerank/make_fixtures.py` generates everything under `data/` with a fixed seed
(`numpy.default_rng(0)`), so re-running it rewrites the same bytes:

- 150 genes, 40 of them positive (the label set in `data/labels.txt`).
- `data/constraint.tsv` (`mis_z`, `pli`) and `data/expression.tsv` (`dev_expr`) -- two
  feature axes, both carrying real signal about the label, with `expression` deliberately
  weaker than `constraint` so the ablation produces a non-trivial, not-quite-conclusive
  result rather than two identical rows.
- `data/cohort.tsv` + `data/cohort.yaml` -- the cohort manifest. It carries the prior
  statistic (`minp`, a burden-style p-value) and the three case/control architecture
  columns (`n_case_var`, `conc`, `driver_af`) the advisory audit reads, in one table,
  because a cohort manifest can only point at a single gene-level table. The fixture
  generator deliberately seeds every branch of the audit -- 4 genes with too few case
  variants, 4 with a recurrent shared driver, 3 with a common driver -- so the flag columns
  in the output are never empty.

Everything here is invented. Read the numbers below as a smoke test that the pipeline
works, not as evidence about any biology.

## The config

`config.yaml` names the cohort manifest, lists the two feature axes (the first,
`constraint`, is the ablation baseline), and points at the label file:

```yaml
name: rerank-example
cohort: examples/rerank/data/cohort.yaml
features:
  - name: constraint
    path: examples/rerank/data/constraint.tsv
  - name: expression
    path: examples/rerank/data/expression.tsv
labels:
  path: examples/rerank/data/labels.txt
```

Paths resolve against the process working directory, not against `config.yaml` itself, so
run the command from the repository root.

## Running it

```bash
poetry install --extras ml
eval "$(poetry env activate)"

mkdir -p out
hvantk rerank -c examples/rerank/config.yaml -o out/rerank_scores.tsv
```

This prints:

```text
Wrote 150 genes -> out/rerank_scores.tsv
Credibility: ROBUST=60 / INTERMEDIATE=90 (over 150 scored genes)
Audit: 11 gene(s) FLAGGED for review (advisory; ranking not overridden).
       Reasons: {'insufficient_data': 4, 'recurrent_variant': 4, 'common_driver': 3}
Per-axis ablation (delta-AUC over constraint):
    family   auc   d_lo  d_md  d_hi
constraint 0.846  0.000 0.000 0.000
expression 0.864 -0.022 0.018 0.066
```

A committed reference lives at `examples/rerank/results/rerank_scores.tsv`; `diff` your
output against it to confirm your environment reproduces it:

```bash
diff out/rerank_scores.tsv examples/rerank/results/rerank_scores.tsv
```

## Reading the output

`out/rerank_scores.tsv` has one row per gene: `gene`, `prior_stat` (the cohort's burden
p-value, carried through, not a model feature), `score` (the calibrated credibility
score), `score_percentile`, `tier`, `verdict`, `flag`, `flag_reason`, and `y` (the label).
See [Output columns](../tools/rerank.md#output-columns) for the full reference.

Two things worth noticing on this run:

- **The audit is advisory.** 11 genes carry a non-empty `flag_reason`
  (`insufficient_data`, `recurrent_variant` or `common_driver`), but a flagged gene keeps
  its `score` and its place in the ranking -- the audit flags genes for a human to look at,
  it never changes the model's answer.
- **The ablation's confidence interval spans zero.** `expression` adds `+0.018` AUC over
  `constraint` alone (`d_md`), but `d_lo` is negative and `d_hi` is positive -- on this small
  a fixture, that gain is not distinguishable from noise. That is the honest, expected
  outcome for 150 genes with a deliberately weak second axis, not a bug in the example.

## Feature selection (Python API only)

The shipped `config.yaml` has no `selection:` key, so a CLI run uses every column of every
axis -- fine for the four columns here, unworkable for a real matrix such as dbNSFP's ~45
correlated predictors. The Python API's `Config.selection` turns on a declarative, within-axis
selection stage (`Config` is a frozen dataclass, so derive a changed copy with
`dataclasses.replace` rather than assigning to an existing one):

```python
import dataclasses

from hvantk.algorithms.rerank.engine import rerank_arms
from hvantk.algorithms.rerank.selection import SelectionPolicy

config = dataclasses.replace(
    config,
    selection=SelectionPolicy(),
    feature_provenance={"mis_z": frozenset(), "pli": frozenset(), "dev_expr": frozenset()},
    label_provenance=frozenset(),  # this example's labels are not derived from a curated DB
)

arms = rerank_arms(config)
arms["clean"].metrics.auc
```

See the full walkthrough -- the circularity arms, the RFECV wrapper, the nesting guarantee
-- in [`examples/rerank/README.md`](https://github.com/bigbio/hvantk/tree/main/examples/rerank#feature-selection-python-api-only).
The CLI does not expose this yet.

## See also

- [Re-rank](../tools/rerank.md) -- the full CLI/API reference: every option, the config
  schema, and the blocked-CV, multiplicity and seed-sweep controls.
- [External cohorts](../tools/cohort.md) -- the cohort manifest format `data/cohort.yaml`
  here is an example of.
