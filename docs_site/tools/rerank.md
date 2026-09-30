# ReRank: Gene Re-ranking by Multi-Omic Credibility

`hvantk rerank` re-ranks a gene burden or prior statistic against one or more gene-keyed feature axes -- constraint, expression, or any other numeric table keyed by gene -- by fitting a calibrated, out-of-fold gradient-boosted model. It assigns credibility tiers and flags cohort artifacts for review (advisory only: the audit never overrides the ranking). The feature tables and the cohort's prior statistic are supplied by the caller through a declarative YAML config; the engine itself is disease-agnostic.

## Overview

A per-axis result -- "axis X adds +0.03 AUC over the constraint baseline" -- is only as informative as the questions that have been asked of it. Three questions recur, and each is answered by its own module rather than by inspection of the output table alone:

| Risk | Question | Answered by |
|---|---|---|
| Circularity | Was this feature trained on the same evidence the label was derived from? | `provenance.py` (`Config.feature_provenance` / `label_provenance`, `rerank_arms`) |
| Presence leakage | Does a column's *missingness*, not its value, predict the label? | `leakage.py` (`Config.leakage`, `LeakagePolicy`) |
| Multiplicity | Several axes were offered and the best one reported -- could that delta have arisen by chance? | `nulls.py` (`Config.nulls`, `NullConfig`, `--n-perm`) |

This page documents the multiplicity control on the right of that last row, plus two further controls that ship alongside it and ask a related but distinct question -- not *is the delta real*, but *is the delta an artefact of one arbitrary modelling choice*:

- **Paralogue-blocked cross-validation** -- does the delta survive when a gene family is not allowed to straddle a cross-validation fold?
- **Multi-seed evaluation** -- does the delta survive a change of the one cross-validation partition every earlier run held fixed?

All three are opt-in CLI flags, and every one of them defaults to today's behaviour: a bare `hvantk rerank -c config.yaml -o out.tsv` is unaffected by any of the flags below unless you pass them.

## Quick Start

```bash
hvantk rerank -c config.yaml -o out.tsv
```

This is the minimal invocation: a config naming a cohort manifest, one or more feature tables, and a label gene list (see **The config YAML** below). It writes a ranked TSV to `out.tsv` and prints a short console summary:

```text
Wrote 1500 genes -> out.tsv
Credibility: ROBUST=210 / INTERMEDIATE=1290 (over 1500 scored genes)
Audit: 3 gene(s) FLAGGED for review (advisory; ranking not overridden). Reasons: {...}
Per-axis ablation (delta-AUC over constraint):
    family   auc  d_lo  d_md  d_hi
constraint 0.706 0.000 0.000 0.000
expression 0.744 0.041 0.058 0.083
```

The ablation table is the per-axis result: out-of-fold AUC of baseline+axis, and a paired gene-resampling bootstrap interval (`d_lo`/`d_md`/`d_hi`) on the delta over the baseline alone. See **Output columns** for the full column reference, including the columns the controls below add.

## The config YAML

`-c`/`--config` reads exactly five top-level keys; any other key is rejected before any work starts (a leftover `prior:` block from before the cohort-manifest migration gets its own hint in the error message).

| Key | Required | Meaning |
|---|---|---|
| `name` | no | A label for the run (default `"rerank"`); has no effect on scoring. |
| `cohort` | yes | Path to a cohort manifest YAML -- the single source of the prior statistic and, when it declares one, the columns the case/control architecture audit reads. |
| `features` | yes | A list of `{name, path}` gene-keyed tables (`.parquet`/`.tsv`/`.csv`: a `gene` column plus one or more numeric feature columns). The **first** entry is the baseline every other axis is measured against. |
| `labels` | yes | `{path}` to a file of gene symbols, one per line -- the positive set. |
| `min_label_coverage` | no | Minimum fraction of label-positive genes that must join the feature matrix (default `0.5`); below it the run aborts rather than silently scoring a mismatched gene universe. |

```yaml
name: my-rerank-run
cohort: cohort.yaml
features:
  - name: constraint
    path: constraint.parquet
  - name: expression
    path: expression.tsv
labels:
  path: labels.txt
min_label_coverage: 0.5
```

`cohort.yaml` is a separate, smaller manifest: the gene-key column, the prior statistic and its direction, and optionally which columns feed the architecture audit.

```yaml
name: my-cohort
key: symbol
key_column: gene
table: cohort.tsv
prior:
  column: p
  direction: lower_is_better
```

This section describes only the fields a rerank config needs, not a full reference for the cohort manifest itself (see `hvantk/algorithms/cohort/spec.py` for the complete `CohortManifest` shape: `key`, `key_column`, `table`, `prior.{column,direction}`, optional `cohort_axes`/`labels`) -- no dedicated cohort-manifest reference page exists under `docs_site/` yet.

## Controls

Three controls, six flags, all opt-in:

| Flag | Controls | Default |
|---|---|---|
| `--n-perm INTEGER` | Multiplicity | `0` (off) |
| `--null-out PATH` | Multiplicity | none |
| `--blocks PATH` | Blocked CV | none (unblocked) |
| `--max-block-frac FLOAT` | Blocked CV | `0.10` |
| `--seed INTEGER` | Seeds | `42` |
| `--seed-sweep INTEGER` | Seeds | `1` |

### Multiplicity -- `--n-perm`, `--null-out`

Offering several feature axes and reporting the best one's delta-AUC is itself a selection procedure, and the ablation table's bootstrap interval has nothing to say about it -- it describes one axis in isolation. `--n-perm` builds two permutation nulls that do:

- **Per-axis** -- permute the label, refit the baseline and baseline+*this* axis, record the delta. Correct when the axis was pre-specified.
- **Selected-maximum** -- permute, refit, add *each* offered axis, record the *largest* delta. Correct when the axis is reported *because* it came top over the other candidates. Only this one is a multiplicity correction: the maximum over several candidate axes is stochastically larger than any single one of them, so a per-axis null can sit near zero while the selected-maximum null does not.

Both nulls refit the baseline on every permutation -- scoring a permuted-label model against the unpermuted baseline would measure the permutation, not the axis. P-values use `p = (1 + #{null >= observed}) / (1 + n_perm)` (Phipson & Smyth 2010): a naive `(null >= observed).mean()` reports `p = 0` for a value no permutation reached, which a finite permutation set can never actually license. `p_selected_max` is the single-step max-T (Westfall-Young) adjusted p-value for that axis's observed delta. The observed deltas themselves come from the null's own scorer, unrounded -- **not** the ablation table's `d_md` (a bootstrap *median*, rounded to 3 decimal places: a different statistic, computed a different way).

> **What null hypothesis this tests.** Permuting the label tests the *global* null that no feature, baseline included, carries label information -- not the *conditional* null that this one axis adds nothing given the baseline already in the model. With an informative baseline, the permutation null of the delta is wider than the delta's true sampling spread under the conditional null, so the test is valid but conservative: failing to clear the selected-maximum null is weak evidence that the axis adds nothing, not proof that it does not.

That conservatism assumes genes are exchangeable under the label permutation. Paralogues share sequence, constraint, expression and disease status, so when labels and features cluster by gene family that assumption breaks and the *default* (unblocked) null becomes **anti-conservative** instead -- understating its own spread -- regardless of whether `--blocks` happens to be set. Pass `--blocks` whenever that clustering is doubtful; when it is set, the null permutes labels **by block** (whole blocks swapped among blocks of the same size, then shuffled within each block -- Winkler et al. 2015, *NeuroImage*), which is what blocking requires of the null too. A block whose size is unique in the universe has no same-size sibling to swap with, so it is only ever shuffled within itself -- a real power caveat for that block, not a bug.

The null is tied to one **control setting** -- the leakage and selection policies, the provenance arm, the exact baseline and candidate columns, the fold count actually used, and which blocking (if any) the folds used -- and refuses to answer about a delta computed under a different one, raising `ControlSettingMismatch` rather than silently comparing incomparable numbers. The ablation (and so the null) always scores 5 folds regardless of `Config.folds`, which governs only the headline score. The null's `observed` deltas are computed at `--seed` alone (never averaged over `--seed-sweep`); see **Seeds** below for why the two are reported separately.

```bash
hvantk rerank -c config.yaml -o out.tsv --n-perm 200 --null-out null_summary.tsv
```

adds, after the ablation table:

```text
Multiplicity: 1 axis searched over 200 permutations; selected-maximum null median +0.0410
      axis  observed  selmax_median  p_selected_max
expression     0.061          0.041            0.07
```

`--null-out` writes the full per-axis summary (see **Output columns**); without it, the console shows the header line and the four columns above only. `--null-out` requires `--n-perm`: passing `--null-out` without `--n-perm` is rejected before any work starts (`Error: --null-out describes a permutation null; pass --n-perm N too.`), as is a negative `--n-perm` (`Error: --n-perm must be >= 0; got -5.`, for example) -- both a plain `click` failure, never a traceback. Via the Python API, `NullConfig()` alone (with no `n_perm=`) defaults to 200 permutations; the CLI's own default is `0` (off) so that a plain `hvantk rerank` run costs nothing extra.

### Blocked CV -- `--blocks`, `--max-block-frac`

Every number `hvantk rerank` produces by default comes from plain stratified 5-fold cross-validation. Gene families break the independence that relies on: paralogues share sequence, constraint, expression pattern and disease status, so a family split across train and test lets the model partly recognise a relative rather than generalise. `--blocks` switches to cross-validation folds that keep each gene family whole.

```bash
hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_complete_set.txt
```

`--blocks` takes an HGNC complete-set TSV with `symbol`/`gene_group` columns (and, optionally, `status` -- when present, only `Approved` rows are kept). Genes are matched against it by their **approved HGNC symbol only**; a gene absent from the table, or known to your universe only by an alias or a previous symbol, becomes its own unconstrained singleton block rather than an error -- unless *none* of the universe's genes match, which almost always means the table's gene keys disagree with your universe's (e.g. Ensembl IDs against HGNC symbols) and *is* an error. A match rate below 50% is logged as a warning, which a CLI user may not see the logger output for -- the console's `Blocks:` line below always prints the matched count too, precisely so this degradation has a visible signal.

Blocking uses each gene's **first-listed** HGNC gene group -- in HGNC's own list order, since HGNC does not designate a "primary" one -- rather than connected components over the (multi-membership, pipe-separated) `gene_group` field. Group membership chains transitively, so a union-find closure over shared membership can collapse a large fraction of a gene universe into a single component; blocking on it would then let `StratifiedGroupKFold` place nearly the whole universe in one fold, and the pooled out-of-fold AUC would no longer estimate the same quantity as the unblocked run. The tell-tale of that mistake is a *blocked* AUC scoring **higher** than the random-fold one -- correct blocking is strictly harder than random folds and cannot do that, so treat it as a red flag rather than a happy surprise. First-listed grouping keeps blocks small; the residual leak it accepts -- two genes sharing only a later-listed group -- is a far smaller error than an incomparable estimator.

**The ceiling is a hard abort, not a warning.** `--max-block-frac` (default `0.10`) caps the largest block's share of the universe; going over it fails the run instead of producing a plausible-looking number, because the failure mode is otherwise silent -- the output looks fine, and a counterintuitive AUC is the only tell, and a reader has no particular reason to question it:

```text
$ hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_complete_set.txt
Error: the largest paralogue block (group 'OR6') holds 620/1500 units (41.3%), over the 0.1 ceiling (max_block_frac). StratifiedGroupKFold must place a whole block in one fold, so one fold would hold that entire block and the pooled out-of-fold AUC would not estimate the same quantity as the unblocked run. Raise max_block_frac deliberately if you accept that, or restrict the universe.
```

No output file is written when this happens. `--max-block-frac` only means anything alongside `--blocks`; passing it alone is rejected before any work starts (`Error: --max-block-frac only applies to paralogue-blocked folds; pass --blocks too (or drop it).`). Raise it deliberately (`--max-block-frac 0.5`, say) if a dominant block is expected and accepted for your universe.

Two further mechanics worth knowing about: the grouped splitter is a **correctly stratified** one, not scikit-learn's `StratifiedGroupKFold(shuffle=True)` -- as of scikit-learn 1.7.2 that shuffles the per-group class-count rows before balancing folds but then assigns each group to a fold under its *original*, unshuffled index, so groups stay whole but the folds stop actually being stratified; `hvantk` instead relabels group IDs with a seeded random permutation and calls the unshuffled (correctly balancing) path underneath, which randomises the same tie-breaks the upstream shuffle intended without inheriting its defect, and stays correct even if that upstream behaviour is later fixed. Separately, a class confined to too few blocks to cover the requested fold count is refused with a clear error rather than silently producing a fold with no positive (or no negative) training example, which would otherwise surface only as an opaque scikit-learn warning, or several calls later as "found array with 0 sample(s)".

Blocking is threaded everywhere a fold boundary matters: the outer cross-validation of both the headline score and the ablation, the multi-seed sweep, and -- as described above -- the permutation null. The one documented exception is `CalibratedClassifierCV`'s own inner calibration split, which takes only an integer fold count and so cannot see the blocks; a family can therefore straddle *that* split within one already-blocked outer training fold. The per-axis deltas that are the scientific claim never go through it -- they use the uncalibrated, fully-blocked scorer in `evaluator.py`. Custom API scorers passed directly to `permutation_deltas` that carry no `.groups` attribute are not checked against `blocks=` at all -- only scorers built by the library's own `oof_scorer` declare one.

The console's `Blocks:` line (shown whenever `--blocks` did not abort) reports the blocking, including the matched-gene count that catches a silent key mismatch:

```text
Blocks: 37 paralogue block(s); largest 6 (4.0%); 58 gene(s) in a block of size > 1; 1,438/1,500 genes matched the gene-group table; digest a1b2c3d4e5f6
```

### Seeds -- `--seed`, `--seed-sweep`

`--seed` (default `42`, the historic value) is the one seed that drives the cross-validation partition (both the headline score and the ablation), the gradient-boosted estimator's own randomness, the bootstrap resample, and the permutation null's base seed. Changing it produces a different, equally valid partition -- it is not a tuning knob, and the default is unchanged so every result produced before `--seed` existed still reproduces exactly.

```bash
hvantk rerank -c config.yaml -o out.tsv --seed-sweep 10
```

The gene-resampling bootstrap (`d_lo`/`d_md`/`d_hi`) asks how a delta would move if the gene *sample* moved. It cannot see a second variance component: which genes land in which cross-validation fold. `--seed-sweep N` recomputes each axis's delta under `N` consecutive cross-validation seeds (`--seed`, `--seed` + 1, ...) and reports the **envelope** -- the union of the bootstrap interval and the across-seed range -- as `d_lo_env`/`d_hi_env`, beside the untouched `d_lo`/`d_md`/`d_hi`. The same sweep width travels under a different name at each layer: the CLI's `--seed-sweep`, `Config.seed_sweep`, `Evaluator.evaluate`'s `n_seeds` parameter, and the ablation table's own `n_seeds` column all refer to the identical count -- so a reader who encounters more than one of these names is not looking at two different settings:

```text
Per-axis ablation (delta-AUC over constraint):
    family   auc  d_lo  d_md  d_hi  d_lo_env  d_hi_env  n_seeds
constraint 0.706 0.000 0.000 0.000     0.000     0.000       10
expression 0.744 0.041 0.058 0.083     0.033     0.091       10
```

The envelope is a **union, not a calibrated interval** over both components: the two are not independent draws from one distribution, and combining them properly would need a nested design nobody has run. `d_lo_env`/`d_hi_env` are that union rounded to 3 decimal places, exactly like `d_lo`/`d_hi`. The union is guaranteed never to be narrower than the bootstrap interval alone, but **not** guaranteed to be strictly wider -- that depends on the data: when the across-seed range happens to fall entirely inside `[d_lo, d_hi]`, `d_lo_env == d_lo` and `d_hi_env == d_hi` and no widening is visible, which is itself informative (it says the CV partition was not, on this data, a meaningful source of extra spread), not a sign that the sweep did nothing. It exists so a reader cannot mistake one lucky (or unlucky) partition for a measured effect -- a single cross-validation partition can land anywhere within that spread. With the default `--seed-sweep 1`, the three extra columns are **absent** from the ablation table entirely, so a plain `hvantk rerank` run is byte-for-byte what it always was.

The permutation null (`--n-perm`, above) already recomputes the stratified partition for every permutation, so it already contains partition noise as one of its own components. The envelope is a separate report on the *observed* delta and is never folded into the null on top of that -- doing so would count partition variance twice. For the same reason, the null's `observed` deltas are always computed at `--seed` alone, never averaged over the swept seeds.

`SelectionPolicy.seed` (used only when a feature-selection policy is configured) is independent of `--seed` by design -- it seeds the selection wrapper's own inner cross-validation, not the outer partition `--seed` controls. Turning on `--seed-sweep N` multiplies the ablation's per-axis cost by roughly `N` (one refit per axis per seed); the baseline is fit once per seed and shared across every axis rather than refit per axis, so the extra cost scales with the number of candidate axes times `N`, not more.

## Output columns

### Ranked table (`-o`/`--output`)

| Column | Meaning |
|---|---|
| `gene` | The gene identifier, in the cohort's key space. |
| `prior_stat` | The cohort's prior statistic for this gene (from the cohort manifest). |
| `score` | The calibrated model probability -- the credibility score. `NaN` for genes excluded from scoring (e.g. `extra_flagged_genes`). |
| `score_percentile` | Percentile rank of `score` among scored genes (0-100). |
| `tier` | The credibility tier assigned from `score` alone. |
| `verdict` | A coarser credibility verdict derived from `tier` (the console summary counts `ROBUST` vs `INTERMEDIATE`). |
| `flag` | Whether the audit flagged this gene for review. Advisory: never changes `score`, `tier` or rank. |
| `flag_reason` | Why `flag` is set (empty when it is not) -- e.g. an architecture/QC problem such as too few case-carrying variants, when the cohort manifest declares the columns the audit needs. |
| `y` | `1` if the gene is in the label set, else `0`. |

### Per-axis ablation (console)

| Column | Meaning |
|---|---|
| `family` | The axis name, or the baseline axis (always `d_lo`/`d_md`/`d_hi` `== 0.0` there, by construction). |
| `auc` | Out-of-fold AUC of baseline+axis (uncalibrated). |
| `d_lo`, `d_md`, `d_hi` | 2.5th/50th/97.5th percentile of a paired, gene-resampling bootstrap on the delta-AUC of this axis over the baseline alone. |
| `d_lo_env`, `d_hi_env` | *Only when `--seed-sweep > 1`.* The envelope described above: the union of `[d_lo, d_hi]` and the across-seed delta range, rounded to 3 decimal places. Never narrower than `[d_lo, d_hi]`; not guaranteed to be strictly wider (depends on the data). |
| `n_seeds` | *Only when `--seed-sweep > 1`.* The sweep width actually used -- the same count as the CLI's `--seed-sweep`, `Config.seed_sweep` and `Evaluator.evaluate`'s `n_seeds` parameter. |

### Null summary (`--null-out`, or the console's abbreviated view)

| Column | Meaning |
|---|---|
| `axis` | Candidate axis name. |
| `n_perm` | Permutations actually present in this null. |
| `n_candidates` | Number of axes searched -- the width the selected-maximum null was built over. |
| `null_mean`, `null_sd`, `null_p95` | Mean / SD / 95th percentile of *this axis's own* per-axis permutation-delta distribution. |
| `selmax_median`, `selmax_p95` | Median / 95th percentile of the selected-maximum null -- identical across every axis's row, since it is one null shared by all of them. |
| `observed` | The real-label delta-AUC for this axis, from the null's own scorer -- not the ablation's `d_md`. |
| `p_per_axis` | P-value of `observed` against this axis's own per-axis null. Correct only if the axis was pre-specified. |
| `p_selected_max` | The multiplicity-corrected p-value of `observed` against the selected-maximum null -- the number that answers "could the best of the N offered axes have arisen by chance?" |

```text
      axis  n_perm  n_candidates  null_mean  null_sd  null_p95  selmax_median  selmax_p95  observed  p_per_axis  p_selected_max
expression     200             1     0.0015    0.021     0.039          0.041       0.079     0.061       0.065            0.07
```

### Blocks summary (console only)

| Field | Meaning |
|---|---|
| Block count | How many paralogue blocks the universe was split into. |
| Largest block, and its share | Size and fraction of the universe held by the single largest block -- compared against `--max-block-frac`. |
| Genes in a block of size > 1 | How many genes actually share a block with at least one paralogue (the rest are unconstrained singletons). |
| Genes matched | How many of the universe's genes were found in the `--blocks` table at all -- the signal for a silent key mismatch. |

There is no file output for the blocking report; it is printed to the console only. From the Python API it is `RerankResult.blocks` (a `BlockReport`), available whether or not the run also used `--n-perm`.

## Chunked runs on a cluster

There is no `--chunk`/`--n-chunks` CLI flag: a CLI process that computed one chunk and wrote only a coarse summary could never be combined with its siblings, so a user fanning a permutation null out across `k` array tasks would be left with `k` p-values and no way to merge them. Chunking a null is an **API-only** feature, built for exactly that fan-out, and it pickles cleanly because `NullDistribution` (and the `ControlSetting` it carries) are plain dataclasses of NumPy arrays, tuples and other pickle-safe values.

Each array task builds the *same* `Config` used everywhere else, but with a `nulls=` that names its own slice of the permutation index:

```python
import dataclasses
import pickle

from hvantk.algorithms.rerank import Config, NullConfig, rerank

# `base_config` is the Config every task in the array shares -- everything the same
# except `nulls` -- build it exactly as you would for a single-process run, with
# `nulls=None`.
base_config: Config = ...

N_CHUNKS = 20
N_PERM = 2000
SEED = 42  # identical across every task in the array

def run_chunk(chunk_index: int) -> None:
    cfg = dataclasses.replace(
        base_config,
        nulls=NullConfig(n_perm=N_PERM, chunk=chunk_index, n_chunks=N_CHUNKS, seed=SEED),
    )
    res = rerank(cfg)
    with open(f"null_chunk_{chunk_index:03d}.pkl", "wb") as fh:
        pickle.dump(res.nulls, fh)
```

A final step loads every chunk's pickle and merges them into one `NullDistribution`:

```python
import glob
import pickle

from hvantk.algorithms.rerank import NullDistribution

parts = []
for path in sorted(glob.glob("null_chunk_*.pkl")):
    with open(path, "rb") as fh:
        parts.append(pickle.load(fh))

merged = NullDistribution.merge(parts)
summary = merged.summary()  # axis, n_perm, ..., observed, p_per_axis, p_selected_max
summary.to_csv("null_summary.tsv", sep="\t", index=False)
```

Read before relying on this:

- Every chunk must share `NullConfig.seed` and `NullConfig.n_perm` -- the base seed and the *planned* (whole-run) permutation count, not the chunk's own slice size; `merge` refuses chunks built from different ones.
- Permutation `i`'s seed is `NullConfig.seed + i`, so two base seeds closer together than `n_perm` share permutations (`seed=42` and `seed=43` share `n_perm - 1` of their draws) -- move an independent replicate's seed by at least `n_perm`.
- Each array task reruns the **whole** `rerank()` call -- scoring and the ablation, not only its slice of the null -- so this fan-out only pays for itself when the permutation null, not the scoring itself, dominates the wall-clock.
- Pickled `NullDistribution` objects are tied to the `hvantk` version that wrote them; do not merge pickles written by different versions.
- `merge` raises rather than silently combining when any chunk: offered different candidate axes; was built under a different control setting (different leakage or selection policy, baseline, fold count, or **block digest** -- every chunk must have used the same `--blocks` table); carries `observed` deltas that disagree with the others' (a free check that every chunk ran on the same data, seed and estimator); or overlaps another chunk's permutation indices.
