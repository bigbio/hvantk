# ReRank: Gene Re-ranking by Multi-Omic Credibility

`hvantk rerank` re-ranks a gene burden or prior statistic against one or more gene-keyed feature axes -- constraint, expression, or any other numeric table keyed by gene -- by fitting a calibrated, out-of-fold gradient-boosted model. It assigns credibility tiers and flags cohort artifacts for review (advisory only: the audit never overrides the ranking). The feature tables and the cohort's prior statistic are supplied by the caller through a declarative YAML config; the engine itself is disease-agnostic.

## Overview

A per-axis result -- "axis X adds +0.03 AUC over the constraint baseline" -- is only as informative as the questions that have been asked of it. Three questions recur, and each is answered by its own module rather than by inspection of the output table alone:

| Risk | Question | Answered by |
|---|---|---|
| Circularity | Was this feature trained on the same evidence the label was derived from? | `provenance.py` (`Config.feature_provenance` / `label_provenance`, `rerank_arms`) |
| Presence leakage | Does a column's *missingness*, not its value, predict the label? | `leakage.py` (`Config.leakage`, `LeakagePolicy`) |
| Multiplicity | Several axes were offered and the best one reported -- could that delta have arisen by chance? | `nulls.py` (`Config.nulls`, `NullConfig`, `--n-perm`) |

The **Controls** section below documents the multiplicity control on the right of that last row, plus two further controls that ship alongside it and ask a related but distinct question -- not *is the delta real*, but *is the delta an artefact of one arbitrary modelling choice*:

- **Paralogue-blocked cross-validation** -- does the delta survive when a gene family is not allowed to straddle a cross-validation fold?
- **Multi-seed evaluation** -- does the delta survive a change of the one cross-validation partition every earlier run held fixed?

All three are opt-in CLI flags, and every one of them defaults to today's behaviour: a bare `hvantk rerank -c config.yaml -o out.tsv` is unaffected by any of the flags below unless you pass them. Circularity and presence leakage -- the other two rows of the table above -- are documented further down, in **Python-only controls**: they are `Config` fields with no CLI flag yet.

## Quick Start

```bash
hvantk rerank -c config.yaml -o out.tsv
```

This is the minimal invocation: a config naming a cohort manifest, one or more feature tables, and a label gene list (see **The config YAML** below). It writes a ranked TSV to `out.tsv` and prints a short console summary. This excerpt comes from one small synthetic cohort built purely for illustration -- a couple hundred genes, three feature axes, no real biology or gene names. The blocked walk-through further down (**Multiplicity**, **Blocked CV**, **Seeds**) uses a second, separate synthetic fixture of its own -- same shape, plus a paralogue-block structure -- introduced where it first appears, so its numbers agree with each other section to section without needing to match this unblocked one:

```text
Wrote 200 genes -> out.tsv
Credibility: ROBUST=80 / INTERMEDIATE=120 (over 200 scored genes)
Audit: 0 gene(s) FLAGGED for review (advisory; ranking not overridden). Reasons: {}
Per-axis ablation (delta-AUC over constraint):
      family   auc   d_lo   d_md   d_hi
  constraint 0.706  0.000  0.000  0.000
  expression 0.738 -0.020  0.032  0.082
conservation 0.658 -0.089 -0.050 -0.009
```

The ablation table is the per-axis result: out-of-fold AUC of baseline+axis, and a paired gene-resampling bootstrap interval (`d_lo`/`d_md`/`d_hi`) on the delta over the baseline alone -- here `expression` adds a modest, not-quite-conclusive gain (the interval straddles zero) while `conservation` makes things worse (a genuine illustration of an axis that should not be added, not a hypothetical one). See **Output columns** for the full column reference, including the columns the controls below add.

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
| `--max-block-frac FLOAT` | Blocked CV | `0.1` |
| `--seed INTEGER` | Seeds | `42` |
| `--seed-sweep INTEGER` | Seeds | `1` |

### Multiplicity -- `--n-perm`, `--null-out`

Offering several feature axes and reporting the best one's delta-AUC is itself a selection procedure, and the ablation table's bootstrap interval has nothing to say about it -- it describes one axis in isolation. `--n-perm` builds two permutation nulls that do:

- **Per-axis** -- permute the label, refit the baseline and baseline+*this* axis, record the delta. Correct when the axis was pre-specified.
- **Selected-maximum** -- permute, refit, add *each* offered axis, record the *largest* delta. Correct when the axis is reported *because* it came top over the other candidates. Only this one is a multiplicity correction: the maximum over several candidate axes is stochastically larger than any single one of them, so a per-axis null can sit near zero while the selected-maximum null does not.

Both nulls refit the baseline on every permutation -- scoring a permuted-label model against the unpermuted baseline would measure the permutation, not the axis. P-values use `p = (1 + #{null >= observed}) / (1 + n_perm)` (Phipson & Smyth 2010): a naive `(null >= observed).mean()` reports `p = 0` for a value no permutation reached, which a finite permutation set can never actually license. `p_selected_max` is the single-step max-T (Westfall-Young) adjusted p-value for that axis's observed delta. The observed deltas themselves come from the null's own scorer, unrounded -- **not** the ablation table's `d_md` (a bootstrap *median*, rounded to 3 decimal places: a different statistic, computed a different way). They are, however, the SAME quantity as the ablation table's own `auc` column difference for that axis: the null's scorer and `_raw_oof` run the identical columns, seed, folds and blocks, so `observed` and `auc[axis] - auc[baseline]` agree up to `auc`'s 3-decimal rounding. In this run, `observed = 0.032381` for `expression` below, against `0.807 - 0.775 = 0.032` from the **Seeds** table further down -- checkable against each other precisely because both come from the one run this whole blocked walk-through uses.

> **What null hypothesis this tests.** Permuting the label tests the *global* null that no feature, baseline included, carries label information -- not the *conditional* null that this one axis adds nothing given the baseline already in the model. With an informative baseline, the permutation null of the delta is wider than the delta's true sampling spread under the conditional null, so the test is valid but conservative: failing to clear the selected-maximum null is weak evidence that the axis adds nothing, not proof that it does not.

That conservatism assumes genes are exchangeable under the label permutation. Paralogues share sequence, constraint, expression and disease status, so when labels and features cluster by gene family that assumption breaks and the *default* (unblocked) null becomes **anti-conservative** instead -- understating its own spread. That clustering is a property of the data, not of a flag: omitting `--blocks` does not remove it, it only removes the correction for it. Pass `--blocks` whenever that clustering is doubtful; when it is set, the null permutes labels **by block** (whole blocks swapped among blocks of the same size, then shuffled within each block -- Winkler et al. 2015, *NeuroImage*), which is what blocking requires of the null too. A block whose size is unique in the universe has no same-size sibling to swap with, so it is only ever shuffled within itself -- a real power caveat for that block, not a bug.

The null is tied to one **control setting** -- the leakage and selection policies, the provenance arm, the exact baseline and candidate columns, the fold count actually used, the CV seed (`--seed`), and which blocking (if any) the folds used -- and refuses to answer about a delta computed under a different one, raising `ControlSettingMismatch` rather than silently comparing incomparable numbers. The ablation and the null both score at `Config.folds` folds (5 by default) -- the same fold count the headline score uses, kept in lockstep end to end rather than a separate, fixed ablation fold count. The null's `observed` deltas are computed at `--seed` alone (never averaged over `--seed-sweep`); see **Seeds** below for why the two are reported separately.

```bash
hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_complete_set.txt --n-perm 20 --seed-sweep 3 --null-out null_summary.tsv
```

(This example also passes `--blocks` and `--seed-sweep`, each covered in its own section below. This is the first of three commands run against ONE separate synthetic fixture built for this blocked walk-through -- 200 genes, the same three feature axes as **Quick Start**, plus a paralogue-block structure -- so that every number here, in **Blocked CV** and in **Seeds** below is drawn from that one consistent run.) It adds, after the ablation table:

```text
Multiplicity: 2 axes searched over 20 permutations; selected-maximum null median +0.0426
Setting: arm='all' leakage=None selection=None folds=5 seed=42 block_digest='7cd6dbdd71eb' baseline=['constraint'] candidates={'expression': 1, 'conservation': 1}
        axis  observed  selmax_median  p_selected_max
  expression  0.032381       0.042649        0.571429
conservation -0.028333       0.042649        0.904762
```

The `Setting:` line is `ControlSetting.describe()`: the leakage and selection policies, the provenance arm, the fold count, the CV seed, the baseline/candidate columns and -- whenever `--blocks` is set -- the digest of the blocking the folds actually used, so a p-value can always be traced back to exactly what it was computed under.

`--null-out` writes the full per-axis summary (see **Output columns**); without it, the console shows the summary line, the `Setting:` line, and the four columns above only. `--null-out` requires `--n-perm`: passing `--null-out` without `--n-perm` is rejected before any work starts (`Error: --null-out describes a permutation null; pass --n-perm N too.`), as is a negative `--n-perm` (`Error: --n-perm must be >= 0; got -5.`, for example) -- both a plain `click` failure, never a traceback. `--n-perm` also needs at least one candidate axis besides the baseline: a config naming only a baseline has nothing to search over, and the CLI rejects it before writing any output the same way (`Error: --n-perm needs at least one feature axis besides the baseline; this config has only ['constraint'].`, for example) rather than silently skipping the correction. Via the Python API, `NullConfig()` alone (with no `n_perm=`) defaults to 200 permutations; the CLI's own default is `0` (off) so that a plain `hvantk rerank` run costs nothing extra.

**Cost.** Each permutation refits the baseline and every candidate axis -- roughly `n_perm x (1 + axes)` five-fold fits, and a selection policy (RFECV especially) multiplies that further. On the order of 10-20 seconds per permutation for about 1,000-4,000 genes and 9-17 axes on a 12-core laptop, i.e. tens of minutes to around an hour of wall-clock for `n_perm = 200`. Progress is logged at `INFO`, below the CLI's default `WARNING` level -- run `hvantk -v rerank ...` to see it as it happens.

### Blocked CV -- `--blocks`, `--max-block-frac`

Every number `hvantk rerank` produces by default comes from plain stratified 5-fold cross-validation. Gene families break the independence that relies on: paralogues share sequence, constraint, expression pattern and disease status, so a family split across train and test lets the model partly recognise a relative rather than generalise. `--blocks` switches to cross-validation folds that keep each gene family whole.

```bash
hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_complete_set.txt
```

(The same blocked synthetic fixture introduced in **Multiplicity** above -- so its `Blocks:` line below, which depends only on the universe and the blocks table, is identical to the one that run also prints.)

`--blocks` takes an HGNC complete-set TSV with `symbol`/`gene_group` columns (and, optionally, `status` -- when present, only `Approved` rows are kept). Genes are matched against it by their **approved HGNC symbol only**; a gene absent from the table, or known to your universe only by an alias or a previous symbol, becomes its own unconstrained singleton block rather than an error -- unless *none* of the universe's genes match, which almost always means the table's gene keys disagree with your universe's (e.g. Ensembl IDs against HGNC symbols) and *is* an error. A match rate below 50% is logged as a `WARNING`, which reaches the console by default (`hvantk` sends `WARNING` and above to stderr unless verbosity is turned down) -- but the quickest check either way is the matched count the console's `Blocks:` line below always prints too.

Blocking uses each gene's **first-listed** HGNC gene group -- in HGNC's own list order, since HGNC does not designate a "primary" one -- rather than connected components over the (multi-membership, pipe-separated) `gene_group` field. Group membership chains transitively, so a union-find closure over shared membership can collapse a large fraction of a gene universe into a single component; blocking on it would then let `StratifiedGroupKFold` place nearly the whole universe in one fold, and the pooled out-of-fold AUC would no longer estimate the same quantity as the unblocked run. The tell-tale of that mistake is a *blocked* AUC scoring **materially** higher than the random-fold one, beyond the across-seed spread `--seed-sweep` reports -- blocking is harder only *on average*, so a single correctly blocked run can still land above the random-fold AUC by chance, and a difference within that spread is not itself suspicious. First-listed grouping keeps blocks small; the residual leak it accepts -- two genes sharing only a later-listed group -- is a far smaller error than an incomparable estimator.

**The ceiling is a hard abort, not a warning.** `--max-block-frac` (default `0.1`) caps the largest block's share of the universe; going over it fails the run instead of producing a plausible-looking number, because the failure mode is otherwise silent -- the output looks fine, and a counterintuitive AUC is the only tell, and a reader has no particular reason to question it:

```text
$ hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_dominant_family.txt
Error: the largest paralogue block (group 'Family_Huge') holds 40/200 units (20.0%), over the 0.1 ceiling (max_block_frac). StratifiedGroupKFold must place a whole block in one fold, so one fold would hold that entire block and the pooled out-of-fold AUC would not estimate the same quantity as the unblocked run. Raise max_block_frac deliberately if you accept that, or restrict the universe.
```

No output file is written when this happens. `--max-block-frac` only means anything alongside `--blocks`; passing it alone is rejected before any work starts (`Error: --max-block-frac only applies to paralogue-blocked folds; pass --blocks too (or drop it).`). Raise it deliberately (`--max-block-frac 0.5`, say) if a dominant block is expected and accepted for your universe.

Two further mechanics worth knowing about: the grouped splitter is a **correctly stratified** one, not scikit-learn's `StratifiedGroupKFold(shuffle=True)` -- as of scikit-learn 1.7.2 that shuffles the per-group class-count rows before balancing folds but then assigns each group to a fold under its *original*, unshuffled index, so groups stay whole but the folds stop actually being stratified; `hvantk` instead relabels group IDs with a seeded random permutation and calls the unshuffled (correctly balancing) path underneath, which randomises the same tie-breaks the upstream shuffle intended without inheriting its defect, and stays correct even if that upstream behaviour is later fixed. Separately, a class confined to too few blocks to cover the requested fold count is refused with a clear error rather than silently producing a fold with no positive (or no negative) training example, which would otherwise surface only as an opaque scikit-learn warning, or several calls later as "found array with 0 sample(s)".

Blocking is threaded everywhere a fold boundary matters: the outer cross-validation of both the headline score and the ablation, the multi-seed sweep, and -- as described above -- the permutation null. There are three documented exceptions to this blocking: `CalibratedClassifierCV`'s own inner calibration split, which takes only an integer fold count and so cannot see the blocks; a family can therefore straddle *that* split within one already-blocked outer training fold. Similarly, with a selection policy that uses the RFECV wrapper (`SelectionPolicy(wrapper="rfecv")`), RFECV's inner `StratifiedKFold` inside each training fold is likewise unblocked (it is opt-in). The third applies only when a fit sees more than 10,000 rows: the estimator leaves `HistGradientBoostingClassifier`'s `early_stopping="auto"`, which scikit-learn (1.7.2) switches on for `n_samples > 10,000`, and each fit then holds out a class-stratified random 10% validation split drawn with the run's seed and blind to the blocks -- a family can straddle *that* split too, and `--seed-sweep` then also varies which rows every fit holds out, on top of the fold assignment. The estimator is left as is deliberately (changing it would change every historic result); below that size no such split exists. The per-axis deltas that are the scientific claim never go through the first two -- they use the uncalibrated, fully-blocked scorer in `evaluator.py` (fully blocked below 10,000 rows; the early-stopping split applies to it exactly as to the headline). Custom API scorers passed directly to `permutation_deltas` that carry no `.groups` attribute are not checked against `blocks=` at all -- only scorers built by the library's own `oof_scorer` declare one.

The per-axis bootstrap interval (`d_lo`/`d_md`/`d_hi`) is blocked the same way, and for the same reason. With `--blocks` set, `_boot_ci` resamples whole paralogue blocks instead of genes: each draw picks as many block ids as there are blocks, uniformly and with replacement, and takes every member gene of the ids it picked (a cluster bootstrap) -- a size-6 family rides in or out of a draw as one unit, exactly as it rides into one fold or the other. Resampling genes independently here while the folds are blocked would silently understate the interval by exactly the amount blocking exists to correct for, which is why this is not a separate flag: it follows `--blocks` automatically. Without `--blocks`, the interval is unchanged -- genes resampled independently, bit-identical to every run before this.

The console's `Blocks:` line (shown whenever `--blocks` did not abort) reports the blocking, including the matched-gene count that catches a silent key mismatch:

```text
Blocks: 89 paralogue block(s); largest 4 (2.0%); 160 gene(s) in a block of size > 1; 200/200 genes matched the gene-group table; digest 7cd6dbdd71eb
```

### Seeds -- `--seed`, `--seed-sweep`

`--seed` (default `42`, the historic value) is the one seed that drives the cross-validation partition (both the headline score and the ablation), the gradient-boosted estimator's own randomness, the bootstrap resample, and the permutation null's base seed. Changing it produces a different, equally valid partition -- it is not a tuning knob, and the default is unchanged so every result produced before `--seed` existed still reproduces exactly.

Doubling as the null's base seed has one consequence worth knowing before you rely on it: permutation `i` is drawn with seed `--seed + i` (see **Chunked runs** below), so bumping `--seed` by 1 is not an independent replicate of the null -- it reuses `n_perm - 1` of the previous run's permutations (`--seed 42` and `--seed 43` at `--n-perm 20` share 19 of their 20 draws). An independent replicate needs a base seed at least `n_perm` away from the first.

```bash
hvantk rerank -c config.yaml -o out.tsv --blocks hgnc_complete_set.txt --seed-sweep 3
```

(The same blocked synthetic fixture as **Multiplicity** and **Blocked CV** above -- shown
with `--blocks` because the envelope and the block bootstrap are both best illustrated
together on blocked folds; the envelope works the same way without blocking too.)

The bootstrap interval (`d_lo`/`d_md`/`d_hi`) asks how a delta would move if the *sample* moved -- genes resampled independently when unblocked, or whole paralogue blocks resampled together when `--blocks` is set (see **Blocked CV** above). Either way it cannot see a second variance component: which genes (or blocks) land in which cross-validation fold. `--seed-sweep N` recomputes each axis's delta under `N` consecutive cross-validation seeds (`--seed`, `--seed` + 1, ...) and reports the **envelope** -- the union of the bootstrap interval and the across-seed range -- as `d_lo_env`/`d_hi_env`, beside the untouched `d_lo`/`d_md`/`d_hi`. The same sweep width travels under a different name at each layer: the CLI's `--seed-sweep`, `Config.seed_sweep`, `Evaluator.evaluate`'s `n_seeds` parameter, and the ablation table's own `n_seeds` column all refer to the identical count -- so a reader who encounters more than one of these names is not looking at two different settings:

```text
Per-axis ablation (delta-AUC over constraint):
      family   auc   d_lo   d_md  d_hi  d_lo_env  d_hi_env  n_seeds
  constraint 0.775  0.000  0.000 0.000     0.000     0.000        3
  expression 0.807 -0.011  0.033 0.088    -0.011     0.088        3
conservation 0.747 -0.064 -0.027 0.011    -0.064     0.011        3
```

The envelope is a **union, not a calibrated interval** over both components: the two are not independent draws from one distribution, and combining them properly would need a nested design nobody has run. `d_lo_env`/`d_hi_env` are that union rounded to 3 decimal places, exactly like `d_lo`/`d_hi`. The union is guaranteed never to be narrower than the bootstrap interval alone, but **not** guaranteed to be strictly wider -- that depends on the data: when the across-seed range happens to fall entirely inside `[d_lo, d_hi]`, `d_lo_env == d_lo` and `d_hi_env == d_hi` and no widening is visible, which is itself informative (it says the CV partition was not, on this data, a meaningful source of extra spread), not a sign that the sweep did nothing. It exists so a reader cannot mistake one lucky (or unlucky) partition for a measured effect -- a single cross-validation partition can land anywhere within that spread. With the default `--seed-sweep 1`, the three extra columns are **absent** from the ablation table entirely, so a plain `hvantk rerank` run is byte-for-byte what it always was.

The permutation null (`--n-perm`, above) already recomputes the stratified partition for every permutation, so it already contains partition noise as one of its own components. The envelope is a separate report on the *observed* delta and is never folded into the null on top of that -- doing so would count partition variance twice. For the same reason, the null's `observed` deltas are always computed at `--seed` alone, never averaged over the swept seeds.

`SelectionPolicy.seed` (used only when a feature-selection policy is configured) is independent of `--seed` by design -- it seeds the selection wrapper's own inner cross-validation, not the outer partition `--seed` controls. Turning on `--seed-sweep N` multiplies the ablation's per-axis cost by roughly `N` (one refit per axis per seed); the baseline is fit once per seed and shared across every axis rather than refit per axis, so the extra cost scales with the number of candidate axes times `N`, not more.

## Python-only controls: selection, circularity and presence leakage

Three further controls exist only on `Config`: `hvantk rerank` does not yet take a
`selection:`/`leakage:` key in `config.yaml`, nor any matching CLI flag. `Config` is a
frozen dataclass, so set them with `dataclasses.replace(config, field=value)` rather than
assigning to an existing instance.

### Feature selection -- `Config.selection`

`Config.selection: Optional[SelectionPolicy]` (default `None`, meaning off: every column of
every axis is used, exactly as before this existed). When set, `engine.rerank` filters each
axis down to an informative, non-redundant subset -- **within that axis only**, independently
of every other axis -- re-run inside every cross-validation fold on the training slice alone.
Within-axis is what keeps a per-axis delta-AUC meaningful: it compares each axis's best subset
against the baseline's best subset, so a wide axis is not penalised merely for carrying
redundant columns; fold-nested is what stops the reported delta-AUC from being inflated by a
selection that has already seen the held-out labels.

Three steps, in this order (`hvantk/algorithms/rerank/selection.py`):

| `SelectionPolicy` field | Default | Meaning |
|---|---|---|
| `univariate` | `"auc"` | `"auc"` scores each column by rank-based \|AUC - 0.5\| and keeps those surviving within-axis BH-FDR at `q`; `"none"` skips this step. |
| `q` | `0.10` | BH-FDR level for the univariate filter, within axis. |
| `redundancy` | `"spearman"` | Greedy, best-first: drop a column whose \|Spearman rho\| against an already-kept column reaches `redundancy_max`. `"none"` skips this step. |
| `redundancy_max` | `0.75` | The Spearman-rho cutoff above. |
| `wrapper` | `"none"` | `"none"` runs only the two filters above. `"rfecv"` adds RFECV (RandomForest-backed) on what survives them. |
| `wrapper_estimator` | `"random_forest"` | The only value the schema allows; not a knob anything else reads. |
| `inner_folds` | `3` | RFECV's own inner `StratifiedKFold` count. |
| `seed` | `42` | Independent of `Config.seed`: seeds the wrapper's inner CV only, never the outer partition. |

`wrapper="rfecv"` is off by default on measured evidence, not on principle: across the
cohorts it was tried on, its eliminations concentrated where positives were fewest (the
signature of a wrapper fitting inner-CV noise), and it can prune columns from the *ablation
baseline* axis itself, which silently inflates every other axis's delta-AUC. Turn it on
deliberately, for a genuinely wide axis (dbNSFP's ~45 correlated predictors is the motivating
case) with enough positives to trust an inner CV.

A YAML policy validates against `hvantk/resources/schemas/selection_policy.schema.json`
(every field optional; an empty document yields the defaults above) and loads with
`hvantk.algorithms.rerank.selection.load_policy(path) -> (SelectionPolicy, equivalence)`. The
schema also accepts an `equivalence` block -- source-name classes for the circularity check
below -- which `load_policy` pops off and returns separately, since the resolver consumes it
and the filters never see it; omit it to use `provenance.DEFAULT_EQUIVALENCE`.

### Circularity -- `Config.feature_provenance`, `label_provenance`

Statistical selection cannot detect circularity; it rewards it -- a univariate filter ranks a
predictor highly *because* it was trained on genes like the ones in the label set. Circularity
is declared, not computed from the data: `Config.feature_provenance` maps each column to a
`frozenset` of the sources it was trained on (an empty `frozenset()` means "trained on nothing
label-derived"; `None` for the whole dict means every declaration is unavailable), and
`Config.label_provenance` is the `frozenset` of sources the labels themselves derive from (an
explicit `frozenset()` asserts "derived from nothing curated"). Leaving `label_provenance`
`None` while `feature_provenance` is set makes `rerank_arms` raise `ValueError` rather than
run: an undeclared label source would otherwise conflict with nothing, so every circular
predictor would be admitted and the run would look healthy. Source names come from a plugin's
own manifest (`scores: {trained_on: [...]}`, see the
[plugin conventions](https://github.com/bigbio/hvantk/blob/main/hvantk/skills/_conventions/SKILL.md),
section 6) paired with whatever curated database the labels were derived from.

`hvantk.algorithms.rerank.engine.rerank_arms(config)` -- not `rerank()` -- runs the
circularity split, and only when both `config.selection` and `config.feature_provenance` are
set (and, per the paragraph above, `config.label_provenance` is not left `None`); otherwise it
returns `{"all": rerank(config)}` -- a single, unrestricted arm, unchanged from before this
existed. When it does split, it returns `{"clean": RerankResult, "all": RerankResult}`:
`clean` keeps only columns that do not conflict with the label's sources and is **always the
headline**; `all` additionally admits conflicted and undeclared columns, existing only to
measure the circularity channel as a number (`arms["all"].selection.n_conflicted`,
`.n_unknown`) rather than leave it assumed away.

### Presence leakage -- `Config.leakage`

Circularity asks what a predictor was *trained* on; this asks which genes it was *run* on.
`Config.leakage: Optional[LeakagePolicy]` (default `None`, off) bars a column whose
**missingness**, not its value, predicts the label -- the motivating case is a predictor
computed only for a curated, well-studied (and so disease-enriched) protein subset, where the
bare flag "was this score computed for this gene" already carries label information before
the score's own value is read. `LeakagePolicy.q` (default `0.10`) is the BH-FDR level for the
presence/label association; `LeakagePolicy.min_auc` (default `0.55`) is the presence-AUC floor
an association must also clear. Both conditions are required: significance alone flags a
trivial coverage difference once a cohort is large, and effect size alone flags noise in a
small one. Independent of `selection`: the leakage check runs in front of within-axis
selection whenever both are configured, since a column whose missingness tracks the label
would otherwise look, to a purely statistical filter, like a column worth keeping.

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
| `d_lo`, `d_md`, `d_hi` | 2.5th/50th/97.5th percentile of a paired bootstrap on the delta-AUC of this axis over the baseline alone: genes resampled independently when unblocked, or whole paralogue blocks resampled together (a cluster bootstrap) whenever `--blocks` is set. |
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
  expression      20             2   0.007976 0.050086  0.078080       0.042649    0.105902  0.032381    0.333333        0.571429
conservation      20             2   0.001140 0.055931  0.092714       0.042649    0.105902 -0.028333    0.761905        0.904762
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
- `merge` raises rather than silently combining when any chunk: offered different candidate axes; was built under a different control setting (different leakage or selection policy, baseline, fold count, CV seed, or **block digest** -- every chunk must have used the same `--blocks` table); carries `observed` deltas that disagree with the others' (a free check that every chunk ran on the same data, seed and estimator); or overlaps another chunk's permutation indices. It also raises when the chunks together leave a gap in the planned permutations (a chunk that never arrived, such as a pre-empted array task), naming the missing indices; pass `merge(parts, allow_partial=True)` to accept the smaller null deliberately.

## See also

- [External cohorts](cohort.md) -- write the cohort manifest `rerank` reads for the prior
  statistic and the audit columns.
- [Re-rank example](../examples/rerank.md) -- a worked walkthrough on small synthetic
  fixtures, with the commands run and their actual output.
