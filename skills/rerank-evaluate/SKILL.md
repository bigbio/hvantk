---
name: rerank-evaluate
description: 'Use when ranking a cohort''s genes by multi-omic evidence with hvantk rerank and evaluating the ranking against a label set you choose: putting the cohort table, the feature axes and the labels into one gene-symbol space, choosing the feature axes (gene structure, constraint, expression), building a label set from the ClinGen or GenCC datasets, running with paralogue blocks and the case/control audit, and computing ROC and precision-recall metrics. Also when rerank stops with "label/feature join coverage too low", "too few examples for 5-fold", an unknown config key, or a blocks table that matches no gene.'
user-invocable: true
argument-hint: "<cohort.yaml> [--features name=path ...] [--labels clingen|gencc|file] [--blocks hgnc_complete_set.txt]"
metadata:
  stage: analysis
  status: provisional
  scale: node
---

# Rank and evaluate a cohort's genes with hvantk rerank

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract. Jobs run through `hpc-hail-job`.

## 1. Scope and outcome

A rerank config YAML, the ranked gene table `hvantk rerank` writes, and an evaluation
receipt (ROC AUC, average precision, label counts, tier and flag counts) over it. The
cohort's own statistic (its `prior`) is carried along, never modelled: `hvantk rerank`
ranks genes by the feature axes only, and the prior rides along as `prior_stat`. The
universe -- which genes get ranked -- is the union of the feature axis tables, so
restricting what is ranked means restricting the feature tables, not the cohort or the
labels.

Use it once a cohort has a gene-level table (`cohort-attach-burden`) and a gene index
(`gene-index`) to build feature axes from. Not this skill: producing the cohort table
itself, building the datasets a feature axis is drawn from (`dataset-build-at-scale`),
or building the gene index (`gene-index`).

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value. Every later command that reads a profile value starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and image as in `hpc-hail-job` (commit, digest, run notes).
3. One identifier space for the cohort table, every feature axis and the labels:
   `hvantk rerank` does no identifier mapping at all -- a cohort keyed by Ensembl
   `gene_id`, a feature axis keyed by HGNC symbol and a label file of symbols never
   join. Decide the space now (symbols, if `--blocks` is used -- it matches HGNC
   approved symbols only) and map everything into it before writing the config.
   `references/inputs-and-identifiers.md`.
4. The label set's definition: source (ClinGen, GenCC, or a file you already have),
   minimum classification, and which groups -- decided and written down now, not
   improvised mid-run. A label set built from a curated disease database scored
   against that same database's own labels is circular; `references/feature-axes.md`
   and `references/audit-and-leakage.md` say where that line is.
5. The universe: which genes will actually be ranked. It is whatever the feature
   tables cover (their outer join), so decide it by choosing what goes into the
   feature tables -- restricting the cohort table or the label file does not shrink
   the universe by itself.
6. The feature axes: at least one, each a gene-keyed table with no identifier
   collisions across axes (`references/feature-axes.md`); chosen for no circularity
   against the label source (`references/audit-and-leakage.md`).

## 3. Procedure

1. **Labels**, as a job each: build the gene-disease table, extract a gene-set
   collection from it, then flatten to a flat symbol list.
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-clingen --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess clingen:gene-disease --raw-dir "$HVANTK_WORK/raw/clingen" \
       --output "$HVANTK_WORK/<run>/clingen.ht" --skip-download --no-check-drift
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-geneset --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     hvantk genesets clingen --ht "$HVANTK_WORK/<run>/clingen.ht" \
       --min-classification Moderate -o "$HVANTK_WORK/<run>/labels.json"
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-labels --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" labels \
       --geneset "$HVANTK_WORK/<run>/labels.json" \
       --restrict-to <the universe's own symbol table> --restrict-column symbol \
       --out "$HVANTK_WORK/<run>/labels.txt"
   ```
   `--min-classification` is "this level or stronger" (`Moderate` means Moderate,
   Strong or Definitive); omit `--group` to union every ClinGen GCEP group, i.e. any
   disease, or pass one or more `--group` names to narrow it. `genesets gencc` is the
   GenCC equivalent (`--group-by submitter`). Restrict labels to the universe now, not
   by lowering `min_label_coverage` later (Traps).
2. **The cohort table, mapped into symbol space**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-cohort --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" cohort \
       --table <gene_burden.tsv> --key-column gene --gene-index <gene-index.ht> \
       --hgnc-tsv <hgnc_complete_set.txt> --out "$HVANTK_WORK/<run>/cohort_symbols.tsv" \
       --manifest-in <cohort.yaml> --manifest-out "$HVANTK_WORK/<run>/cohort_symbols.yaml"
   ```
   This adds a `symbol` column (via the gene index's own `hgnc_id`, mapped to an
   approved symbol through the HGNC TSV), drops rows with no approved symbol or whose
   symbol is shared by more than one gene id, and reports both counts. With
   `--manifest-in`/`--manifest-out` it also writes a rerank-ready manifest pointing at
   the new table, carrying over the original's `prior` and `cohort_axes` verbatim and
   setting `key`/`key_column` to `symbol`.
3. **Feature axes**, one TSV per axis, each as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-structure --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" features \
       --ht <gene-index.ht> --fields cds_length,n_coding_exons,n_transcripts \
       --hgnc-tsv <hgnc_complete_set.txt> \
       --restrict-to "$HVANTK_WORK/<run>/cohort_symbols.tsv" --restrict-column symbol \
       --out "$HVANTK_WORK/<run>/structure.tsv"
   ```
   `--ht` is any gene-keyed hvantk table (`references/feature-axes.md` for constraint
   and expression axes); `--restrict-to` is what makes the universe equal the cohort's
   tested genes here -- drop it to rank the whole gene index instead. Every axis's
   feature columns must have globally unique names (`references/feature-axes.md`).
4. **The config.** Exactly `name`, `cohort`, `features` (a list of `{name, path}`),
   `labels` (`{path}`), and optionally `min_label_coverage` -- no `prior:` key, no
   other key. `cohort` points at the symbol-keyed manifest from step 2; `features`
   lists each axis TSV from step 3 under the name used in `--out`'s basename.
5. **Run**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-rerank --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk rerank -c "$HVANTK_WORK/<run>/config.yaml" -o "$HVANTK_WORK/<run>/ranked.tsv" \
       --blocks <hgnc_complete_set.txt> --seed 42
   ```
   `--blocks` switches cross-validation to paralogue-blocked folds; it is required by
   this skill's own measurement, not optional. A cohort whose declared axes include
   `n_case_var`, `conc` and `driver_af` (the architecture axis step 2 carries over)
   auto-wires the advisory `CaseControlArchitectureAudit` -- there is no `audit:`
   config key.
6. **Evaluate**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-eval --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/evaluate-ranking.py" \
       --ranked "$HVANTK_WORK/<run>/ranked.tsv" --receipt "$HVANTK_WORK/<run>/eval.json"
   ```
   Record its `auc:` line in the run notes.

## 4. Scale-up

A genome-wide universe (about 20,000 genes, several axes) still fits one node; a
feature axis is one small table per axis, never itself a Hail job once written to
TSV/parquet. Cost grows with `--n-perm` (permutation nulls: roughly
`n_perm * (1 + axes)` five-fold fits; needs at least two feature axes, and `--blocks`
to keep the null itself from being anti-conservative) and `--seed-sweep`; run nulls as
their own job, not folded into the main run's budget. Above 10,000 rows, early
stopping opens a further seeded split inside each fold -- no action needed, just more
wall time.

## 5. Traps

- `label/feature join coverage too low: only a/b (p%) ...`. Cause: the labels are in
  another identifier space than the features, or name genes outside the ranked
  universe. Fix: restrict the label file to the universe (step 1's `--restrict-to`);
  never lower `min_label_coverage` to get past it -- that hides the mismatch rather
  than fixing it.
- A cohort keyed by Ensembl `gene_id` scored against symbol-keyed features (or the
  reverse). Cause: `hvantk rerank` does no identifier mapping at all. Fix: map the
  cohort table to the feature axes' space first (step 2), before writing the config.
- `Blocks: ... 0/... genes matched` or a `DominantBlockError`/`ValueError` from
  `--blocks`. Cause: the ranked units are not HGNC approved symbols (`blocks.py`
  matches the approved `symbol` column only -- not aliases, not Ensembl ids). Fix:
  map to symbols before running; `references/inputs-and-identifiers.md`.
- `too few examples for 5-fold OOF: ...`. Cause: fewer than `folds` (5) positive or
  negative units in the universe. Fix: a broader label set, a bigger universe, or
  fewer folds -- not a smaller `min_label_coverage`.
- The audit silently not running (no `Audit:` line in the log). Cause: the cohort
  manifest's declared axis columns do not cover all of `n_case_var`, `conc` and
  `driver_af` -- a near-miss axis (e.g. one named `architecture` with only two of the
  three) falls back to `NoAudit` with a logged warning, not an error. Fix: carry the
  original cohort's architecture axis through step 2 unchanged
  (`references/audit-and-leakage.md`).
- Duplicate symbols after mapping. Cause: two gene ids map to the same approved
  symbol. Fix: `prepare-rerank-inputs.py` drops both sides of a shared symbol and
  reports the count; the cohort loader refuses a table with more than one row per
  gene key outright (`duplicate row(s) for gene key(s) ...`) if any slipped through.
- `--n-perm` with one feature axis: rejected before any work starts (the CLI checks
  `len(features) >= 2`) -- the first axis is always the ablation baseline, so a
  one-axis config has nothing left to search a null over.
- `--max-block-frac` without `--blocks`: rejected outright (`--max-block-frac only
  applies to paralogue-blocked folds`).
- A leftover `prior:` key in the config. Cause: a config written before the
  cohort-manifest migration. Fix: delete it -- the prior now lives entirely in the
  cohort manifest's own `prior:` block, referenced via `cohort:`.
- The console's overall credibility numbers not matching a column of the per-axis
  ablation table. Cause: the ablation's per-axis AUCs are UNCALIBRATED subset fits
  (baseline alone vs. baseline+axis), a different quantity from the full model's
  calibrated, out-of-fold `score` column. Fix: compute the overall AUC from the ranked
  table (`roc_auc_score(y, score)`), never from the ablation table.
- A circular axis quietly inflating AUC (built from the same curated database the
  labels come from). `references/audit-and-leakage.md`.

## 6. Verification

1. `config.yaml` parses with exactly `name`, `cohort`, `features`, `labels` and
   optionally `min_label_coverage`; `features` lists at least one axis TSV that
   exists.
2. `ranked.tsv` carries exactly `gene, prior_stat, score, score_percentile, tier,
   verdict, flag, flag_reason, y`; every label-positive gene restricted into the
   universe is a `y = 1` row; `roc_auc_score(y, score)` reproduces the run's own
   console AUC (there is no printed overall-AUC line -- compute it from the table).
3. The console's `Audit: N gene(s) FLAGGED ... Reasons: {...}` matches the table's own
   `flag`/`flag_reason` counts; `Blocks: ... digest <12 hex>` names how many genes
   matched the gene-group table (most of the universe, for an HGNC-approved-symbol
   universe).
4. `evaluate-ranking.py`'s receipt: `n_genes`, `n_positive`, `auc`, `average_precision`,
   the positive count in the top decile, tier and verdict counts, flag counts by
   reason, and `ranked_sha256`.
5. Run notes carry `commit:`, `container:`, one `job:` line per job submitted, and
   `auc: <value>` from step 6's receipt.
