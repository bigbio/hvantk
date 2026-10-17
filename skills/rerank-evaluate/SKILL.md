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

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract. Jobs run through `hpc-hail-job`.

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

1. Site profile:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and its commit: `git -C "$HVANTK_REPO" rev-parse HEAD` and
   `git -C "$HVANTK_REPO" status --porcelain`; compare with the commit the user expects.
3. Image and its digest: `readlink -f "$HVANTK_SIF"` and
   `sha256sum "$(readlink -f "$HVANTK_SIF")"`.
4. Inputs and outputs: the cohort table, the gene index and the raw ClinGen/HGNC files are
   readable, and the output location is writable, from compute nodes -- under
   `$HVANTK_WORK` or similar, never home.
5. One identifier space for the cohort table, every feature axis and the labels:
   `hvantk rerank` does no identifier mapping at all -- a cohort keyed by Ensembl
   `gene_id`, a feature axis keyed by HGNC symbol and a label file of symbols never join.
   Decide the space now (symbols, if `--blocks` is used -- it matches HGNC approved
   symbols only) and map everything into it before writing the config.
   `references/inputs-and-identifiers.md`.
6. The label set's definition: source (ClinGen, GenCC, or a file you already have),
   minimum classification, and which groups -- decided and written down now, not
   improvised mid-run. A label set built from a curated disease database scored against
   that same database's own labels is circular; `references/feature-axes.md` and
   `references/audit-and-leakage.md` say where that line is.
7. The universe: which genes will actually be ranked. It is whatever the feature tables
   cover (their outer join), so decide it by choosing what goes into the feature tables --
   restricting the cohort table or the label file does not shrink the universe by itself.
8. The feature axes: at least one, each a gene-keyed table with no identifier collisions
   across axes (`references/feature-axes.md`); chosen for no circularity against the label
   source (`references/audit-and-leakage.md`).

## 3. Procedure

1. **Download the raw files, on the login node, through the image** (compute nodes
   have no egress):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk download hgnc --output "$HVANTK_WORK/raw/hgnc/hgnc_complete_set.txt"
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk download clingen --output-dir "$HVANTK_WORK/raw/clingen"
   ```
   `download hgnc` exits 1 if `--output` already exists (`Use --overwrite to replace
   existing file`) -- skip it if the file is already there, rather than retrying;
   `download clingen` writes a dated file per run and does not refuse an existing
   directory.
2. **The cohort table, mapped into symbol space**, as a job -- this decides the
   universe, so it comes before labels and before the feature axes that will be
   measured against it:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-cohort --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" cohort \
       --table <gene_burden.tsv> --key-column gene --gene-index <gene-index.ht> \
       --hgnc-tsv <hgnc_complete_set.txt> --out "$HVANTK_WORK/out/rerank/cohort_symbols.tsv" \
       --manifest-in <cohort.yaml> --manifest-out "$HVANTK_WORK/out/rerank/cohort_symbols.yaml"
   ```
   This adds a `symbol` column (via the gene index's own `hgnc_id`, mapped to an
   approved symbol through the HGNC TSV), sorts by it, and drops rows with no
   approved symbol or whose symbol is shared by more than one gene id (reported, and
   refused outright below `--min-mapping-rate`, default 0.9). With
   `--manifest-in`/`--manifest-out` it also writes a rerank-ready manifest pointing at
   the new table with an ABSOLUTE `table:` path, carrying over the original's `prior`
   and `cohort_axes` verbatim and setting `key`/`key_column` to `symbol`.
3. **Feature axes**, one TSV per axis, each as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-structure --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" features \
       --ht <gene-index.ht> --fields cds_length,n_coding_exons,n_transcripts \
       --hgnc-tsv <hgnc_complete_set.txt> \
       --restrict-to "$HVANTK_WORK/out/rerank/cohort_symbols.tsv" --restrict-column symbol \
       --out "$HVANTK_WORK/out/rerank/structure.tsv"
   ```
   `--ht` is any gene-keyed hvantk table (`references/feature-axes.md` for constraint
   and expression axes); `--restrict-to` is what makes the universe equal the cohort's
   tested genes here -- drop it to rank the whole gene index instead. Every axis's
   feature columns must have globally unique names (`references/feature-axes.md`).
4. **Labels**, as a job each: build the gene-disease table, extract a gene-set
   collection from it, then flatten to a flat symbol list restricted to the universe
   step 2 fixed.
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-clingen --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess clingen:gene-disease --raw-dir "$HVANTK_WORK/raw/clingen" \
       --output "$HVANTK_WORK/out/rerank/clingen.ht" --skip-download --no-check-drift
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-geneset --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     hvantk genesets clingen --ht "$HVANTK_WORK/out/rerank/clingen.ht" \
       --min-classification Moderate -o "$HVANTK_WORK/out/rerank/labels.json"
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-labels --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/prepare-rerank-inputs.py" labels \
       --geneset "$HVANTK_WORK/out/rerank/labels.json" \
       --restrict-to "$HVANTK_WORK/out/rerank/cohort_symbols.tsv" --restrict-column symbol \
       --out "$HVANTK_WORK/out/rerank/labels.txt"
   ```
   `--min-classification` is "this level or stronger" (`Moderate` means Moderate,
   Strong or Definitive); omit `--group` to union every ClinGen GCEP group, i.e. any
   disease, or pass one or more `--group` names to narrow it. `genesets gencc` is the
   GenCC equivalent (`--group-by submitter`). Restrict labels to the universe now, not
   by lowering `min_label_coverage` later (Traps).
5. **The config.** Exactly `name`, `cohort`, `features` (a list of `{name, path}`),
   `labels` (`{path}`), and optionally `min_label_coverage` (never below 0.5 -- Traps)
   -- no `prior:` key, no other key. `cohort` points at the symbol-keyed manifest from
   step 2; `features` lists each axis TSV from step 3 under the name used in
   `--out`'s basename.
6. **Run**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-rerank --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk rerank -c "$HVANTK_WORK/out/rerank/config.yaml" -o "$HVANTK_WORK/out/rerank/ranked.tsv" \
       --blocks <hgnc_complete_set.txt> --seed 42
   ```
   `--blocks` switches cross-validation to paralogue-blocked folds: an HGNC-symbol
   universe can always be blocked, and an unblocked AUC on a gene universe that
   includes paralogue families overstates how well the model would generalise to a
   gene it has never seen a relative of (`references/inputs-and-identifiers.md`). A
   cohort whose declared axes include `n_case_var`, `conc` and `driver_af` (the
   architecture axis step 2 carries over) auto-wires the advisory
   `CaseControlArchitectureAudit` -- there is no `audit:` config key.
7. **Evaluate**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-eval --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/rerank-evaluate/templates/evaluate-ranking.py" \
       --ranked "$HVANTK_WORK/out/rerank/ranked.tsv" --receipt "$HVANTK_WORK/out/rerank/eval.json"
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
  universe -- note that `b` here is the FULL label count, so a label outside the
  universe still counts against coverage, it is not simply ignored
  (`references/inputs-and-identifiers.md`). Fix: restrict the label file to the
  universe (step 4's `--restrict-to`); never lower `min_label_coverage` to get past
  it -- that hides the mismatch rather than fixing it.
- A cohort keyed by Ensembl `gene_id` scored against symbol-keyed features (or the
  reverse). Cause: `hvantk rerank` does no identifier mapping at all. Fix: map the
  cohort table to the feature axes' space first (step 2), before writing the config.
- `ValueError: none of the N gene(s) in this universe appear in the gene-group
  table`, raised by `--blocks` before any `Blocks:` line ever prints. Cause: the
  ranked units are not HGNC approved symbols (`load_gene_groups` matches the
  approved `symbol` column only -- not aliases, not Ensembl ids). Fix: map to
  symbols before running; `references/inputs-and-identifiers.md`.
- A logged warning `gene-group table matched only k/N genes (X%) in this universe`
  (below 50%), with the run otherwise completing. Cause: the same id-space mismatch
  as above, just not total -- the unmatched genes became unconstrained singletons,
  which quietly degrades "blocked" toward unblocked. Fix: the same mapping as above;
  a near-100% match rate is the expectation for an HGNC-symbol universe.
- `DominantBlockError: the largest paralogue block ... holds .../... (...%), over
  the ... ceiling`. Cause: one gene family exceeds `--max-block-frac` (default 0.1)
  of the universe -- plausible on a SMALL universe, where one family can easily be
  a tenth of everything ranked, not just on a genome-wide one. Fix: raise
  `--max-block-frac` deliberately (and say so in the run notes), or accept a smaller
  universe is simply family-dominated.
- A blocked-fold abort from `_grouped_splits`: `only N paralogue block(s) for 5
  folds`, `blocked fold i has no test example`, `blocked fold i has no training
  example at all`, or `blocked fold i has no {positive,negative} training example`.
  Cause: the paralogue blocks split too unevenly across 5 folds for this universe's
  size or label balance. Fix: a bigger universe or a broader label set; lowering the
  fold count below 5 is Python-API-only (`Config.folds` -- there is no `folds:` key
  in a CLI config.yaml at all).
- `too few examples for 5-fold OOF: ...`. Cause: fewer than 5 positive or negative
  units in the universe. Fix: a broader label set or a bigger universe (through the
  CLI); a smaller fold count is Python-API-only, same as above -- never a smaller
  `min_label_coverage`.
- `Audit: 0 gene(s) FLAGGED for review ... Reasons: {}` when you expected the
  architecture audit to flag something. Cause: that line prints the same whether
  `NoAudit` ran or the real audit found nothing (`references/audit-and-leakage.md`).
  Fix: carry the original cohort's architecture axis through step 2 unchanged.
- Duplicate symbols after mapping. Cause: two gene ids map to the same approved
  symbol. Fix: `prepare-rerank-inputs.py` drops both sides of a shared symbol and
  refuses outright (exit 1) if that -- or any other unmapped-symbol reason --
  pushes the mapped fraction below `--min-mapping-rate` (default 0.9); the cohort
  loader also refuses a table with more than one row per gene key outright
  (`duplicate row(s) for gene key(s) ...`) if any slipped through regardless.
- Every `prior_stat` empty, with `hvantk rerank` otherwise completing normally (no
  error). Cause: the cohort table's key space differs from the feature axes' (e.g.
  mapped into symbols inconsistently, or not mapped at all on one side) -- the left
  join onto the universe then matches nothing. Fix: map both through the SAME
  `prepare-rerank-inputs.py` call chain so they agree by construction; its
  `--min-mapping-rate` floor (step 2) is what catches a wholesale mismatch before
  it ever reaches `rerank`.
- `--n-perm` with one feature axis: rejected before any work starts (the CLI checks
  `len(features) >= 2`) -- the first axis is always the ablation baseline, so a
  one-axis config has nothing left to search a null over.
- `--max-block-frac` without `--blocks`: rejected outright (`--max-block-frac only
  applies to paralogue-blocked folds`).
- A leftover `prior:` key in the config. Cause: a config written before the
  cohort-manifest migration. Fix: delete it -- the prior now lives entirely in the
  cohort manifest's own `prior:` block, referenced via `cohort:`.
- A cohort manifest's `table:` not found, or resolving to someone else's file.
  Cause: `table:` is resolved against the READING process's current working
  directory, never the manifest file's own location -- and a Slurm job's working
  directory is its log directory, never wherever the manifest and its table happen
  to sit. Fix: always write an ABSOLUTE `table:` path (`prepare-rerank-inputs.py
  cohort --manifest-out` already refuses a relative `--out`, for the same reason).
- A computed AUC (`roc_auc_score(y, score)` over the ranked table) differing from
  the numbers in the console's per-axis ablation table. Cause: the ablation's
  per-axis AUCs are UNCALIBRATED subset fits (baseline alone vs. baseline+axis), a
  different quantity from the full model's calibrated, out-of-fold `score` column
  -- and the `Credibility: ROBUST=a / INTERMEDIATE=b` line is a count of genes in
  two tiers, not an AUC at all, so it is not comparable to either. Fix: compute the
  overall AUC from the ranked table only; there is no printed overall-AUC line to
  compare against (`evaluate-ranking.py` computes it).
- A circular axis quietly inflating AUC (built from the same curated database the
  labels come from). `references/audit-and-leakage.md`.

## 6. Verification

1. `config.yaml` parses with exactly `name`, `cohort`, `features`, `labels` and
   optionally `min_label_coverage`; `features` lists at least one axis TSV that
   exists.
2. `ranked.tsv` carries exactly `gene, prior_stat, score, score_percentile, tier,
   verdict, flag, flag_reason, y`; every label-positive gene restricted into the
   universe is a `y = 1` row. There is no printed overall-AUC line on the console
   at all (Traps) -- `evaluate-ranking.py`'s own `auc` (`roc_auc_score(y, score)`
   over the ranked table) is the number to read and to copy into the run notes.
3. The console's `Audit: N gene(s) FLAGGED ... Reasons: {...}` matches the table's own
   `flag`/`flag_reason` counts; `Blocks: ... digest <12 hex>` names how many genes
   matched the gene-group table (most of the universe, for an HGNC-approved-symbol
   universe).
4. `evaluate-ranking.py`'s receipt: `n_genes`, `n_positive`, `auc`, `average_precision`,
   the positive count in the top decile, tier and verdict counts, flag counts by
   reason, and `ranked_sha256`.
5. Run notes carry `commit:`, `container:`, one `job:` line per job submitted, and
   `auc: <value>` from step 7's receipt.
