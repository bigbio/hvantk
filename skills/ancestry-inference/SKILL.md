---
name: ancestry-inference
description: 'Use when samples of a cohort must be assigned a continental ancestry with hvantk: building a labelled reference MatrixTable from 1000 Genomes, running the PCA and random-forest pipeline as a Slurm job, and reading the predictions and their confidence. Also when the pipeline refuses for too few shared variants, mixed genome builds or a missing label column, or when many samples come back unassigned.'
user-invocable: true
argument-hint: "<query MatrixTable> [--reference <labelled MatrixTable>] <output directory>"
metadata:
  stage: cohort
  status: provisional
  scale: node
---

# Ancestry inference against a labelled reference

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract. Jobs run through `hpc-hail-job`.

## 1. Scope and outcome

Each query sample gets a predicted continental ancestry, a probability, and principal
component coordinates, from a model trained on a reference whose samples carry known
labels; the run also reports how many samples stayed unassigned below the confidence
threshold and the cross-validated accuracy on the reference. The outcome is a Hail
table and a TSV of predictions, a statistics file, and optionally an HTML report.

Use it with any reference that shares a genome build with the query; the common case is
the 1000 Genomes high-coverage callset with its super-population labels, which this
skill builds. Not this skill: fine-scale or admixed ancestry (the model assigns one
label), relatedness, and sex inference.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and image as in `hpc-hail-job` (commit, digest, run notes).
3. The query: a MatrixTable keyed `(locus, alleles)` with a `GT` entry (or `LGT` and
   `LA`), on a known build. `hvantk hgc` outputs qualify. Note its build and sample count.
4. The reference: a MatrixTable on the same build with a flat column holding the label
   (`super_pop` in the 1000 Genomes build below). If there is none, build one
   (section 3 step 1). A nested struct does not count as a column; flatten it.
5. Shared variants: the pipeline needs 10,000 by default after merging; a query on a
   chromosome slice or an exome needs `--min-shared-variants` lowered deliberately, and
   the result read with that in mind.

## 3. Procedure

1. Build the reference once per build, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/ancestry-inference/templates/ancestry-job.sh" \
     --prepare --vcf "$HVANTK_WORK/raw/1kg/<chromosome vcf>" --samples-table "$HVANTK_WORK/raw/1kg/samples.ht" \
     --out "$HVANTK_WORK/references/1kg_<build>.mt" --cpus 16 --mem-gb 120 --time 0-01:30:00 --name <run>-prepare --wait
   ```
   `templates/prepare-1kg-reference.py` imports the VCF, keeps biallelic SNPs with a
   minor allele frequency above a threshold, thins to a target count, joins the 1000
   Genomes sample table (`hvantk reprocess onek-genomes:samples` builds it, downloading
   the panel) and writes `super_pop` as a flat column. `references/reference-panel.md`
   says where the public files are and how to combine chromosomes.
2. Run the pipeline as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/ancestry-inference/templates/ancestry-job.sh" \
     --query "$HVANTK_WORK/callsets/<name>/<name>.mt" --reference "$HVANTK_WORK/references/1kg_<build>.mt" \
     --out "$HVANTK_WORK/ancestry/<name>" --cpus 16 --mem-gb 64 --time 0-01:00:00 --name <run>-ancestry --wait
   ```
   It runs `hvantk ancestry-inference -q ... -r ... --ancestry-col super_pop -o <out>/predictions.ht --output-dir <out> --export-tsv --generate-report`,
   plus any pipeline flags you pass with `--extra`.
3. Read `<out>/pipeline_stats.json`: `n_shared_variants`, `n_variants_after_filter`,
   `cv_accuracy`, `n_query_assigned`, `n_query_unassigned`. `references/reading-results.md`
   says what each means and when to distrust a run.
4. Read `<out>/predictions.tsv`: `s`, `predicted_ancestry`, `ancestry_probability`,
   the per-population probabilities and the PCs. Unassigned samples carry a missing
   prediction; do not fill them in by hand.
5. Write the run notes: commit, image digest, job ids, the reference used (path and its
   build), the shared-variant count, the accuracy, and the assigned fraction.

## 4. Scale-up

- Whole-genome references: prepare per chromosome and union the MatrixTables, or
  prepare from a few chromosomes; ancestry at the continental level needs tens of
  thousands of common, well-spread variants, not millions.
- LD pruning is the slow step on large inputs; thinning the reference first keeps it
  short. `--skip-ld-pruning` is for already-pruned inputs only.
- Memory grows with samples times variants in the merged table: 64 GB covers thousands
  of samples against a 20,000-variant reference; a dense exome needs more.
- `--n-pcs` beyond 20 rarely helps at the continental level; `--n-pcs-classify` below
  `--n-pcs` keeps the classifier on the informative components.

## 5. Traps

- "fewer than 10000 shared variants". Cause: a slice, an exome, or a reference on the
  other build. Fix: same build; a richer reference; or `--min-shared-variants` lowered
  deliberately, and the run notes say so.
- "reference genome mismatch". Cause: GRCh37 query against a GRCh38 reference (or the
  reverse). Fix: a reference on the query's build; no liftover inside this skill.
- "column super_pop not found". Cause: the label sits inside a struct or under another
  name. Fix: a flat column (`--ancestry-col` names it).
- Many samples unassigned. Cause: `--min-prob` 0.75 with a reference that does not
  cover the cohort's ancestry, or admixed samples. Fix: read the probabilities; lower
  the threshold only with a reason written down; never force a label.
- `cv_accuracy` of 1.0 on real data. Cause: too few reference samples per population,
  or a synthetic reference. Fix: at least ten samples per label, a real panel.
- The documentation names a `--log-level` flag. Cause: a stale page. Fix: verbosity is
  `hvantk -v ancestry-inference ...`.
- `RuntimeError` about scikit-learn. Cause: running outside the image without the `ml`
  extra. Fix: the image has it; jobs run inside it.
- Related samples in the reference. Cause: the full 3202-sample 1000 Genomes release
  includes trios. Fix: the sample table covers the 2504 unrelated samples; the join
  drops the rest.

## 6. Verification

1. `<out>/predictions.ht/_SUCCESS` and `<out>/predictions.tsv` exist; the TSV has the
   `s`, `predicted_ancestry` and `ancestry_probability` columns and one row per query
   sample.
2. `<out>/pipeline_stats.json` reports `n_shared_variants` at or above the minimum you
   accepted, a `cv_accuracy` you believe (section 5), and `n_query_assigned` plus
   `n_query_unassigned` summing to the query's sample count.
3. When the query has samples of known ancestry, their predictions agree; the eval
   measures exactly this on held-out 1000 Genomes samples.
4. The run notes carry the commit, image digest, job ids, the reference and its build,
   and the three numbers from step 2.
