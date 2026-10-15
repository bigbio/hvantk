---
name: cohort-attach-burden
description: 'Use when turning a case/control cohort into a gene-level table hvantk can join and rank: running the gene burden test on a cleaned MatrixTable (hvantk cohort burden), writing or checking the cohort manifest (hvantk cohort validate), and attaching the cohort''s gene table to the gene index or the composed annotation matrix (hvantk cohort attach) as Slurm jobs, or attaching a gene-level table you already have. Also when attach stops with "mapped ... below the required", "duplicate row(s) for gene key", "missing declared column(s)" or "needs --hgnc", or burden rejects multi-allelic rows or an undefined case/control field.'
user-invocable: true
argument-hint: "<cohort.mt | gene-table.tsv> [--layer1 gene-index.ht] [--key gene_id|hgnc_id|symbol]"
metadata:
  stage: analysis
  status: provisional
  scale: node
---

# Attach a cohort's gene-level burden to the gene index

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

Two entry points, depending on what you start from:

- **A cleaned, annotated case/control MatrixTable.** Run the burden test first
  (`hvantk cohort burden`), which also emits a ready-to-edit cohort manifest.
- **A gene-level table you already have** (an external cohort's own burden or
  association output). Write the manifest by hand and start from validation.

The outcome either way: the gene-level table (a TSV, from burden or supplied), the
cohort manifest describing it, the table attached onto a gene-keyed Layer-1 matrix with
its coverage report, and a receipt.

Not this skill: QC and variant-level annotation of the MatrixTable (`joint-genotyping`,
or your own annotation pipeline) -- burden assumes a clean, already-annotated input;
building the gene index or the composed Layer-1 matrix (`gene-index`); ranking genes
once a cohort is attached (`rerank-evaluate`).

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value.
2. Checkout and image: confirm the commit and the image digest, as every skill does.
3. The MatrixTable contract, if you are starting from one (`references/mt-contract.md`):
   split bi-allelic rows, a gene field in one identifier space
   (`gene_id`/`hgnc_id`/`symbol`) and unversioned, a route field (a string or array of
   strings -- the variant classes burden tests separately and then picks the best of),
   and a boolean case/control column defined for every sample. Check it with
   `inspect-cohort-mt.py` as a job *before* spending a burden job on an MT that will
   just reject it.
4. The Layer-1 table you will attach onto (the composed annotation matrix, or the gene
   index itself -- `attach --layer1` accepts any table keyed by `gene_id` that spans
   the genes you care about) and its path.
5. Whether an HGNC table is needed: `--hgnc <hgnc:lookup .ht>` is required when the
   cohort's `key` is `hgnc_id` or `symbol`, or when the manifest declares `labels`
   (label gene sets are always symbol-keyed, regardless of the cohort's own key).
6. If you are starting from an existing table instead of running burden: confirm you
   have the right to use it, and which column is the prior statistic and in which
   direction (`lower_is_better`/`higher_is_better`) -- never guess the direction.

## 3. Procedure

1. **Inspect the MatrixTable** (job, only if starting from one):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-inspect --cpus 4 --mem-gb 16 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/cohort-attach-burden/templates/inspect-cohort-mt.py" \
       --mt "$HVANTK_WORK/cohort.mt" --gene-col gene_id --route-col route --arm-col is_case
   ```
2. **Run burden** (job):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-burden --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
     hvantk cohort burden --mt "$HVANTK_WORK/cohort.mt" --gene-col gene_id --route-col route \
       --arm-col is_case --key gene_id --carrier-mode het --mtc bh \
       --output "$HVANTK_WORK/out/gene_burden.tsv" --emit-manifest "$HVANTK_WORK/out/cohort.yaml"
   ```
   Carrier modes are explained in `references/mt-contract.md`. `--mtc` adds `p_adj`
   (Benjamini-Hochberg or Bonferroni); the prior column (`minp`) always stays the raw,
   uncorrected min-p -- `rerank` and the architecture audit read `minp` directly.
3. **Edit the manifest** (`references/cohort-manifest.md`): the emitted stub names the
   output's filename stem; give it a real `name`. Check `key`/`key_column`,
   `min_mapping_rate`, the `architecture` axis (`n_case_var`/`conc`/`driver_af` --
   this is what `rerank-evaluate`'s case/control audit reads), and add `labels` if you
   plan to rerank. Starting from an existing table instead: write this file by hand,
   pointing `table` at it and naming the right prior column and direction.
4. **Validate, on the login node** (no job needed -- pure Python, header-only):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk cohort validate --cohort "$HVANTK_WORK/out/cohort.yaml"
   ```
   This only reads the table's header -- it proves the declared columns exist, not
   that the rows are clean (see Traps).
5. **Attach** (job):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-attach --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
     hvantk cohort attach --cohort "$HVANTK_WORK/out/cohort.yaml" --layer1 "$HVANTK_WORK/gene-index.ht" \
       --output "$HVANTK_WORK/out/attached.ht" [--hgnc "$HVANTK_WORK/hgnc-lookup.ht"]
   ```
   Copy the job's printed `attach NAME: n/m genes tested -> ...` line into your run
   notes -- it is the one place the tested-gene count is reported in plain text.
6. **Verify** (job):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/cohort-attach-burden/templates/verify-attached.py" \
       --attached "$HVANTK_WORK/out/attached.ht" --report "$HVANTK_WORK/out/attached.ht.report.json" \
       --burden "$HVANTK_WORK/out/gene_burden.tsv" --cohort "$HVANTK_WORK/out/cohort.yaml" \
       --receipt "$HVANTK_WORK/out/attached.receipt.json"
   ```

## 4. Scale-up

Burden's cost grows with samples times variants. For a callset too large for one node,
run burden **per chromosome** (genes never span a chromosome boundary) **without
`--mtc`**, concatenate the resulting TSVs, then apply the multiple-testing correction
exactly once over every gene -- a per-chunk Benjamini-Hochberg is wrong (its denominator
is the total test count). Size: 8 cores / 32 GB is enough for a few thousand exomes;
budget a fat node for a genome-scale callset. Give each chromosome's job its own
node-local scratch. `attach` and `verify-attached.py` are both small and run on the
full gene universe regardless of cohort size -- they never need to be chunked.

## 5. Traps

- Burden refuses multiallelic rows outright (`N multiallelic row(s): ... run
  hl.split_multi_hts upstream`). Cause: the clean-input contract requires split,
  bi-allelic variants. Fix: `hl.split_multi_hts` before handing the MT to this skill.
- `case/control column 'x' must be a defined boolean for every sample`, or `must be
  binary: found N case / M control`. Cause: the arm column has missing values, or is
  all-case/all-control. Fix: define it for every sample upstream; sample QC is out of
  scope here.
- `validate` prints OK, then `attach` fails with `duplicate row(s) for gene key(s)
  ...`. Cause: `validate` only reads the table's header, never a data row -- it cannot
  see duplicate keys. Fix: aggregate the table to exactly one row per gene before
  attaching.
- `cohort table ... has no column 'symbol' (key_column)`. Cause: `key` names the
  identifier *space* (`gene_id`/`hgnc_id`/`symbol`); `key_column` names the *column* in
  the table that actually holds those values, and real tables often differ (a
  symbol-keyed table commonly names that column `gene`). Fix: set `key_column`
  explicitly when it is not literally named after `key`.
- `cohort table ... is missing declared column(s): ...`, or `mapped 0.0x of its
  identifiers`. Cause (the latter): a versioned Ensembl id (`ENSG...\.5`) -- the
  mapper never strips versions -- or `key` naming a different identifier space than
  the column actually holds. Fix: unversion gene ids upstream; double-check `key`
  against the column's real content.
- `needs --hgnc <hgnc:lookup .ht>` on a `gene_id`-keyed cohort. Cause: the manifest
  declares `labels`, which are always symbol-keyed regardless of the cohort's own key.
  Fix: pass `--hgnc` whenever `labels` is declared, not only for a symbol/`hgnc_id`
  cohort.
- The emitted manifest's `prior.column`/`cohort_axes` point at columns a hand-edited
  table does not have. Cause: `--emit-manifest`'s defaults describe burden's own
  output (`minp`, `n_case_var`/`conc`/`driver_af`); an external table rarely matches.
  Fix: re-point `prior.column` at the real statistic, and drop `cohort_axes` entirely
  for a table with no variant-level architecture columns -- do not emit it with the
  wrong names.
- Untested genes read back as zero downstream. Cause: `attach` leaves every cohort
  column *missing* for a gene absent from the cohort table -- `cohort_tested` is the
  only signal of "never tested"; nothing is ever imputed to zero. Fix: filter or
  branch on `cohort_tested`, never treat a missing value as zero.
- A declared column named `cohort_tested`, `label`, or a name that collides with an
  existing Layer-1 field. Cause: `cohort_tested` (always) and `label` (when `labels`
  is declared) are reserved output names `attach` writes itself. Fix: rename the
  cohort column; the error names exactly which name collided and why.
- `attach` silently overwrites a previous run's output. Cause: `--output` is written
  unconditionally, no overwrite flag needed or offered. Fix: give each run its own
  output path if you need to keep the old one.
- `--mtc` applied to a per-chromosome chunk. Cause: Benjamini-Hochberg's correction
  depends on the total number of tests; computing it per chunk under-corrects. Fix:
  concatenate first, correct once, as in Scale-up.
- `hvantk: command not found`, or a `python ...` or `hvantk ...` invocation that looks
  like it ran but touched the wrong checkout. Cause: the console script and the right
  Python (3.10, with Hail) live inside the image, not on the login node. Fix: route
  every command through `run-in-image.sh` or `submit-hail-job.sh`; never
  `python -m hvantk` from another directory.

## 6. Verification

`verify-attached.py`'s receipt (and stdout) carries: `n_genes` and `n_tested` (from the
attach report); `mapping_rate` and `prior_non_null_rate` (also from the report);
the burden table's own gene count, its column list, its minimum `minp`, and how many
genes have `p_adj < 0.05`; the attached table's key and which declared columns are
present; `burden_sha256` (the burden TSV's own file hash) and
`attached_content_sha256` (a digest over the *tested* rows' declared columns only).

Expect: `n_genes` equal to the Layer-1 table's own gene count (every spine gene is a
row, tested or not); `n_tested` at most the burden table's row count (a gene can be
tested without mapping onto the spine, e.g. a stale symbol, but never the reverse);
the mapping rate at or above the manifest's own `min_mapping_rate` (`attach` would
already have raised `MappingRateError` otherwise, so this is a sanity echo, not a new
check). Tie this back to the `.provenance.json` sidecar of whichever Layer-1 table you
attached onto, if you need to know what built it.
