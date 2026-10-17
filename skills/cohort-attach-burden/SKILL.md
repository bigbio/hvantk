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

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract.

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

1. Site profile:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and its commit: `git -C "$HVANTK_REPO" rev-parse HEAD` and
   `git -C "$HVANTK_REPO" status --porcelain`; compare with the commit the user expects.
3. Image and its digest: `readlink -f "$HVANTK_SIF"` and
   `sha256sum "$(readlink -f "$HVANTK_SIF")"`.
4. Inputs and outputs: confirm where they live and that the output location is writable
   from compute nodes, under `$HVANTK_WORK` or similar, never home.
5. The MatrixTable contract, if you are starting from one (`references/mt-contract.md`):
   split bi-allelic rows, a gene field in one identifier space
   (`gene_id`/`hgnc_id`/`symbol`) and unversioned, a route field (a string or array of
   strings -- the variant classes burden tests separately and then picks the best of), and
   a boolean case/control column defined for every sample. Check it with
   `inspect-cohort-mt.py` as a job *before* spending a burden job on an MT that will just
   reject it.
6. The Layer-1 table you will attach onto (the composed annotation matrix, or the gene
   index itself -- `attach --layer1` accepts any table keyed by `gene_id` that spans the
   genes you care about) and its path.
7. Whether an HGNC table is needed: `--hgnc <hgnc:lookup .ht>` is required when the
   cohort's `key` is `hgnc_id` or `symbol`, or when the manifest declares `labels` (label
   gene sets are always symbol-keyed, regardless of the cohort's own key).
8. If you are starting from an existing table instead of running burden: confirm you have
   the right to use it, and which column is the prior statistic and in which direction
   (`lower_is_better`/`higher_is_better`) -- never guess the direction.

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
       --output "$HVANTK_WORK/out/cohort/gene_burden.tsv" --emit-manifest "$HVANTK_WORK/out/cohort/cohort.yaml"
   ```
   Carrier modes are explained in `references/mt-contract.md`. `--mtc` adds `p_adj`
   (Benjamini-Hochberg or Bonferroni); the prior column (`minp`) always stays the raw,
   uncorrected min-p -- `rerank` ranks genes on the prior column directly, while the
   separate case/control architecture audit never looks at `minp` at all -- it reads
   only `n_case_var`, `conc` and `driver_af` (`references/cohort-manifest.md`).
3. **Edit the manifest** (`references/cohort-manifest.md`): the emitted stub names the
   output's filename stem; give it a real `name`. Check `key`/`key_column`,
   `min_mapping_rate`, the `architecture` axis (`n_case_var`/`conc`/`driver_af` --
   this is what `rerank-evaluate`'s case/control audit reads). A `labels` block is
   optional and does NOT feed `rerank-evaluate` -- that reads labels only from its
   own rerank config's `labels.path` (`references/cohort-manifest.md`). Declaring
   `labels` here only makes `attach` add its own `label` bool column to
   `attached.ht` (for inspection) and makes `--hgnc` required. Starting from an
   existing table instead: write this file by hand, pointing `table` at it and
   naming the right prior column and direction.
4. **Validate, on the login node** (no job needed -- pure Python, header-only):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk cohort validate --cohort "$HVANTK_WORK/out/cohort/cohort.yaml"
   ```
   This only reads the table's header -- it proves the declared columns exist, not
   that the rows are clean (see Traps).
5. **Attach** (job):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-attach --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
     hvantk cohort attach --cohort "$HVANTK_WORK/out/cohort/cohort.yaml" --layer1 "$HVANTK_WORK/out/gene-index.ht" \
       --output "$HVANTK_WORK/out/cohort/attached.ht" [--hgnc "$HVANTK_WORK/out/hgnc-lookup.ht"]
   ```
   The job prints `attach NAME: n/m genes tested -> ...`; write `n` and `m` into
   your run notes as their own line, `attached: n/m`, not the printed line
   verbatim -- it is the one place the tested-gene count is reported in plain text.
6. **Verify** (job):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/cohort-attach-burden/templates/verify-attached.py" \
       --attached "$HVANTK_WORK/out/cohort/attached.ht" --report "$HVANTK_WORK/out/cohort/attached.ht.report.json" \
       --burden "$HVANTK_WORK/out/cohort/gene_burden.tsv" --cohort "$HVANTK_WORK/out/cohort/cohort.yaml" \
       --layer1 "$HVANTK_WORK/out/gene-index.ht" --receipt "$HVANTK_WORK/out/cohort/attached.receipt.json"
   ```

## 4. Scale-up

Burden's cost grows with samples times variants. For a callset too large for one node,
run burden **per chromosome** (genes never span a chromosome boundary) **without
`--mtc`**, concatenate the resulting TSVs (plain `pandas.concat`, no deduplication
needed -- each gene appears in exactly one chromosome's chunk), then apply the
multiple-testing correction exactly once over every gene -- a per-chunk
Benjamini-Hochberg is wrong (its denominator is the total test count). Do this with
`hvantk.algorithms.burden.fet.apply_mtc(df, method)` (reads the frame's `minp` column,
adds `p_adj`; `method` is `"bh"` or `"bonferroni"`) in a short script run as a job
through `run-in-image.sh` or `submit-hail-job.sh`, e.g.:

```python
import pandas as pd
from hvantk.algorithms.burden.fet import apply_mtc

df = pd.concat(pd.read_csv(p, sep="\t") for p in chunk_paths)
apply_mtc(df, "bh").to_csv(concatenated_output, sep="\t", index=False)
```

Then rewrite `--emit-manifest`'s own `cohort.yaml` so its `table:` points at the
*concatenated* TSV, not at whichever per-chromosome chunk happened to be last --
`--emit-manifest` (or a hand-written manifest) only ever names the one file it was
given. Size: 8 cores / 32 GB is enough for a few thousand exomes; budget a fat node
for a genome-scale callset. Give each chromosome's job its own node-local scratch.
`attach` and `verify-attached.py` are both small and run on the full gene universe
regardless of cohort size -- they never need to be chunked.

## 5. Traps

- Burden refuses multiallelic rows outright: `N multiallelic row(s): the burden op
  requires split, bi-allelic variants. Run hl.split_multi upstream.` (the exact message
  the code prints -- note it says `hl.split_multi`, not the `_hts` variant). Cause: the
  clean-input contract requires split, bi-allelic variants. Fix: run
  `hl.split_multi_hts` before handing the MT to this skill -- it also downcodes `GT`
  the way burden's carrier logic expects, which the bare `hl.split_multi` does not.
- `case/control column 'x' must be a defined boolean for every sample; N/M are missing
  (sample QC is upstream).` Cause: the arm column has missing values for some samples.
  Fix: define it for every sample upstream; sample QC is out of scope here.
- `case/control column 'x' must be binary: found N cases / M controls.` Cause: the arm
  column is all-case or all-control (one of the two counts is zero). Fix: the cohort
  needs both arms represented; this is not something to work around by adjusting a
  threshold.
- A non-boolean arm column crashes instead with a raw Hail error:
  `count_where: parameter 'condition': expected expression of type bool, found
  <StringExpression of type str>` (or whatever type the column actually has). Cause:
  `--arm-col` must be a genuine boolean field, not a string like `"case"`/`"control"`
  or an integer `0`/`1` -- burden never coerces it. Fix: cast it to `bool` upstream, or
  catch this earlier with `inspect-cohort-mt.py`, which checks the arm column's dtype
  before any aggregation touches it.
- `validate` prints OK, then `attach` fails with `duplicate row(s) for gene key(s)
  ...`. Cause: `validate` only reads the table's header, never a data row -- it cannot
  see duplicate keys. Fix: aggregate the table to exactly one row per gene before
  attaching.
- A gene field holding empty strings becomes a gene of its own in burden. Cause:
  burden drops variants whose gene is missing (`hl.is_defined`), and an empty string
  is defined, so every variant with a blank gene is counted under a gene keyed `""`.
  Fix: set blank values to missing upstream (`hl.or_missing(hl.len(g) > 0, g)`);
  `inspect-cohort-mt.py`'s undefined-gene count counts only missing values. A cohort
  TABLE is different: `cohort attach` refuses a blank or `NA` key cell outright
  (since 0.6.0, a clean `UsageError` naming the row count and `key_column`), as
  rerank's pandas loader already does.
- `cohort table ... has no column 'symbol' (key_column)`. Cause: `key` names the
  identifier *space* (`gene_id`/`hgnc_id`/`symbol`); `key_column` names the *column* in
  the table that actually holds those values, and real tables often differ (a
  symbol-keyed table commonly names that column `gene`). Fix: set `key_column`
  explicitly when it is not literally named after `key`.
- `cohort table ... is missing declared column(s): ...`. Cause: a column named in
  `prior.column` or a `cohort_axes` entry does not exist in the table at all --
  commonly a hand-edited manifest that still names burden's own column names
  (`minp`/`n_case_var`/...) over a table that never had them. Fix: match the declared
  column names to what the table actually has, or fix the table.
- `X mapped 0 of N identifiers onto the spine; nothing to prepare. First unmapped:
  ...` (every single identifier failed to map) or, a milder version,
  `X mapped 0.42 of its identifiers, below the required 0.90. Unmapped examples: ...`
  (some mapped, but below `min_mapping_rate`). Cause, either way: a versioned Ensembl
  id (`ENSG...\.5`) -- the mapper never strips versions -- or `key` naming a different
  identifier space than the column actually holds. Fix: unversion gene ids upstream;
  double-check `key` against the column's real content.
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
- A declared column such as `conc`/`driver_af` reads back with Hail type `str`
  instead of the expected `float64` (`verify-attached.py` flags this: it prints each
  declared column's type and warns on a non-numeric one). Cause: in hvantk 0.5.0 and
  earlier, `cohort burden` can write a `NaN` statistic as a blank cell (one variant's
  `driver_af`, or every row's `mean_score_case` when `--score-field` is unset), and
  `cohort attach`'s `hl.import_table` call only recognised the literal string `"NA"`
  as missing -- a blank cell made the whole column fail numeric type inference and
  fall back to `str`. Fixed in 0.6.0 (`attach` now also treats a blank cell as
  missing). Fix: update to 0.6.0 or later, or, on an older one, rewrite blank cells in
  the table to the literal string `NA` before attaching.
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
  like it ran but touched the wrong checkout. Cause and fix: see
  `workflow-conventions`'s Runtime rules (section 4) -- route every command through
  `run-in-image.sh` or `submit-hail-job.sh`; never bare `hvantk` or `python -m hvantk`.

## 6. Verification

`verify-attached.py`'s receipt (and stdout) carries: `n_genes` and `n_tested` (from the
attach report); `mapping_rate` and `prior_non_null_rate` (also from the report);
the burden/gene table's own row count, its column list, the manifest's own prior
column and its minimum value, and how many rows have `p_adj < 0.05` (burden-specific;
null for an externally-supplied table with no such column); the attached table's key,
which declared columns are present with their Hail type, and which of those types are
non-numeric (a warning, not a failure -- see the blank-cell trap above); the
`--layer1` path you passed it (recorded for traceability only -- nothing the attached
table itself carries names which Layer-1 input it came from); `burden_sha256` (the
burden/gene table's own file hash) and `attached_content_sha256` (a digest over the
*tested* rows' declared columns only).

Expect: `n_genes` equal to the Layer-1 table's own gene count (every spine gene is a
row, tested or not); `n_tested` at most the burden table's row count (a gene can be
tested without mapping onto the spine, e.g. a stale symbol, but never the reverse);
the mapping rate at or above the manifest's own `min_mapping_rate` (`attach` would
already have raised `MappingRateError` otherwise, so this is a sanity echo, not a new
check); every declared column numeric, none of them `str`.

None of this skill's own outputs carry a `.provenance.json` sidecar -- unlike most
hvantk build outputs, neither `annotate spine`, `annotate compose` nor
`cohort attach` stamps one (all three use a plain Hail `.write()`). To trace what the
Layer-1 table you attached onto was actually built from, use whichever receipt that
table's own producer leaves instead: the `gene-index` skill's own describe/receipt
step if you attached onto the gene index directly, or `annotate compose`'s own
`<output>.manifest.json` sidecar if you attached onto a composed Layer-1 matrix.
`attach`'s own output (`<output>.report.json`, default path) is the only record this
skill's own steps leave; that is what `verify-attached.py` reads.
