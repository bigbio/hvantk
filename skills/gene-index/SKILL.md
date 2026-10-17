---
name: gene-index
description: 'Use when building the gene index that gene-level annotations, cohort tables and rankings join onto (hvantk annotate spine) from the Ensembl gene-structure and HGNC datasets on a Slurm cluster: staging and building the two inputs, choosing the biotype, checking the gene count and the HGNC mapping rate, writing a receipt, and adding annotation axes onto it with hvantk annotate prepare and compose. Also when the gene index comes out empty, the mapping rate is far below what you expect, or a build stops with "HGNC join changed the row count" or "needs --hgnc".'
user-invocable: true
argument-hint: "[--biotype protein_coding|all|<biotype>] [--out DIR]"
metadata:
  stage: analysis
  status: provisional
  scale: node
---

# Build the gene index

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract.

## 1. Scope and outcome

Building the gene index -- the Hail table every gene-level annotation, cohort-burden
table and ranking joins onto -- with `hvantk annotate spine`, from the Ensembl
gene-structure dataset (`ensembl-gene:structure`) and the HGNC dataset (`hgnc:lookup`),
on a Slurm cluster. The outcome is a Hail table keyed on `gene_id`, one row per gene of
the chosen biotype, carrying the gene name, chromosome, start and end, biotype, CDS
length, coding-exon and transcript counts, the MANE Select transcript, the HGNC id and
HGNC gene groups; plus the gene count and HGNC mapping rate the build command reports;
plus a receipt (`templates/describe-gene-index.py`) binding a content digest of the gene
index to the two input datasets' own `.provenance.json` sidecar digests, because the
gene index itself carries no sidecar of its own.

Use this before `cohort-attach-burden` (which takes the gene index, or a matrix composed
on top of it, as `cohort attach --layer1`), before exporting features for
`rerank-evaluate`, and before `annotate prepare`/`compose` (adding annotation axes onto
the gene index -- `references/annotation-axes.md`).

Not this skill: building `ensembl-gene:structure` or `hgnc:lookup` at scale alongside
other datasets (`dataset-build-at-scale` covers the array-job machinery for that; this
skill's own procedure below covers exactly the two builds the gene index needs);
authoring a new annotation-source plugin (a resource spec under
`hvantk/skills/<provider>/SKILL.md`).

## 2. Preflight

1. Site profile:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and its commit: `git -C "$HVANTK_REPO" rev-parse HEAD` and
   `git -C "$HVANTK_REPO" status --porcelain`; compare with the commit the user expects.
3. Image and its digest: `readlink -f "$HVANTK_SIF"` and
   `sha256sum "$(readlink -f "$HVANTK_SIF")"`.
4. Inputs and outputs: confirm where the two raw directories and the three output paths
   (the two built datasets and the gene index) live, and that they are writable from
   compute nodes. The procedure below writes them under `$HVANTK_WORK/out`; if the user
   passes `--out DIR`, use that directory in place of `$HVANTK_WORK/out` everywhere a
   command below names it.
5. Ensembl release: the checkout pins one release in `hvantk/resources/ensembl_release.py`
   (`ENSEMBL_RELEASE`, GRCh38 only -- `annotate spine` takes no genome argument). Confirm
   it is the release you intend, and that `ensembl-gene:structure` and `hgnc:lookup` were
   both built by THIS checkout -- a stale build from an older pin or an older HGNC dump
   joins silently against the wrong gene set. `builder_commit` on their `.provenance.json`
   sidecars does not help here: it is `null` whenever a dataset is built through
   `hvantk reprocess` (the normal path), so same-checkout is established by the run notes'
   own `commit:` line and the job ids that produced each input, not by the sidecar
   (`references/inputs-and-release.md`).
6. Biotype: default `protein_coding`; `all` keeps every biotype; anything else filters to
   that literal `gene_biotype` string with no validation against the input, so a typo
   (e.g. `protien_coding`) builds a table of zero rows rather than failing. If unsure
   what values exist, list them off the built `ensembl-gene:structure` table first, as a
   job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-biotypes --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python -c "import hail as hl; from hvantk.core.utils.hail_context import init_hail; init_hail(quiet=True); ht = hl.read_table('$HVANTK_WORK/out/ensembl-gene-structure.ht'); print(sorted(ht.aggregate(hl.agg.collect_as_set(ht.gene_biotype))))"
   ```
   rather than guessing.

## 3. Procedure

1. **Stage the two raw inputs on the login node**, through the image:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk download ensembl-structure --raw-dir "$HVANTK_WORK/raw/ensembl-gene"
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk download hgnc --output "$HVANTK_WORK/raw/hgnc/hgnc_complete_set.txt"
   ```
   The two downloaders do not behave the same on a re-run: `ensembl-structure` silently
   reuses the cached pinned file, but `hvantk download hgnc` (and the equivalent
   `reprocess hgnc:lookup ... --skip-build` download stage) exits 1 with "File already
   exists" if its target is already there, even if a previous run left it empty. Test
   `[[ -s FILE ]]` yourself before calling it, or pass `--overwrite`
   (`--plugin-arg overwrite=true` through `reprocess`) to force a re-download either way.
2. **Build the two input datasets, each as a job** (resource notes and sidecar format in
   `references/inputs-and-release.md`):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-ensembl --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess ensembl-gene:structure --raw-dir "$HVANTK_WORK/raw/ensembl-gene" \
       --output "$HVANTK_WORK/out/ensembl-gene-structure.ht" --skip-download --no-check-drift
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-hgnc --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess hgnc:lookup --raw-dir "$HVANTK_WORK/raw/hgnc" \
       --output "$HVANTK_WORK/out/hgnc-lookup.ht" --skip-download --no-check-drift
   ```
   or fan both out as one array through `dataset-build-at-scale` instead of two
   submissions.
3. **Build the gene index as a job**, and record its two reported values in your run
   notes:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-gene-index --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk annotate spine --gene-table "$HVANTK_WORK/out/ensembl-gene-structure.ht" \
       --hgnc "$HVANTK_WORK/out/hgnc-lookup.ht" --output "$HVANTK_WORK/out/gene-index.ht" \
       [--biotype <biotype>]
   ```
   `submit-hail-job.sh --wait` prints the job's log path on stderr (`log: ...`); the
   command's own two lines -- `spine: N genes -> OUT` and `HGNC mapping rate: X%` -- are
   in that log. Write their VALUES into your run notes as their own lines, `genes: N`
   and `mapping_rate: X%`, not the printed lines themselves verbatim.
4. **Describe it and write the receipt**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-describe --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/gene-index/templates/describe-gene-index.py" \
       --gene-index "$HVANTK_WORK/out/gene-index.ht" \
       --gene-table "$HVANTK_WORK/out/ensembl-gene-structure.ht" \
       --hgnc "$HVANTK_WORK/out/hgnc-lookup.ht" \
       --receipt "$HVANTK_WORK/out/gene-index.receipt.json"
   ```
5. **Optional: add annotation axes.** Run `annotate prepare` once per source, then
   `annotate compose` to left-join every prepared axis onto the gene index into the
   Layer-1 matrix `cohort-attach-burden` takes. Each is its own job.
   `references/annotation-axes.md` has the spec format, a worked example, and the
   `--hgnc` rule.

## 4. Scale-up

The gene index itself is small: about 20,000 rows for `protein_coding` (the default),
roughly three times that with `--biotype all` -- one job of 4 cores and 16 GB covers it
with room to spare. The GTF parse inside `ensembl-gene:structure`'s build is the
heaviest of the three jobs (a pure-pandas parse of the whole gene-model file before the
`hl.Table.from_pandas` conversion); budget it a few GB of memory, not dozens.
`hgnc:lookup`'s build is lighter still. Downstream, `annotate prepare`'s cost follows the
SOURCE, not the gene index: a `gene_id`-keyed metrics table is cheap, a `variant`-keyed
aggregate needs a variant-scale job of its own; `compose` is one more light job
regardless of how many axes it joins. Rebuild the gene index whenever the Ensembl pin
(`hvantk/resources/ensembl_release.py`) moves or the HGNC dump is re-fetched, and then
rebuild everything joined onto it in turn -- the gene index is the one place identifier
reconciliation happens, so a stale copy silently reintroduces the mismatches it exists to
remove.

## 5. Traps

- Zero genes with a reported 100% HGNC mapping rate. Cause: `--biotype` is a typo (e.g.
  `protien_coding`) that matches nothing, and the mapping rate is defined as `1.0` for an
  empty table. Fix: list the input's actual `gene_biotype` values first and pass one of
  them, or `all`.
- `hl.read_table` errors on a `.gtf.gz` or `.txt` path. Cause: `annotate spine` takes the
  two BUILT Hail tables, never the raw GTF or HGNC TSV. Fix: build
  `ensembl-gene:structure` and `hgnc:lookup` first (step 2), and pass their `.ht` outputs.
- An existing gene index silently replaced. Cause: `annotate spine` has no overwrite
  guard -- it always writes with `overwrite=True`. Fix: write every build to a new output
  path, never reuse one you want to keep.
- `HGNC join changed the row count (a -> b); the lookup is not unique per gene_id`.
  Cause: on this skill's only input, this guard cannot fire by construction -- if it
  does, `spine.py`'s own dedup step changed (`references/inputs-and-release.md`).
  Fix: report it; a hand-edited input is not the explanation.
- A mapping rate far below the usual protein-coding norm (upper-90s%). Cause: either a
  stale or wrong-vintage HGNC table, or you built with `--biotype all`, where most
  non-coding biotypes have no HGNC record at all. Fix: compare like with like -- rebuild
  HGNC and compare the SAME biotype filter across runs.
- No `.provenance.json` next to `gene-index.ht`. Cause: by design -- `annotate spine`
  writes with a raw `ht.write()`, not `hvantk.core.io.save`, so it never stamps a
  sidecar. Fix: nothing to fix; `templates/describe-gene-index.py`'s receipt is what
  records the gene index's own content digest and the two inputs' sidecar digests
  instead.
- A GRCh37 cohort failing to join onto the gene index. Cause: there is no GRCh37 gene
  index -- the Ensembl pin is GRCh38 only. Fix: lift the cohort over to GRCh38 before
  joining; never try to build a GRCh37 gene index.
- A source with versioned Ensembl ids (e.g. `ENSG00000139618.5`) failing to map onto the
  gene index. Cause: the gene index's `gene_id` carries no version suffix, and nothing in
  `annotate prepare`'s mapper strips one. Fix: strip the trailing `\.[0-9]+$` from the
  source's ids before mapping.
- `entry '<source>' needs --hgnc <hgnc:lookup .ht>` from `annotate prepare`. Cause: the
  spec entry's `key` is `hgnc_id`, `symbol` or `uniprot_id` (or a `variant` aggregate not
  targeting `gene_id`), none of which can resolve onto `gene_id` without the HGNC lookup.
  Fix: pass `--hgnc hgnc-lookup.ht` (the same build from step 2).
- `compose` rejecting a `--prepared` option. Cause: either the spec declares an axis with
  no matching `--prepared axis=path` (missing), or a `--prepared` names an axis the spec
  does not declare (unknown) -- `compose` is strict on both. Fix: pass exactly one
  `--prepared` per spec `layer1` axis, spelled the same as the spec's `axis:` label.
- `hvantk: command not found`, or a write landing in the wrong checkout. Cause and
  fix: see `workflow-conventions`'s Runtime rules (section 4) -- route every
  `hvantk`/`python` invocation through `run-in-image.sh` or `submit-hail-job.sh`,
  never bare.
- A small, consistent number of protein-coding genes (a few dozen) never joining onto a
  cohort's VCF, even though their biotype and chromosome both look right. Cause: those
  genes sit on unplaced or unlocalised scaffolds (`KI270...`, `GL000...`), not on
  `1`-`22`/`X`/`Y`/`MT`, and a VCF naming its contigs `chrUn_KI270728v1`-style never
  matches the gene index's own contig spelling for them. Fix: expected and small; do not
  treat it as a join bug unless the count is much larger than a couple dozen genes.

## 6. Verification

`templates/describe-gene-index.py`'s receipt and stdout carry `n_genes`, `mapping_rate`
(full precision, not the build command's rounded percentage), per-chromosome counts,
`n_without_hgnc` plus up to ten example genes with no HGNC record, the biotypes actually
present, `content_sha256` (a digest over every row and field, since the gene index has no
sidecar of its own), and the two input datasets' `.provenance.json` sha256 digests.
Expect: exactly one biotype present when you filtered to one (`protein_coding` or your
literal choice), several when you built with `all`; chromosomes `1`-`22`, `X`, `Y`, `MT`,
plus a small number of protein-coding genes (a few dozen) on unplaced or unlocalised
scaffolds (`KI270...`, `GL000...` -- Traps); nearly every protein-coding gene mapped to
an HGNC record (a rate well below that means a stale or mismatched HGNC build --
`references/inputs-and-release.md`). Cross-check that the run notes' `genes:` and
`mapping_rate:` lines (the VALUES from the build command's `spine: N genes` and `HGNC
mapping rate: X%` output) agree with the receipt's `n_genes` and `mapping_rate`.
