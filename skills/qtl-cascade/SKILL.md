---
name: qtl-cascade
description: "Use when running hvantk's molecular QTL cascade for one or more tissues: building the eQTL and pQTL tables with the gtex-eqtl and pqtl plugins, classifying shared variants (eQTL-mediated, discordant, eQTL-only, pQTL-only), colocalising eQTL and pQTL signals per gene, producing the HTML report, or colocalising a GWAS locus with eQTLs (hvantk qtlcascade gwas-coloc). Also when coloc returns no genes, a tissue filter yields zero rows, the run stops with \"Both eqtl_allpairs_ht and pqtl_allpairs_ht are required\" or \"No eQTL parquet files matched\", or the pQTL build asks for hgnc_ht."
user-invocable: true
argument-hint: "--tissues T[,T...] [--eqtl-raw DIR --pqtl-raw DIR] [--out DIR]"
metadata:
  stage: analysis
  status: provisional
  scale: node
---

# Run the QTL cascade (eQTL -> pQTL)

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

Building, per tissue, the eQTL and pQTL Hail Tables with the `gtex-eqtl` and `pqtl`
plugins, joining them into the QTL cascade (`eqtl_mediated`, `discordant`,
`eqtl_only`, `pqtl_only`), colocalising eQTL and pQTL signals per gene with the ABF
method, and producing a gene-level summary, a coloc results table, and a self-contained
HTML report. Optionally, `hvantk qtlcascade gwas-coloc` colocalises one GWAS locus
against an eQTL dataset by remote tabix — a separate, egress-needing step this skill
covers as a procedure step and a reference, not as a measured outcome.

Use this skill once the eQTL and pQTL raw files exist on disk (BYO or downloaded
separately) and you want the cascade, the coloc table, and the report. Not this skill:
acquiring GTEx or pQTL raw files (both BYO; the resource specs
`hvantk/skills/gtex_eqtl/SKILL.md` and `hvantk/skills/pqtl/SKILL.md` cover acquisition
and raw-format detail) or designing a fine-mapping model (`gwas-coloc`'s optional
SuSiE confirmation only; `references/gwas-coloc.md`).

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value. Every later command that reads a profile value starts with the same `source … &&`.
2. Checkout and image: confirm the commit and the image digest the user expects, as in every other skill.
3. The contract's four steps, in order: (a) the eQTL table (one plugin build per source), (b) the pQTL table (one plugin build, needing an HGNC table), (c) the cascade join plus coloc plus gene summary (`hvantk qtlcascade run`, or the component commands), (d) the HTML report (part of `run`, or `hvantk qtlcascade report` standalone).
4. Genome build: every eQTL source but `eqtlgen` is GRCh38; `eqtlgen` is GRCh37. Never cascade an `eqtlgen` table against a GRCh38 pQTL table — the join is on `(locus, alleles, gene_id)` and a build mismatch joins to nothing, silently (`references/inputs.md`).
5. Tissue names: the cascade's tissue filter (`--tissue`/`--tissues`) is an exact string match against both tables' `tissue` field. GTEx sub-tissue names (`Liver`, `Colon_Transverse`, `Heart_Left_Ventricle`, ...) and Fang pQTL tissue names (`Colon`, `Heart`, `Liver`, `Lung`, `Thyroid`) are NOT the same vocabulary — some Fang names split across several GTEx sub-tissues. Confirm which GTEx sub-tissue(s) correspond to the pQTL tissue you are using before you build either table (`references/inputs.md`).
6. All-pairs vs significant pairs: coloc needs every variant-gene pair in the regional window, not only genome-wide-significant ones. Confirm the eQTL and pQTL raw files you are building from are all-pairs files, and where to get them when they are not (`references/inputs.md`); a significant-pairs-only table silently produces near-empty coloc windows instead of an error.
7. An HGNC table: the pQTL builder maps Fang gene symbols to Ensembl gene IDs through an HGNC lookup table (built separately with `hvantk reprocess hgnc:lookup`), because the cascade join needs Ensembl IDs on both sides. Confirm one exists or build it first (step 1 of the procedure).
8. Outputs: confirm the output directory is writable from compute nodes and does not already hold a cascade run you want to keep (`--overwrite` replaces it).

## 3. Procedure

1. **Build the HGNC table, if none exists yet** (job, through `submit-hail-job.sh`):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-hgnc --cpus 2 --mem-gb 8 --time 0-00:20:00 --wait -- \
     hvantk reprocess hgnc:lookup --raw-dir <hgnc raw dir> --output <hgnc>.ht \
       --skip-download --no-check-drift
   ```
2. **Build the eQTL table(s)** (one job per tissue/source; all pairs, p_threshold=0 —
   `cascade` applies its own `--eqtl-p` later, so the same table serves both the
   cascade input and the coloc input):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-eqtl --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess gtex-eqtl:eqtls --raw-dir <eqtl raw dir> --output <eqtl>.ht \
       --skip-download --no-check-drift \
       --plugin-arg source=gtex_v8 --plugin-arg tissue=<Tissue> --plugin-arg p_threshold=0
   ```
   `source` is one of `gtex_v11` (parquet), `gtex_v8` (TSV), or `eqtlgen` (TSV,
   GRCh37, whole-genome in one file — no `--plugin-arg tissue=` needed for it).
3. **Build the pQTL table** (one job; all pairs — no `p_threshold` plugin-arg at all;
   `hgnc_ht` is required unless you deliberately want a non-cascade, symbol-keyed
   table):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-pqtl --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess pqtl:metrics --raw-dir <pqtl raw dir> --output <pqtl>.ht \
       --skip-download --no-check-drift \
       --plugin-arg tissue=<Tissue> --plugin-arg hgnc_ht=<hgnc>.ht
   ```
4. **Dry-run the cascade plan** on the login node, through `run-in-image.sh` (no Hail
   work happens; it only prints the plan):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk qtlcascade run --eqtl-ht <eqtl>.ht --pqtl-ht <pqtl>.ht \
       --eqtl-allpairs <eqtl>.ht --pqtl-allpairs <pqtl>.ht \
       --tissues <Tissue> -o <out dir> --dry-run
   ```
5. **Run the cascade** (job; omit `--eqtl-allpairs`/`--pqtl-allpairs` together to skip
   coloc — both or neither, never one):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-cascade --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
     hvantk qtlcascade run --eqtl-ht <eqtl>.ht --pqtl-ht <pqtl>.ht \
       --eqtl-allpairs <eqtl>.ht --pqtl-allpairs <pqtl>.ht \
       --tissues <Tissue> -o <out dir> \
       [--eqtl-p 5e-8] [--pqtl-p 5e-8] [--window-kb 500] [--no-plots] [--no-report] [--overwrite]
   ```
   Component commands (`cascade`, `coloc`, `report`), the output layout and the gene
   summary columns are in `references/components-and-outputs.md`.
6. **Summarise the run** (job): `summarize-cascade.py` reads the per-tissue outputs
   and prints and records the class counts, the coloc summary and digests of the two
   TSVs:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-summary --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/qtl-cascade/templates/summarize-cascade.py" \
       --out-dir <out dir> --tissue <Tissue> --receipt <out dir>/summary.json
   ```
   Copy its printed counts into the run notes.
7. **Optional: colocalise a GWAS locus** (`hvantk qtlcascade gwas-coloc`) on a node
   with egress — it streams FinnGen and the eQTL Catalogue by remote tabix, so it
   cannot run on a compute node with no egress. `references/gwas-coloc.md` covers its
   options, the fine-mapping inputs, the LD reference, and using the GWAS Catalog
   table to choose a lead locus.

## 4. Scale-up

All-pairs eQTL and pQTL tables are large (on the order of 10^8 rows for a GTEx
tissue): build each tissue's table as its own job, 8-16 cores and 64 GB or more for
the v11 parquet source (Spark-backed), node-local scratch. Several tissues can share
one `qtlcascade run --tissues a,b,c` invocation, or one job per tissue when sizes
differ a lot; a multi-tissue run writes a cross-tissue heatmap alongside each
tissue's own outputs when at least two tissues complete. Coloc runs gene by gene over
the variants in `--window-kb` of the lead eQTL variant, so its memory follows the
window size and the local variant density, not the full table size. Plot and report
generation are cheap and run inline. `--overwrite` replaces an existing output
directory wholesale; without it, a rerun into the same directory fails on the
Hail Table checkpoint.

## 5. Traps

- Zero cascade rows, or far fewer than expected, for a tissue. Cause: the eQTL and
  pQTL tables use different tissue vocabularies and the join's tissue filter is exact
  `==`. Fix: confirm the GTEx sub-tissue name(s) that correspond to the pQTL tissue
  before building either table (`references/inputs.md`); there is no built-in mapping
  applied automatically.
- Coloc returns no genes, or H4 values that look meaningless, even though the cascade
  itself has rows. Cause: one or both allpairs tables are significant-pairs-only, so
  the regional window around the lead eQTL variant has fewer than 2 shared variants
  per gene (`run_coloc_per_gene` skips a gene below 2). Fix: rebuild from an all-pairs
  source; `references/inputs.md` says where to get one per eQTL source.
- `Both eqtl_allpairs_ht and pqtl_allpairs_ht are required for colocalization (or omit
  both to skip coloc).` Cause: `--eqtl-allpairs`/`--pqtl-allpairs` were given for only
  one side. Fix: pass both, or neither (to run the cascade without coloc).
- The pQTL build errors asking for `hgnc_ht`, or raises because no symbol mapped.
  Cause: Fang pQTL rows are keyed by gene symbol, and the cascade join needs Ensembl
  IDs; the builder refuses to build a cascade-incompatible, symbol-keyed table unless
  told to. Fix: pass `--plugin-arg hgnc_ht=<hgnc:lookup .ht path>`; or, only for a
  deliberately non-cascade use, `--plugin-arg no_gene_map=true` (that table will never
  join an eQTL table, which is keyed by Ensembl ID).
- `No eQTL parquet files matched (tissue=T)` or the pQTL equivalent for TSV files.
  Cause: the tissue filter matches the filename prefix before the first dot, exactly
  — a spelling or case difference, or a non-pQTL file sharing the same extension in
  the raw directory, silently fails or picks the wrong file. Fix: check the raw
  directory holds one file per tissue of the right kind and that `--plugin-arg
  tissue=` spells the prefix exactly as named.
- eQTLGen (GRCh37) joined against a GRCh38 pQTL or eQTL table: the cascade runs
  without error and produces a near-empty or empty result. Cause: the join key
  includes `locus`, and a locus on the wrong reference genome never equals one on the
  other. Fix: never mix `eqtlgen` with a GRCh38 table in one cascade run.
- `'JavaPackage' object is not callable` building a `gtex_v11` (parquet) table. Cause:
  a Spark session was started before Hail initialised one of its own. Fix: build
  through `hvantk reprocess`, never through a hand-started `SparkSession` — the
  builder calls `init_hail()` before touching Spark for exactly this reason.
- The HTML report differs between two runs that otherwise reproduced the same
  numbers. Cause: its footer always carries the generation timestamp. Fix: compare
  the gene summary and coloc TSVs (and `summarize-cascade.py`'s digests), never the
  report file's bytes, when checking reproducibility.
- `gwas-coloc` fails with a connection error. Cause: it streams FinnGen and the eQTL
  Catalogue by remote tabix on every call; a compute node with no egress cannot reach
  either. Fix: run it on a node with egress (`references/gwas-coloc.md`); there is no
  local-file substitute for the GWAS side, though the eQTL side accepts a local tabix
  path.
- A stray `.tbi` file appears in whatever directory `gwas-coloc` (or anything else
  calling `fetch_eqtl_region`/`fetch_finngen_region`) was run from. Cause: pysam's
  remote tabix downloads the remote file's index into the current working directory,
  not a temp location. Fix: run it from a scratch directory you control, never from
  the checkout.
- `fine-mapping requires gwas_N and eqtl_N (or pass --no-fine-map)`. Cause:
  `--fine-map` is the default and needs both sample sizes to run SuSiE. Fix: pass
  `--gwas-n`/`--eqtl-n`, or `--no-fine-map` to keep only the single-variant ABF
  ranking.
- An existing output directory. Cause: `qtlcascade run` and `cascade`/`coloc` refuse
  to overwrite a checkpoint by default. Fix: `--overwrite`.

## 6. Verification

After a real run, for each tissue: the per-tissue directory holds `cascade.ht` (with
`_SUCCESS`, keyed `(locus, alleles, gene_id)`), `gene_summary.ht` and
`gene_summary.tsv`, and — when coloc ran — `coloc_results.tsv`; the run directory
holds `qtlcascade_report.html` naming every tissue analysed. `summarize-cascade.py`
(step 6) prints, per tissue: the cascade class counts (they sum to the cascade
table's row count), the number of cascade genes (`gene_summary.tsv`'s row count), how
many of them have a complete cascade, how many genes were coloc-tested and how many
colocalised (`H4 > 0.8`, at most the tested count), the top five genes by `H4`, and
the sha256 of `gene_summary.tsv` and `coloc_results.tsv`. Copy these into the run
notes next to `commit:`, `container:` and the `job:` lines — they are the receipt a
later run can be compared against.
