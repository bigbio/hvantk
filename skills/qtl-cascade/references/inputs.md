# QTL cascade inputs: sources, layouts, tissues, genome builds

## eQTL sources (`gtex-eqtl:eqtls`, `build_eqtl_associations`)

`--plugin-arg source=` selects one of three, dispatched inside the same builder
(`hvantk/skills/gtex_eqtl/builder.py`); all three end up with the same output schema
(`locus`, `alleles`, `gene_id`, `beta`, `se`, `p_value`, `af`, `maf`, `tissue`,
`source`, `is_cis`), keyed `(locus, alleles, gene_id)`.

| `source` | Raw format | Reference genome | Notes |
| --- | --- | --- | --- |
| `gtex_v11` (default) | Per-tissue parquet (`spark.read.parquet`) | GRCh38 | `af` read directly; `maf = min(af, 1-af)`. Needs Hail initialised before Spark touches the file — see Traps in `SKILL.md`. |
| `gtex_v8` | Per-tissue gzipped TSV (`hl.import_table`) | GRCh38 | Carries `maf` directly, not `af`; both fields are set equal to it (v8 loses ALT-direction information). |
| `eqtlgen` | One gzipped TSV, whole-genome | **GRCh37** | No `tissue`/`source` filename dispatch; `tissue` is hardcoded `"Blood"`. Provides a Z-score, not a `beta`/`se` pair (`beta`/`se` are null). Never cascade this against a GRCh38 pQTL table — the join key includes `locus`, and nothing on GRCh37 equals anything on GRCh38. |

`--plugin-arg tissue=<name>` (for `gtex_v11`/`gtex_v8`) selects one tissue's file by
filename prefix before the first dot (`scan_tissue_files`,
`hvantk/core/utils/qtl_helpers.py`) — matching is exact and case-sensitive; a
misspelled or miscased tissue name raises `No eQTL parquet/TSV files matched`, not a
silent empty table.

`--plugin-arg p_threshold=` filters by `p_value <= p_threshold`; the builder's own
default is `5e-8` for every source (it is one shared parameter applied after the
per-source dispatch), so omitting it still gives you a significant-pairs-only table.
Pass `--plugin-arg p_threshold=0` explicitly to keep every variant-gene pair. Build
the all-pairs table once and let `qtlcascade cascade`/`run` apply its own `--eqtl-p`
at cascade time; the same table then also serves as the `--eqtl-allpairs` coloc
input.

Gene IDs are version-stripped (`ENSG00000000003.15` -> `ENSG00000000003`) for every
source, so they join cleanly against HGNC, Ensembl, and the pQTL table's mapped IDs.

## pQTL source (`pqtl:metrics`, `build_pqtl_metrics`)

Only `source=gtex_fang` (the default) is implemented. Raw format: Fang et al. (2025)
allpairs files, **space**-delimited gzip, columns `gene_name SNP CHR BP A1 NMISS BETA
STAT P` (no SE column — the builder derives `se = |BETA / STAT|`, and drops any row
with `STAT == 0` first, since that SE would be undefined). Real files are named
`<Tissue>.allpairs_nobsGE72.txt.gz`; the same `scan_tissue_files` filename-prefix
convention applies, so always pass `--plugin-arg tissue=<Tissue>` when the raw
directory might hold other `.txt.gz`/`.tsv.gz` files.

**Gene-symbol mapping is required for a cascade-compatible table.** Fang rows carry a
gene *symbol* (`gene_name`); the eQTL side is keyed by Ensembl ID. Pass
`--plugin-arg hgnc_ht=<path to an hgnc:lookup .ht>` (built separately with `hvantk
reprocess hgnc:lookup`); the CLI translates this into the builder's `gene_catalog`
parameter (a `GeneCatalogStreamer`) before the build runs — the builder itself never
takes a bare path. Approved symbols map first; a symbol that is not currently
approved is then looked up among previous symbols, then aliases, accepting a match
only when exactly one gene lists it there. Within one tissue, if more than one source
symbol would land on the same Ensembl ID, the one approved symbol among them wins and
the others keep their own source symbol as `gene_id` instead (logged, not dropped).
Pass `--plugin-arg no_gene_map=true` only for a deliberately non-cascade,
symbol-keyed table — it will never join an eQTL table.

`--plugin-arg p_threshold=` is accepted but has no default (omit it entirely to keep
every pair); there is no separate "all pairs" vs "significant pairs" Fang file the way
GTEx ships both — the one allpairs file the builder reads already carries every pair,
so simply not passing `p_threshold` is enough.

## All-pairs vs significant-pairs — where coloc's input actually comes from

Colocalisation needs every variant-gene pair in the regional window around the lead
eQTL variant, not only the pairs that cleared a significance threshold — `coloc`
skips a gene outright below 2 shared variants, and a significant-pairs-only table
rarely has that many within one window. Where to get an all-pairs table per source:

- `gtex_v8`/`gtex_v11`: the anonymous public GTEx bucket that this skill's downloader
  story covers ships **significant pairs only** for the per-tissue files this skill
  reads (the v11 all-pairs archive is one multi-gigabyte file per release; the v8
  all-pairs folder is empty in that bucket). Treat any GTEx-sourced table built from
  that bucket as significant-pairs-only unless you independently confirm otherwise.
- The eQTL Catalogue (study `QTS000015` = GTEx, by default) serves **all-variant**
  per-tissue summary statistics over a genomic window, served by remote tabix, and is
  what `hvantk.algorithms.qtlcascade.gwas_coloc.fetch_eqtl_region(dataset, chrom,
  start, end, study=...)` reads — `dataset` is an eQTL Catalogue dataset id (e.g.
  `QTD000266` for GTEx liver gene expression) rather than a GTEx tissue name. It
  returns `{gene_id: [((pos, ref, alt), beta, se, pval), ...]}` for the window
  (`gene_id` already version-stripped; rows with `se <= 0` already dropped; no `maf`
  in this payload). This is a region-at-a-time, network-needing path (`gwas-coloc`'s
  own machinery, `references/gwas-coloc.md`), not a bulk downloader — use it to build
  a small, genuinely-all-pairs slice when the GTEx bucket only offers significant
  pairs.
- `eqtlgen`: distributes only its own cis-eQTL file; there is no separate
  "significant" vs "all" split to worry about for this source.
- Fang pQTL: the one allpairs file per tissue already is all-pairs (see above).

## Tissue vocabularies do not line up

GTEx eQTL sub-tissue names (`Liver`, `Colon_Transverse`, `Colon_Sigmoid`,
`Heart_Left_Ventricle`, `Heart_Atrial_Appendage`, ...) and Fang pQTL tissue names
(`Colon`, `Heart`, `Liver`, `Lung`, `Thyroid`) are different vocabularies. A mapping
exists (`FANG_TISSUE_EQTL_MAPPING`, `hvantk/skills/pqtl/shared/constants.py`) but
**nothing in the codebase applies it** — the cascade's tissue filter is a plain `==`
against one string on both tables. `Liver`, `Lung` and `Thyroid` happen to use the
same name on both sides; `Colon` and `Heart` do not (each spans two GTEx sub-tissues).
When they differ, decide yourself which GTEx sub-tissue(s) correspond to the pQTL
tissue, and build (or filter) the eQTL table under that same tissue label before
joining — there is no shortcut flag for this.

## Genome builds

Every source above is GRCh38 except `eqtlgen` (GRCh37). The cascade join is on
`(locus, alleles, gene_id)`; a locus on the wrong reference genome never equals one on
the other, so a build mismatch does not error — it silently joins to nothing (or to
far fewer rows than expected). Confirm both inputs to one cascade run share a genome
build before running it.
