# `gwas-coloc`: GWAS locus x eQTL colocalisation

A separate pipeline from the eQTL/pQTL cascade above: it colocalises one FinnGen GWAS
locus against one eQTL Catalogue dataset, ranks the genes in the window by ABF `PP4`,
and (by default) confirms the top effector with SuSiE/coloc.susie fine-mapping.
Covered here as a procedure step and reference.

## Egress — read this before putting it in a job

`gwas-coloc` streams **both** sides by remote tabix on every call: FinnGen R10
summary statistics and the eQTL Catalogue dataset (`fetch_finngen_region` /
`fetch_eqtl_region`, `hvantk/algorithms/qtlcascade/gwas_coloc.py`). It needs egress
every time it runs, with no local-file substitute for the FinnGen side. The eQTL side
alone accepts a local tabix path or URL in place of a dataset id (`--eqtl`). Run this
on a node that has egress (typically the login node, through `run-in-image.sh`); a
compute node with no egress fails with a connection error, not a helpful message
about the data being missing. The `qtlcascade` tool manifest
(`hvantk/tools/qtl/qtlcascade_cli.tool.yaml`) declares `network: false` for the
WHOLE command group — true for `cascade`/`coloc`/`run`/`report`, wrong for this one
subcommand. Do not trust that manifest field for `gwas-coloc`.

pysam's remote tabix (what both fetch functions use) downloads the remote file's
index (a `.tbi`) into the **current working directory**, not a temp location. Run
`gwas-coloc` from a scratch directory you control — otherwise a stray `.tbi` lands
wherever the command happened to be run from (the checkout, if that is the cwd).

## Command

```
hvantk qtlcascade gwas-coloc --endpoint FINNGEN_ENDPOINT --chrom C --lead POS \
  --eqtl QTD_ID_or_URL_or_local_tabix [--eqtl-study QTS000015] \
  [--window-kb 500] [--gene ENSG...] \
  [--fine-map/--no-fine-map] [--gwas-n N] [--eqtl-n N] \
  [--superpop EUR] [--ld-cache-dir DIR] [--ld-vcf LOCAL.vcf] \
  -o OUT_DIR
```

- `--endpoint`: a FinnGen R10 endpoint code (e.g. `I9_AF`), not a free-text phenotype
  name.
- `--chrom`/`--lead`: GRCh38 chromosome (no `chr` prefix) and 1-based lead position;
  the window is `±window_kb` around it.
- `--eqtl`: an eQTL Catalogue dataset id (e.g. `QTD000251`), a full URL, or a local
  tabix-indexed file — `eqtl_catalogue_url()` passes a value through unchanged if it
  already looks like a URL/path.
- `--gene`: confirm a specific ENSG instead of the ABF-top gene; the verdict and
  `goi_PP4_abf` are then about *that* gene even if it ranks lower than another.
- `--fine-map` (default on) needs `--gwas-n` and `--eqtl-n` — without both, the run
  stops with `fine-mapping requires gwas_N and eqtl_N (or pass --no-fine-map)` before
  touching the network. `--no-fine-map` keeps only the single-variant ABF ranking.
- `--superpop`/`--ld-cache-dir`/`--ld-vcf`: the LD reference for fine-mapping. Default
  is a 1000 Genomes high-coverage GRCh38 panel (unrelated founders of the chosen
  super-population), downloaded and cached; `--ld-vcf` points at a local VCF instead
  (uses every sample in it) for an offline or reproducible run — still needs the GWAS
  and eQTL sides fetched by tabix, so this alone does not make the whole command
  egress-free.

## Minimum shared variants differ from the cascade's own coloc

A gene needs at least `min_snps` (default 20, `DEFAULT_COLOC_MIN_SNPS`) GWAS/eQTL
variants matched by position and allele set within the window to be tested at all —
a much higher bar than the eQTL/pQTL cascade's `coloc` command, which only skips a
gene below 2 shared variants. A locus with a narrow, sparse window can test zero
genes even though both sides returned data.

## Choosing a lead locus: the GWAS Catalog table

`gwas-catalog:associations` (downloader exists; not consumed by any `qtlcascade`
command) is a plain table of published GWAS associations, keyed `(locus, alleles)`
with `alleles = [risk_allele, "N"]` — the ALT allele is ALWAYS the sentinel `"N"`
(`hvantk/skills/gwas_catalog/builder.py:151`; the catalog's risk-allele rows never
carry a real ALT, not just sometimes), and there is no standard error column. Use
it to pick a `--chrom`/`--lead` pair for a trait of
interest (e.g. the most significant published variant for a phenotype you then look
up a matching FinnGen endpoint for) — it is a locus-selection aid, not an input
`gwas-coloc` itself reads.

## Output

Per call: `coloc_<endpoint>_<chrom>_<lead>.tsv` (columns `gene_id`, `gene`, `n_snp`,
`eqtl_min_p`, `PP3`, `PP4`, sorted by `PP4` descending) and
`report_<endpoint>_<chrom>_<lead>.json`, both under `-o`/`--output-dir`. The JSON
report carries `run_date` (today's date) — never diff two reports for
reproducibility on that basis alone. Its top-level keys: `endpoint`, `region`,
`genome_build` (always `"GRCh38"`), `run_date`, `gwas` (source, `N`,
`min_p_in_region`), `eqtl` (source, `N`), `coloc_params` (`W1`, `W2`, `p1`, `p2`,
`p12`, `window_kb`, `min_snps`), `results` (`n_genes_tested`, `top_effector`,
`top_PP4`, `gene_of_interest`, `goi_PP4_abf`, `coloc_table_tsv`), `fine_map` (`null`
when `--no-fine-map`, otherwise `available`, `n_variants`, `credible_sets_gwas`,
`credible_sets_eqtl`, `coloc_susie_PP4`, LD-consistency diagnostics, `ld_panel`,
`note`), `verdict`, `provenance_notes`, `software`.

`verdict` is one of: `NO COLOC (...)` (no gene of interest resolved, or its ABF PP4
is below the 0.5 gate), `SUGGESTIVE (ABF only; ...)` (ABF clears the gate but
fine-mapping was not run or was unavailable), `INCONCLUSIVE (...)` (fine-mapping ran
but produced no PP4), `CONFIRMED (ABF + fine-mapping agree)`, or `REFUTED (...)`
(fine-mapping disagrees with the ABF ranking — a likely single-variant-ABF artefact
or LD confound). The verdict's PP4 gate (0.5) is deliberately looser than the
eQTL/pQTL cascade's own colocalised threshold (`H4 > 0.8`) — coloc.susie tends to be
more conservative than single-variant ABF.
