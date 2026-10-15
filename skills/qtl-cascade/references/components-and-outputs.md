# QTL cascade: components, options, classes, and output layout

`hvantk qtlcascade` (`hvantk/tools/qtl/qtlcascade_cli.py`) has four subcommands.
`run` chains `cascade` + `coloc` + the gene summary + the report; the component
commands below are for driving one stage at a time.

## `cascade` — the join and classification

```
hvantk qtlcascade cascade --eqtl-ht E --pqtl-ht P -o OUT.ht \
  [--tissue T] [--eqtl-p 5e-8] [--pqtl-p 5e-8] [--overwrite]
```

Both `E` and `P` must be keyed `(locus, alleles, gene_id)`. Each table is filtered by
its own threshold first (`p_value <= --eqtl-p` / `<= --pqtl-p`; `0` keeps every row —
note `cascade`'s own CLI default is `5e-8` for both, same as the builders'), then
outer-joined on `(locus, alleles, gene_id)`. The triple key (gene_id included) means a
variant that is an eQTL for one gene and a pQTL for a different gene at the same locus
is never joined — only matching `gene_id` counts as a cascade candidate. If neither
`--tissue` nor a tissue filter was applied upstream and either table carries a
`tissue` field, the join logs a warning: it may cross-multiply rows across tissues.

Each joined row gets a `cascade_class`:

| Class | Condition |
| --- | --- |
| `eqtl_mediated` | both present, `eqtl_beta * pqtl_beta > 0` (same direction) |
| `discordant` | both present, opposite direction |
| `eqtl_only` | eQTL present, pQTL absent |
| `pqtl_only` | pQTL present, eQTL absent |

Output fields (on top of the `(locus, alleles, gene_id)` key): `eqtl_beta`,
`eqtl_se`, `eqtl_pvalue`, `pqtl_beta`, `pqtl_se`, `pqtl_pvalue`, `cascade_class`,
`attenuation_ratio` (`1 - |pqtl_beta / eqtl_beta|`, only for `eqtl_mediated` rows;
missing otherwise), `tissue` (the `--tissue` value given, or missing if none).
Checkpointed to `-o`/`--output`; refuses to overwrite without `--overwrite`.

## `coloc` — colocalisation (ABF)

```
hvantk qtlcascade coloc --eqtl-allpairs E --pqtl-allpairs P \
  --cascade-genes genes.txt -o OUT.tsv [--tissue T] [--window-kb 500]
```

`--cascade-genes` is a file, one gene_id per line — normally the cascade genes with
*both* eQTL and pQTL evidence (`eqtl_mediated` + `discordant`; `pipeline.py`'s own
`_stage_coloc` builds this set that way). For each gene: inner-join the two allpairs
tables on `(locus, alleles, gene_id)` (optionally filtered to `--tissue` first),
take the window within `--window-kb` of the position of the **lead eQTL variant**
(lowest eQTL p-value in that gene), skip the gene if fewer than 2 variants fall in
that window, and run single-variant ABF (Wakefield 2009 / Giambartolomei et al. 2014)
over the window:

- Priors: `p1 = p2 = 1e-4` (causal for one trait only), `p12 = 1e-5` (causal for
  both), prior effect-size variance `W = 0.04`.
- Five posteriors `H0`-`H4` per gene (`H4` = shared causal variant between the two
  traits — true cascade; `H3` = both associated, different causal variants — an LD
  artefact). `H4 > 0.8` is the conventional colocalised threshold
  (`DEFAULT_COLOC_H4_THRESHOLD`, `hvantk/algorithms/qtlcascade/constants.py`).
- `H4` asks only whether the two traits share a causal variant in the window — it
  does **not** care whether the effect directions agree. A `discordant`
  cascade-classified gene can still colocalise with a high `H4`; a concordant
  (`eqtl_mediated`) gene is not guaranteed to.

Output TSV columns: `gene_id`, `tissue` (the `--tissue` value, or `"unknown"` if
none), `H0`, `H1`, `H2`, `H3`, `H4`, `n_variants`.

## `run` — the full per-tissue pipeline

```
hvantk qtlcascade run --eqtl-ht E --pqtl-ht P \
  [--eqtl-allpairs EA --pqtl-allpairs PA] \
  [--constraint-ht C] [--disease-genes-ht D] [--tissues a,b,c] \
  -o OUT_DIR [--eqtl-p 5e-8] [--pqtl-p 5e-8] [--window-kb 500] \
  [--no-plots] [--no-report] [--overwrite] [--dry-run]
```

`--eqtl-allpairs`/`--pqtl-allpairs` are both-or-neither (`Both eqtl_allpairs_ht and
pqtl_allpairs_ht are required for colocalization (or omit both to skip coloc).`).
`--constraint-ht` (a gnomAD constraint table keyed by `gene_id`, field
`oe_lof_upper`) and `--disease-genes-ht` (any table keyed by `gene_id`) are optional
overlays added to the gene summary. `--tissues` (comma-separated) runs one tissue at
a time and writes a combined report; omit it to run once with no tissue filter.
`--dry-run` prints the plan (tissues, coloc on/off, output dir, the stage list) with
no Hail work. With one tissue filter, `run` prints `Cascade genes: <n>` and
`Class counts: <dict>`; with `--tissues` (even a single name), it prints
`Completed: <n> tissue(s)` and, per tissue, `  <tissue>: <n> genes, <m> colocalised`.

### Per-run output layout

```
OUT_DIR/
  per_tissue/<tissue>/
    cascade.ht             # keyed (locus, alleles, gene_id); see `cascade` above
    gene_summary.ht
    gene_summary.tsv       # same rows, tab-separated
    coloc_results.tsv      # only when coloc ran for this tissue
  plots/
    <tissue>_cascade_classes.png
    <tissue>_coloc_posteriors.png
    cross_tissue_heatmap.png   # only with >= 2 tissues and plots enabled
  qtlcascade_report.html
```

`<tissue>` is sanitised (non `[A-Za-z0-9_-]` characters become `_`); plain tissue
names such as `Liver` are unchanged.

### Gene summary columns

`gene_summary.ht`/`.tsv`, keyed `gene_id`: `n_eqtl_variants`, `n_pqtl_variants`,
`n_concordant` (`eqtl_mediated` rows), `n_discordant`, `best_eqtl_pvalue` (min),
`best_pqtl_pvalue` (min), `has_complete_cascade` (`n_concordant > 0`); `tissues` (set
of tissue values seen, when the cascade table carries one). When `--constraint-ht`
was given: `oe_lof_upper`. When `--disease-genes-ht` was given: `is_disease_gene`.
When coloc ran: `coloc_max_h4` (max `H4` across tested tissues for that gene) and
`coloc_n_tissues` (count of tissues where `H4 > 0.8`).

## `report` — HTML from existing results

```
hvantk qtlcascade report [--gene-summary TSV] [--coloc-results TSV] \
  [--plots-dir DIR] -o OUT.html [--title "..."]
```

Standalone equivalent of `run`'s report step, from already-built TSVs (reads every
`*.png` in `--plots-dir` and embeds it as base64). Section headings, in order, each
only present when its data is: `Overview` (tissues, class counts and their
percentages), `Top Cascade Genes (by concordant count)` (up to 50 rows, by
`n_concordant` when that column exists), `Colocalization Results` (genes tested,
count and percent with `H4 > 0.8`, median `H4`), `Plots`, `Methods` (always present).
The report's footer always carries its generation timestamp
(`Generated on <timestamp> by hvantk qtlcascade`) — never compare report *files*
for reproducibility; compare the TSVs and `summarize-cascade.py`'s digests instead.
