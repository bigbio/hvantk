# The cohort MatrixTable contract

`hvantk cohort burden` fails loud, before spending any real compute, if the input
MatrixTable does not meet this contract (`hvantk/algorithms/burden/aggregate.py`'s
`assert_clean_mt`, and the same field-name checks in
`hvantk/algorithms/burden/checks.py`). `inspect-cohort-mt.py` checks it too, so you
find out before a burden job runs, not after.

## The four requirements

1. **Split, bi-allelic rows.** Every row's `alleles` must have exactly two entries
   (`hl.len(mt.alleles) == 2`). Run `hl.split_multi_hts` upstream -- burden never does
   this for you, and a multiallelic row makes the per-(gene,route) carrier counts
   ambiguous.
2. **A gene field in one identifier space, unversioned.** The row field you pass as
   `--gene-col` must hold values in exactly one of `gene_id` (Ensembl), `hgnc_id`, or
   `symbol` -- whichever you declare with `--key`. No version suffix stripping happens
   anywhere downstream: `ENSG00000139618.15` never maps onto the spine's
   `ENSG00000139618`, silently becoming an unmapped gene rather than an error.
3. **A route field.** The row field you pass as `--route-col` names the variant class
   (e.g. `lof`, `missense`, `synonymous`) each row belongs to -- a plain string, or an
   `array<str>` when a variant contributes to more than one route (it is then counted
   under each one it is tagged with). Burden tests every route for a gene separately
   and reports the single best (minimum-p) route as that gene's result.
4. **A boolean arm column, defined for every sample.** The column field you pass as
   `--arm-col` must be `True` for every case and `False` for every control, with no
   missing values and not all one value -- sample-level QC and case/control
   assignment are upstream of this skill, not something burden infers.

`GT` (the genotype entry field) and `locus`/`alleles` (the variant key) are also
required but are standard MatrixTable fields you are unlikely to be missing by
accident.

## Carrier modes (`--carrier-mode`)

A sample "carries" a gene[,route] when at least one of its genotypes in that
gene[,route] group qualifies, by one of four rules:

| mode | counts a sample as a carrier when... |
| --- | --- |
| `het` (default) | it has at least one heterozygous genotype in the group |
| `hom` | it has at least one homozygous-alt genotype in the group |
| `chet` | it has two or more distinct heterozygous genotypes in the group (a candidate compound-het) |
| `homs_chet` | either `chet` or at least one homozygous-alt genotype (the union) |

The per-*variant* reductions below (`n_case_var`, `conc`, `driver_af`) use a related
but not identical per-variant carrier rule (`hvantk/algorithms/burden/aggregate.py`'s
`_carrier_variant_expr`): `het` and `chet` both count only a heterozygous genotype as
a variant-level carrier; `hom` counts **only** a homozygous-alt genotype (not het);
`homs_chet` is the one mode that counts either (het **or** hom-alt). (`het` vs `chet`
differ only in the *gene-level* qualifying rule above, not in which genotypes count
as carriers at the variant level.)

## What burden computes

For each (gene, route) pair present in the data: a 2x2 table of distinct carrier
counts (case carriers / case non-carriers / control carriers / control non-carriers,
by the chosen carrier mode), a two-sided Fisher exact test (`scipy.stats.fisher_exact`)
giving a p-value and an odds ratio. For each gene, the route with the smallest p-value
wins (ties broken by route name, so the result is deterministic regardless of
collection order); that route's p-value becomes the gene's `minp`, and its
architecture reductions are carried onto the gene's one output row.

Architecture reductions, computed over the winning route's per-variant carrier counts:

- `n_case_var`: how many distinct variants (in the winning route) have at least one
  case carrier.
- `conc`: concentration of case carriers onto the single most-case-carrying variant --
  `max(case carriers of any one variant) / sum(case carriers across all variants)`.
  Close to 1 means the signal is really one recurrent variant, not a rare-variant
  burden.
- `driver_af`: the control-CARRIER frequency (control carriers of that variant /
  total controls -- *not* an allele frequency) of the single variant that drove
  `conc`, i.e. the max-case-carrier variant. Ties on case-carrier count break toward
  the *highest* control-carrier frequency among the tied candidates -- the
  conservative direction, since this value is read downstream as an upper-bound leak
  alarm.
- `mean_score_case` (optional, needs `--score-field`): the mean of that row field over
  variants with at least one case carrier.
- `frac_case_private`: fraction of the winning route's case-carrying variants that
  have *zero* control carriers.

## Output columns (the burden TSV)

`gene, route, minp, odds_ratio, n_case_var, conc, driver_af, mean_score_case,
frac_case_private` always; `p_adj` is added only with `--mtc {bh,bonferroni}`. `minp`
is always the raw, uncorrected p-value -- the correction never overwrites it, it only
adds the extra column. A gene with zero case- or control-defined counts in every
route it appears in can still produce `NaN` in `conc`/`driver_af`/`mean_score_case`/
`frac_case_private` (division by a zero denominator); this is expected and is not
itself a sign of a broken run.
