# Audit and leakage controls

`hvantk rerank` ships three controls under three different names, and only one of
them is "the audit" the CLI prints a line about. They do different jobs and have
different reach (CLI vs. Python-API-only); none of them rewrites the ranking, and
none of them substitutes for choosing axes with no circularity in the first place
(`references/feature-axes.md`).

## 1. The case/control architecture audit -- advisory, CLI, per-gene

`hvantk.algorithms.rerank.audit.CaseControlArchitectureAudit` flags individual genes
for a human to look at again; it never removes, reorders or re-scores anything. It
auto-wires exactly when the cohort manifest's declared axis columns (across all of
`cohort_axes`, not the prior) are a superset of `n_case_var`, `conc`, `driver_af`
(`ARCHITECTURE_AUDIT_COLUMNS`, `has_architecture_columns`); otherwise `NoAudit` runs
and the CLI logs a warning if it finds a near-miss (an axis named `architecture`
that is missing one of the three). There is no `audit:` config key -- the CLI decides
this from the cohort manifest alone, every time.

Flags (advisory, last-write-wins in this order):

| flag | condition |
| --- | --- |
| `insufficient_data` | `n_case_var <= 2` |
| `recurrent_variant` | `conc >= 0.6` and `driver_af > 5e-5` |
| `common_driver` | `driver_af > 1e-3` |

`hvantk cohort burden --emit-manifest` writes exactly this axis
(`hvantk/algorithms/burden/manifest.py` `ARCHITECTURE_AXIS_COLUMNS`), so a cohort
built that way carries the audit automatically once its manifest is carried through
`prepare-rerank-inputs.py cohort --manifest-in/--manifest-out` unchanged.

## 2. Paralogue blocks -- `--blocks`, CLI, cross-validation structure

Covered fully in `references/inputs-and-identifiers.md`. In short: `--blocks
hgnc_complete_set.txt` switches 5-fold CV (and the per-axis bootstrap interval) from
random gene folds to folds that keep each HGNC gene family's first-listed group
together, so a paralogue pair can never straddle train/test. `--max-block-frac`
(default 0.1) aborts the run if one family would dominate the universe -- raise it
deliberately, or restrict the universe, never just to make the abort go away.
`--max-block-frac` without `--blocks` is a flag-combination error, checked before any
file is read.

## 3. Presence-leakage and provenance arms -- Python API only, column admissibility

Neither of these has a CLI flag. Both decide WHICH COLUMNS may contribute to the
headline result, before cross-validation ever runs; both are properties of a
feature/label PAIR, never of a feature alone.

**Presence-leakage** (`hvantk.algorithms.rerank.leakage.LeakagePolicy`) bars a column
whose MISSINGNESS -- not its values -- predicts the label (the motivating case: EVE
scores exist only for a well-studied subset of proteins, and "well studied" already
correlates with "known disease gene"). Set via `Config.leakage`. `Config.audit`
defaults to `NoAudit()` when left unset (`Config.__post_init__`) -- the CLI's
auto-wiring of `CaseControlArchitectureAudit` from the cohort manifest (section 1
above) is CLI-only; a direct `Config(...)` construction gets no audit at all unless
`audit=CaseControlArchitectureAudit()` is passed explicitly:

```python
from hvantk.algorithms.cohort.spec import load_cohort
from hvantk.algorithms.rerank.catalog.builders import genelist_labels, table_axis
from hvantk.algorithms.rerank.config import Config
from hvantk.algorithms.rerank.leakage import LeakagePolicy

cohort = load_cohort("cohort_symbols.yaml")
feats = [table_axis("structure", "structure.tsv")]
labels = genelist_labels("labels.txt")

cfg = Config(
    name="my-run",
    features=feats,
    labels=labels,
    cohort=cohort,
    leakage=LeakagePolicy(),  # defaults: q=0.10, min_auc=0.55
)
```

**Provenance arms** (`hvantk.algorithms.rerank.provenance`, via
`hvantk.algorithms.rerank.engine.rerank_arms`) bar a column TRAINED on the same
curated resource the labels were derived from (REVEL, trained on HGMD/ClinVar, is
badly circular against a ClinGen- or GenCC-derived label). Declare both sides and
let `rerank_arms` compute the conflict -- never hand-classify a column as "clean" by
inspection, since the equivalence classes are the part that is easy to get wrong:

```python
from dataclasses import replace
from hvantk.algorithms.rerank.engine import rerank_arms
from hvantk.algorithms.rerank.selection import SelectionPolicy

cfg = replace(
    cfg,
    selection=SelectionPolicy(),  # provenance arms require a selection policy
    feature_provenance={"revel_score": frozenset({"ClinVar", "HGMD"}), "cds_length": frozenset()},
    label_provenance=frozenset({"ClinGen", "GenCC"}),
)
results = rerank_arms(cfg)  # {"clean": RerankResult, "all": RerankResult}
```

`feature_provenance` maps every model column to the sources it was trained on (an
empty `frozenset()` means "trained on nothing label-derived", which is clean;
omitting a column from the dict is UNDECLARED and is barred from the `clean` arm on
the conservative assumption that it conflicts, same as a computed conflict).
`label_provenance` is mandatory once `feature_provenance` is set -- an undeclared
label source would admit every circular predictor into a `clean` arm that then looks
entirely healthy, which is the one failure mode this control cannot recover from if
left unset. `rerank_arms` raises rather than guessing. The `clean` arm is always the
headline; `all` exists only so the circularity channel is a measured number instead
of an assumption.

## Permutation nulls -- `--n-perm`, CLI, multiplicity correction

A different question from all three controls above: not "is this column
admissible" but "given N axes were searched, does the best one's delta survive
correcting for having tried N of them". Needs at least one feature axis besides the
baseline (checked at config time: `--n-perm` with one axis is rejected outright) and
is anti-conservative without `--blocks` wherever labels cluster by gene family -- the
CLI prints a warning rather than refusing, since an unblocked run is still a valid
(if weaker) null.

Neither the architecture audit (section 1) nor `--blocks` (section 2) says anything
about circularity: both are structural/QC controls, not provenance ones. A feature
axis built, even partially, from the same curated resource a label set comes from
needs the presence-leakage or provenance controls above, or the informal check
(`references/feature-axes.md`'s Circularity section) -- read what an axis is
actually built from before adding it, since nothing here stops a plain CLI config
from pairing a circular axis with its label source silently.
