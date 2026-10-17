# Inputs and identifiers

`hvantk rerank` joins three inputs -- the cohort table, every feature axis, and the
labels -- purely by exact string match on a column literally named `gene`. There is no
identifier mapping anywhere inside the tool (`hvantk/algorithms/rerank/engine.py`,
`hvantk/algorithms/rerank/features.py`): all three must already share one identifier
space before the config is written.

## How the three inputs actually join

- **Features.** `FeatureAssembler.assemble` (`features.py`) outer-joins every feature
  axis table on its `gene` column. This outer join IS the ranked universe -- nothing
  else determines it.
- **Cohort / prior.** `engine.rerank()` loads the cohort manifest's prior column via
  `CohortManifest.key_column` (renamed to `gene`; `hvantk/algorithms/cohort/frame.py`
  `load_prior_frame`) and left-joins it onto the feature universe
  (`matrix.merge(prior, on="gene", how="left")`). The feature universe is the LEFT
  side, so a cohort gene absent from every feature axis never enters the scored
  table at all -- it is not an error, it simply never appears. A feature-axis gene
  with no matching cohort row keeps `prior_stat` as `NaN` (it is still ranked).
- **Labels.** `y = df["gene"].isin(pos)` -- a plain set-membership test against the
  label set, no join at all. A label naming a gene outside the universe still
  COUNTS in `min_label_coverage`'s denominator -- it is not simply ignored, it
  actively lowers the ratio (`engine.py`: `covered = sum(1 for g in pos if g in
  set(matrix["gene"])); covered / len(pos)`, where `len(pos)` is the FULL label
  count). Below the floor, `rerank()` raises the exact message:
  `label/feature join coverage too low: only {covered}/{len(pos)}
  ({covered/len(pos):.0%}) positive units are in the feature matrix — likely a
  gene-symbol/build mismatch`.
- **`--blocks`** (`hvantk/algorithms/rerank/blocks.py` `load_gene_groups`) matches the
  ranked table's `gene` values against an HGNC complete-set TSV's approved `symbol`
  column ONLY -- not `alias_symbol`, not `prev_symbol`, not an Ensembl id. If none of
  the universe's genes appear in that column, `gene_blocks` raises
  (`none of the N gene(s) in this universe appear in the gene-group table`); if the
  match rate is under 50% it logs a warning and treats the unmatched genes as
  unconstrained singletons rather than failing.

The practical consequence: whichever identifier space you choose, `--blocks` forces
that choice to be HGNC approved symbols, because blocking needs the gene-group table
and the gene-group table only speaks symbols. A cohort or feature axis keyed by
Ensembl `gene_id` must be mapped to symbol before `--blocks` can do anything.

## Mapping Ensembl `gene_id` to an approved symbol

The gene index (`hvantk annotate spine`, see the `gene-index` skill) is keyed
`gene_id` and carries `hgnc_id` as a left join -- a protein-coding gene with no HGNC
record keeps a missing `hgnc_id` rather than being dropped. `prepare-rerank-inputs.py`
resolves `gene_id -> hgnc_id -> symbol` through that gene index plus the raw HGNC
complete-set TSV (`hgnc_id`, `symbol`, `status` columns -- the raw download, not the
built `hgnc:lookup` table, whose own builder renames `symbol` to `gene_symbol`
internally; see `hvantk/skills/hgnc/shared/constants.py` `HGNC_GENE_FIELDS`).
Approved-only: a withdrawn or merged HGNC entry's symbol is not a statement about
whichever gene carries that symbol today, the same reasoning
`hvantk/algorithms/rerank/blocks.py`'s own `load_gene_groups` applies.

Two things are dropped, and both are counted and printed, never silently:

- **No approved symbol.** A `gene_id` with no `hgnc_id`, or an `hgnc_id` the HGNC
  table does not carry as `Approved`.
- **A symbol shared by more than one `gene_id`.** Rare, but mapping both sides onto
  the same `gene` value would either silently merge two genes' rows or (worse) raise
  later as a duplicate-key error deep inside `load_cohort_frame`
  (`cohort {name} table has N duplicate row(s) for gene key(s) ...`). Both `gene_id`s
  are dropped from the map, not just one, so the drop is symmetric and the count is
  exact.

Apply this SAME mapping to the cohort table and to every feature axis that needs it
(`prepare-rerank-inputs.py cohort` / `features`) so a `gene_id` that drops out of one
drops out of the other identically -- mapping each input independently, from
different runs of the same lookup, can only ever agree or disagree by construction
luck, not by anything the tool enforces. Both subcommands refuse outright (exit 1)
if the mapped share of their own input's rows falls below `--min-mapping-rate`
(default 0.9, the floor `hvantk cohort attach` already uses for the same kind of
join -- `hvantk.algorithms.annotation.mapping.enforce_rate`'s
`{source} mapped {rate:.2f} of its identifiers, below the required {min_rate:.2f}`),
so a wholesale id-space mismatch is caught here rather than surfacing later as
`hvantk rerank` silently ranking a near-empty, all-`NaN`-prior universe.

Both subcommands also SORT their output by the mapped symbol before writing
(`references/feature-axes.md`'s determinism note: with one feature axis, which a
small run commonly has, nothing downstream re-sorts the universe, so the axis
file's own row order becomes the cross-validation fold assignment).

`cohort --manifest-out` writes the new manifest's `table:` field as an ABSOLUTE
path and refuses a relative `--out` outright: `hvantk` resolves a manifest's
`table:` against whichever process reads the manifest later, never against the
manifest file's own location (`hvantk/algorithms/cohort/frame.py`), and that
reading process is commonly a Slurm job whose working directory is its own log
directory -- a relative path would resolve against the wrong base every time.

## Coverage, symbols and duplicates

- `min_label_coverage` (default 0.5) is the fraction of label-positive units that
  must land in the feature matrix. Restrict the LABEL file to the universe
  (`prepare-rerank-inputs.py labels --restrict-to ... --restrict-column ...`) rather
  than lowering this number -- a low coverage is almost always a space mismatch or a
  universe that is narrower than the label source, and the number exists to catch
  exactly that.
- `FeatureAssembler.assemble` drops duplicate `gene` keys WITHIN one axis table (kept
  first, logged), to avoid a cartesian blow-up on the outer join -- this is silent
  unless you read the log; a mapping step that already guarantees one row per gene
  (as above) never reaches this path.
- `load_cohort_frame` (the cohort table's own loader) refuses outright if the mapped
  table ends up with more than one row for the same `gene` key, or a blank/`NA`-token
  key on any row.

## Labels: from a GeneSetCollection to a flat list

`hvantk genesets clingen`/`gencc` write a `GeneSetCollection` JSON -- several NAMED
groups (by GCEP, by submitter, by disease, ...), never a flat gene list.
`prepare-rerank-inputs.py labels` reads it with the repository's own reader
(`hvantk.core.utils.gene_sets.GeneSetCollection.load`, pure Python, no Hail) and
unions either every group or only the `--group` names given, then optionally
restricts to the universe, then writes one approved symbol per line, sorted. The
genes `ClinGenGeneDiseaseTableStreamer` collects are already `gene_symbol` values from
the built table itself -- no further HGNC mapping is applied to labels by this
skill, only the restriction to the universe.
