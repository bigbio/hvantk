# Changelog

## Unreleased

### Added

- **A multiplicity correction for `hvantk rerank` (`--n-perm`).** The engine shipped the
  circularity and presence-leakage controls but nothing that asked whether a best-of-N
  delta could arise by chance, so "axis X adds +0.02" was not interpretable. Two nulls are
  built: per-axis, and the **selected maximum** over every offered axis — only the second
  is a multiplicity correction, since the maximum over several candidate axes is
  stochastically larger than any single one of them: the selected-maximum null sits above
  zero and can exceed an axis's entire measured gain, which is exactly the multiplicity a
  per-axis-only report would miss. Every permutation refits the baseline, p-values use
  `(1 + #{null >= obs}) / (1 + n_perm)` so a finite permutation set can never license
  `p = 0`, and a null refuses to answer about a delta computed under a different control
  setting. Chunkable through the API (`NullConfig(chunk=, n_chunks=)` +
  `NullDistribution.merge`) for cluster array jobs.
- **Paralogue-blocked cross-validation (`--blocks`, `--max-block-frac`).** Plain
  `StratifiedKFold` let gene families straddle folds, so pooled out-of-fold AUC was
  optimistic wherever paralogues share a label. Blocks come from each gene's first-listed
  HGNC `gene_group`, deliberately not connected components over the multi-membership
  field — that closure can collapse a large fraction of a gene universe into a single
  block and make the grouped AUC incomparable to the ungrouped one (the tell-tale is a
  *blocked* AUC scoring materially higher than the random-fold one, beyond the
  across-seed spread `--seed-sweep` reports — blocking is harder only on average, so a
  single correctly blocked run can still score above random folds by chance). A block
  over the ceiling is a hard abort, because the failure is otherwise silent.
- **A multi-seed evaluation (`--seed-sweep`).** The shipped interval resampled genes only;
  which genes landed in which fold was a second variance component fixed at one hardcoded
  seed. A single cross-validation partition can land anywhere in the across-seed spread,
  and nothing in the API let a user notice. With `--seed-sweep > 1`, the ablation now
  carries `d_lo_env`/`d_hi_env` — the union of the bootstrap interval and the across-seed
  range, never narrower than the interval alone — beside the bootstrap columns.
- **Every shipped dataset is now graded: executable contracts for `pqtl:metrics` (#414),
  `cosmic-cgc:submissions` (#417) and `alphagenome:predictions` (#421).** Each ships a
  fixture, schema and row snapshots and a Hail round-trip test, so the known-incomplete
  ledger in `hvantk/tests/test_plugin_contract_artifacts.py` is empty and all 26
  datasets run the contract. The pQTL and COSMIC fixtures are synthetic but
  format-faithful, because their licences forbid redistributing rows; AlphaGenome's is
  a real subset under the AlphaGenome Output Terms. `peptideatlas:phospho`'s fixture is
  now a synthetic raw build zip, so its round-trip test grades `parse_raw_dir` as well
  as the builder (#415).
- **`THIRD_PARTY_DATA.md` and a fixture-provenance rule (#413).** A fixture that holds
  real third-party data must have an entry there: source, version or retrieval date,
  licence, attribution and modifications (`hvantk/skills/_conventions/SKILL.md` § 9). A
  licence that forbids redistributing rows limits what a fixture may contain, not
  whether one exists. The Expression Atlas, AlphaGenome and dbNSFP fixtures are listed
  so far (#413, #421, #423).

### Changed

- **`hvantk plugins errors` now exits 1 when it lists anything** (was 0). Rows mean
  the registry is missing something, and a script asking `plugins errors` should not
  have to parse the text to learn that.
- **`hvantk drift <dataset>` now scopes load errors to the requested
  dataset/provider**, so an unrelated broken plugin no longer makes it exit 2
  (previously every load error in the registry counted). `--regenerate` exits 2
  when the dataset's own provider or entry point failed to load, or when the
  probe or the write fails (previously an uncaught traceback and exit 1), and
  rejects `--json` as a usage error.
- **`seed` is a `Config` field and a `--seed` flag.** `random_state=42` was hardcoded in
  five places across `evaluator.py`, `reranker.py` and `selection.py`; a test now fails if
  a bare `42` reappears anywhere in `hvantk/algorithms/rerank/`. Defaults are unchanged, so
  every previously produced result reproduces exactly.
- **`evaluator.py` lost its prototype header and its multi-statement lines.** The file still
  opened with a prototype path comment and packed several statements per line;
  it is the module all of the above touches most. No behaviour change — pinned by a
  golden-value test captured before the reformat.
- **`scipy` is a base dependency.** It was declared only in seven extras, yet the UCSC Cell
  Browser plugin imports `scipy.sparse` at module scope on paths a base install reaches
  (`hvantk expression summarize-ucsc`, the `ucsc-cellbrowser` builder); it only ever arrived
  through `anndata`. `requirements.txt` and `environment.yml` already installed it. The
  `cohort` extra stays, now empty, so `hvantk[cohort]` still resolves, and `poetry.lock`
  changes only its `[extras]` table (#376).
- **The `SKILL.md` spec checker requires a real body in every section.** `hvantk plugins
  validate` and the conformance test checked content in two of the nine sections (§6, §9),
  and only for keywords, so an emptied section, or a `TODO`-led stub that named the right
  words, passed. A bare `n/a` is rejected; an explained `N/A ...` is an answer and passes
  (#364).
- **`hvantk --help` describes each command the way `hvantk tools list` does.** 14 of the 18
  top-level descriptions disagreed with the tool manifests; the manifests are now canonical,
  and a test keeps `_LAZY_COMMANDS` equal to them (#304).
- **An explicit `acquisition.mode: download` must declare `lifecycle.download`.** The
  manifest schema accepted a download mode with no downloader, so the mistake surfaced only
  when someone ran `hvantk reprocess`. A downloader that is not written yet is now spelled by
  omitting the `acquisition` block (`hvantk reprocess` then needs `--skip-download`), and a
  test pins which datasets are in that state (#360).
- **Code is now formatted by `ruff format`, and `black` is dropped.** CI fails on
  unformatted files (width 88, black's default). `black` was a dev dependency that no
  workflow ran, which left 262 of 596 files unformatted. The advisory lint (never
  blocking) now also covers bugbear (`B`), blind-except (`BLE`) and bandit (`S`, minus
  `assert` in test trees), plus the preview whitespace rules (#309).
- **`poetry.lock` was regenerated with Poetry 2.3.4.** The previous 2.2.1 lock
  silently dropped the `markers` entry on eleven extras-gated packages (`cycler`,
  `fonttools`, `joblib`, `kiwisolver`, `matplotlib`, `pyparsing`, `scikit-learn`,
  `seaborn`, `threadpoolctl`, `tspex`, `xlrd`), so a base `poetry install` pulled
  all eleven in unconditionally; relocking restores the markers, so a base install
  is now eleven packages lighter. `mypy-extensions`, pulled in only by the
  now-removed `black` dev dependency (#309), drops out of the lock alongside it
  (#374). `click` moves from 8.1.8 to 8.5 in the same lock, so a lock-faithful
  install (the HPC container) gets the `Did you mean` suggestions and the click
  the CI matrix already tests.
- **`h5py` is a base dependency**, the same move `scipy` made in #376: the UCSC
  Cell Browser plugin imports it at module scope on a path a base install reaches
  (`hvantk expression summarize-ucsc`, the `ucsc-cellbrowser` builder), and it had
  arrived only transitively via `anndata` (#363).
- **`hvantk tools errors` now exits 1 when it lists anything**, matching `plugins
  errors`.
- **The `SKILL.md` checker also treats a heading-only body, an HTML comment, an
  empty code fence, a numbered-list `TODO`, or a table of `TODO`s as a
  placeholder**, extending the real-body check above to shapes that passed the
  keyword/emptiness test while carrying no content.
- **`hvantk ptm build` reports the sites mapped from each source**, as
  `Sources: UniProt (N mapped), PeptideAtlas (M mapped)`, from the new
  `PTMBuildResult.source_counts`. A source whose sites all fail to map, usually the
  wrong file, is logged as a warning naming the file.
- **Breaking (ptm). `generate_phase2_report` takes `build_result=` (a
  `PTMBuildResult`) instead of `atlas_result=`**, which now raises `TypeError`. Its
  "Atlas Summary" section is now "PTM Sites Summary", with sites per source, a
  "Mapped TSV" row (was "Combined TSV") and a Hail Table row that reads "not built"
  when no table was built.
- **`PTMBuildConfig.validate()` no longer requires `output_ht`**, which
  `ptm_build_pipeline_core` never writes, and now rejects a negative
  `flanking_codons`. `ptm_build_pipeline`, the one function that writes the table,
  checks `output_ht` (`hvantk.tools.ptm.pipeline.config_errors`) and resolves the
  `uniprot-ptm:sites` plugin before it downloads anything.
- **`alphagenome:predictions` ingests the AlphaGenome SDK's `tidy_scores` parquet instead
  of calling the AlphaGenome API (#421, #424).** It builds one row per variant, with a
  summary struct per (output type, scorer). A missing column, a missing or unknown
  scorer, a malformed variant id, or the same scores given twice (a file passed twice,
  or two scoring runs mixed) fails the build instead of being dropped or counted twice.
  A missing or NaN score leaves out only the statistic it feeds, and dropped rows are
  counted in a warning. `schema_id` is now `alphagenome-v2` and the plugin 0.2.0. The
  old builder called the API, discarded the predictions and returned only its input
  variants. Score variants with the AlphaGenome SDK (`score_variant` →
  `variant_scorers.tidy_scores()`), save the full result as parquet (the builder needs
  `track_strand` and, for splice junctions, `junction_Start`/`junction_End`), then run
  `python -m hvantk reprocess alphagenome:predictions --raw-dir <dir-with-parquet>
  --output <out.ht>`.

### Removed

- **`hvantk ptm atlas` and its Python API (`build_atlas`, `PTMAtlasConfig`,
  `PTMAtlasResult`, `DEFAULT_ATLAS_SOURCES`, module `hvantk.algorithms.ptm.atlas`)
  (#401).** The command ran only the mapping step: it never downloaded UniProt and
  never wrote the Hail Table its required `--output-ht` named, and `--flanking-codons`
  had no effect. Repaired, it would have duplicated `hvantk ptm build`, which runs the
  same mapping plus both missing steps. Migrate with
  `hvantk ptm atlas --uniprot-tsv U --peptideatlas-tsv P -o D --output-ht H` →
  `hvantk ptm build --ptm-tsv U --peptideatlas-tsv P -o D --output-ht H`. `--sources`
  has no replacement: passing a source's TSV is what selects it. `--cptac-tsv`,
  `--gtf-path`, `--overwrite` and `--flanking-codons` keep their names, and `build`
  applies `--flanking-codons` to the table with a default of 5 (`atlas` advertised 7).
  In Python, use `hvantk.tools.ptm.pipeline.ptm_build_pipeline`, or
  `hvantk.algorithms.ptm.ptm_build_pipeline_core` to map without Hail;
  `PTMAtlasConfig.uniprot_tsv` is `PTMBuildConfig.ptm_tsv`, and
  `PTMAtlasResult.combined_tsv` / `.n_sites` / `.sources_used` are
  `PTMBuildResult.mapped_tsv_path` / `.n_mapped` / `.source_counts`.
- **`hvantk.skills.alphagenome.pipelines` and `examples/alphagenome/` (#421)**, the
  pipeline that called the AlphaGenome API (`AlphaGenomePipeline`, `CheckpointManager`,
  `RateLimitedCaller`, `compute_intervals`, `load_config` and their helpers). hvantk no
  longer calls the API; ingest the SDK's `tidy_scores` output instead (see Changed).
- **An unused 11.7 MB copy of the Expression Atlas E-MTAB-6798 files** under
  `hvantk/tests/testdata/raw/expression_atlas/` (#413). The fixture the tests read is
  unchanged.

### Fixed

- **The dbNSFP drift probe could never see a new release.** It watched the legacy Google
  Sites landing page, frozen at v4.9 since 2024. It now reads the dbnsfp.org releases page,
  and the release list stays under `headers`, so a new release still opens its own
  `drift:schema` PR (#371).
- **The uniprot_ptm drift probe recorded `source_version: null` and could never
  see a UniProt release.** It now fingerprints the live response's
  `X-UniProt-Release` header and `X-Total-Results` count, failing closed when
  either is missing rather than hashing only the response's key shape, and
  rejects a result count smaller than the number of results it actually
  received (#370).
- **`hgnc:lookup` promised a column upstream no longer ships.** HGNC dropped
  `location_sortable`; the field list, fixture, snapshots, SKILL.md and drift baseline now
  follow the live dump (#355).
- **The rebuild ledger records the fingerprint bumps accepted by hand in #378
  (dbnsfp) and #381 (hgnc).** Both were accepted with a manual `--regenerate`,
  which writes no ledger row, so `hvantk drift --ledger` could not report that
  artifacts built before those bumps are pending a rebuild. The rule in the plugin
  conventions (every accepted bump gets a row in the same commit) applies to manual
  acceptances as much as to the bot's.
- **Expression Atlas transcript-level builds had non-unique `var_names`.** `var` was keyed by
  the gene id, which repeats once per transcript, and every RNA-seq accession in the shipped
  catalog is transcript-level. `var` is now keyed by the transcript id, with the gene id kept
  as the `Gene ID` column; gene-level exports are unchanged. The generic expression tools
  read `var_names` as gene ids, so for these builds they now see transcript ids, which is
  tracked in #383 (#349).
- **Two catalog URLs did not point at what their entries describe.** `gwas_catalog`'s `url`
  was the rolling `releases/latest/` alias, which has moved past the entry's release; it now
  pins that release's archive. `msigdb`'s `url` is a registration-gated landing page and is
  marked provenance-only. The closed `data_source` enum is now documented for plugin
  authors (#185).
- **The documented `hvantk reprocess` commands for `gevir` and `gwas-catalog` could not
  run.** Both declared a download mode with no downloader; they now omit the block and
  their documented commands pass `--skip-download` (the downloaders are #386).
  `ucsc-cellbrowser:adult-ctx` and `dev-ctx`, summaries derived locally from multi-gigabyte
  collections, are now `byo` (#360).
- **`hvantk hgc vds2mt` and `hgc pipeline` now write the dense MatrixTable's
  columns sorted by sample ID** (previously the VDS's own order), so the sample
  order of every exported VCF is deterministic, and may differ from that of an
  export made before this change.
  `--skip-keying-by-cols` keeps the old VDS order, and its `--help` text now says
  so. A MatrixTable written before this change should be regenerated, or combined
  through `combine_matrix_table_rows(force_sort_cols=True)`, before being unioned with one
  written after it (#368).
- **`_driver_af` breaks a tied case-carrier count toward the driver with the
  highest control-carrier frequency, deterministically.** Ties are common with
  rare variants (many genes have every driver at `cc == 1`), and the previous
  order depended on how Hail happened to collect them. `driver_af` and the
  `common_driver` audit flag can therefore change on a re-run for genes whose top
  drivers tied (#369).
- **`hvantk ptm test` fails with an actionable message naming the `constraint`
  extra when `statsmodels` is missing**, instead of a bare `No module named
  'statsmodels'` traceback that never said which extra fixes it (#362).
- **`rerank()` raises `ValueError` when `Config.nulls` is set but only the
  baseline axis has columns left in an arm** (the provenance arm restriction barred
  every other axis, or their tables contribute no feature column), rather than
  logging a warning and quietly returning `nulls=None`. A requested multiplicity
  correction can no longer go silently missing from the result.
- **`Config.selection` is now type-checked the same way as `leakage`, `nulls` and
  `blocks`, and `SelectionPolicy` rejects an unknown `univariate` / `redundancy` /
  `wrapper` value at construction.** A misspelling such as `wrapper="RFECV"`
  previously disabled that selection stage with no error or warning at all.
- **`hvantk rerank` fails with a message naming the `ml` extra when scikit-learn
  is not installed**, instead of a bare `ModuleNotFoundError` traceback.
- **`hvantk rerank --n-perm` without `--blocks` warns that the null is
  anti-conservative where labels cluster by gene family.** An unblocked null
  permutes labels across families; on family-clustered labels it can return a
  small "multiplicity-corrected" p for an axis that only recognises families. The
  warning goes to stderr, and the `--n-perm` help says the same.
- **`hvantk rerank` checks that the `--output` and `--null-out` directories exist
  before the run**, instead of failing with an `OSError` traceback after scoring
  and the whole permutation null had finished.
- **`NullDistribution.merge` refuses a chunk set with missing permutations.** Every
  chunk agreed on `planned_n_perm`, so a pre-empted array task used to shrink the
  null silently; the merge now names the missing indices, and
  `merge(..., allow_partial=True)` accepts a smaller null deliberately.
- **The permutation null's control setting records the CV seed.** The fold
  partition depends on it, yet two chunks scored under different seeds merged
  silently and a null answered for a delta computed under another seed; both now
  raise `ControlSettingMismatch`.
- **Above 10,000 rows the rerank scorer is not fully blocked, and the docs now say
  so.** `HistGradientBoostingClassifier`'s default `early_stopping="auto"` turns on
  there and holds out a random validation split that ignores the paralogue blocks;
  the estimator is unchanged so historic results reproduce.
- **The Expression Atlas builder raises `ValueError` naming the duplicated ids**
  when `var` would not be uniquely indexed — for example a transcript-level
  export whose id column is not the configured `transcript_id_column` — instead
  of writing an `.h5ad` with duplicate `var_names` that `anndata` only warns
  about.
- **The HGNC builder logs which declared fields are missing from the input
  header**, instead of silently dropping them from the renamed output.
- **Several `hvantk drift` robustness gaps closed.** `DriftResult.status` is
  validated against its four allowed values; a fingerprint is canonicalised to
  its JSON form before comparison, so a probe returning a tuple, `Path` or
  `datetime` no longer reports `drifted` forever after `--regenerate`;
  `--regenerate` refuses to write a stub or placeholder-shaped fingerprint; the
  human-readable `hvantk drift` output now prints each `probe_failed` row's
  reason to stderr; `drift --all --json` emits its load-error rows and exits 2
  when every in-scope dataset failed to bind, instead of printing nothing.
- **`hvantk ptm build` exited 0 when no PTM site mapped**, printing an empty
  `Hail Table:` line and leaving any table from an earlier run at `--output-ht` for
  `ptm annotate` and the others to read. It now fails (exit 1) and says no table was
  written.
- **`hvantk ptm build` reported its own validation errors as a crash** (`PTM build
  failed: 1`, a traceback, `Error: 1`), because the `ctx.exit` that ends it was caught
  by the command's catch-all; an empty `-o` or `--output-ht` now prints the error and
  exits 1.
- **`ptm_build_pipeline_core` checked for the UniProt TSV only after downloading and
  parsing the Ensembl GTF**, and let an empty `ptm_tsv` through; a missing or empty
  `ptm_tsv` now fails before either.
- **`cosmic-cgc:submissions` builds from the current Census export (#417).** The builder
  mapped only the legacy column names; it now also maps the 21 upper-snake-case
  columns of the current export, casts genome coordinates to `int32`, and reads
  plain-gzip input, which Hail refused to load without `force=True`. Multi-valued
  `TISSUE_TYPE` stays a known gap (#419).
- **The dbNSFP fixture is now attributed (#423).** It is real data: the header and the
  first 4,999 chromosome-10 variants of the dbNSFP v4.9a academic branch, which is
  licensed CC BY-NC-ND 4.0 (non-commercial use, attribution, no modified versions), and
  it shipped without the attribution or non-commercial notice that licence requires. It
  now has a `THIRD_PARTY_DATA.md` entry and a `NOTICE.md` beside the fixture and the
  snapshot, and the catalog entry names the licence.
- **The documented COSMIC CGC build command crashed (#424).** `hvantk reprocess
  cosmic-cgc:submissions --raw-dir <dir>` handed the directory to a builder that reads
  one file, and failed with `IsADirectoryError`. The skill page and the data-sources
  guide now give the form that works, `--raw-dir <dir> --intermediate <dir>/<file>
  --skip-parse --skip-download`, as for `gwas-catalog:associations`.
- **PeptideAtlas per-site `n_observations` counted each peptide once per modified form
  (#425, #430).** `parse_peptideatlas_zip` added the parent peptide's total count, which
  also covers its unmodified and other forms, once for every phospho form instead of
  each form's own count. On the real build 202512/606 that inflated 237,428 of 259,932
  sites (82× in total, median 16× per site) and reordered them (Spearman 0.82 against the
  corrected counts). Rebuild any `peptideatlas:phospho` table built before this fix, and
  re-check anything that used its counts. A phospho form without an integer count now
  fails the parse.

## 0.3.1 — 2026-08-30

Reworks the scheduled drift bot. Fewer PRs, each carrying a signal that means
something, and each visible enough to get reviewed.

In its first month the bot opened **34 PRs from 7 of 26 datasets**, all looking alike.
On 2026-08-28 five were found sitting unreviewed for three days and were then merged
without CI having ever run on them.

### Added

- **A rebuild ledger.** `hvantk/resources/drift_ledger.json` records each accepted
  upstream change per dataset. `hvantk drift --ledger` lists datasets whose upstream
  moved since their last rebuild; `hvantk drift --mark-rebuilt <dataset>` clears one.
  A fingerprint bump is not only a test-baseline update — it is also the signal that a
  built artifact may be stale (ClinVar gained ~408 KB of variants across 2026-08), and
  that fact previously survived only in git history.
- **`informational`, a fingerprint block excluded from drift comparison.** Probes can
  record context for human readers — upstream publish dates and similar — without it
  becoming a drift trigger. Sits alongside `fetched_at` and `probe_version` in
  `PROBE_FINGERPRINT_IGNORED_KEYS`.
- **Risk classification and batching.** Drifted datasets are split into *routine*
  (schema signal unchanged) and *schema* (column list or header hash moved).
  All routine ones share one branch and become one PR; each schema change keeps its
  own. Classification defaults to *schema* for anything it cannot read — burying a
  builder-breaking change inside a batch is worse than one extra PR.
- **A weekly probe-health workflow** (`drift-health.yml`) that opens no PRs and only
  files an issue when a probe reports `probe_failed`. Exists because moving
  regeneration to fortnightly doubles the worst-case delay on a *broken* probe, and
  the `cptac` probe had failed silently for months before anyone noticed.
- **A fast path-filtered CI job** (`drift-validate.yml`) gating fingerprint-only PRs in
  ~1m40s instead of the ~25-minute three-version matrix.
- **A monthly promotion workflow** (`drift-promote.yml`) that opens one `dev` → `main`
  PR when the delta is fingerprints-only. It never merges.
- **Stale-PR escalation.** A drift PR open past a full cycle gets one comment on itself
  — never a second PR, which is the failure #262 fixed for force-pushes.

### Changed

- **`hgnc` and `gencc` probes record `Content-Length` as their content signal**, and
  `Last-Modified` moves to `informational`. Both previously carried `Last-Modified` in
  `source_version`, which drift compares, while hashing only the column-header line —
  a *schema* signal. So they recorded no content signal at all, and a byte-identical
  re-publish was indistinguishable from a real update. Across all 8 committed hgnc
  fingerprints from 2026-05-16 to 2026-08-27 the checksum never moved while
  `Last-Modified` moved every time. Both probes now fail closed if the server omits
  `Content-Length`, rather than recording a fingerprint with no content signal — the
  scheduled bot regenerates drifted baselines automatically, so one transient omission
  would otherwise be baked in permanently. `PROBE_VERSION` → 2 for both; baselines
  regenerated. **Anything parsing those files by hand needs updating.**
- **The drift bot authenticates as a GitHub App** rather than `GITHUB_TOKEN`. A PR
  authored by `GITHUB_TOKEN` has its workflow runs parked at `action_required` until a
  human approves them, so no drift PR had ever been CI-tested unattended.
- **Drift PRs open ready for review, never as drafts**, labelled `drift:routine` or
  `drift:schema`, assigned from the plugin's `maintainers:` or the
  `DRIFT_DEFAULT_ASSIGNEE` environment variable, with a classification table leading
  the body ahead of the raw JSON diffs.
- **Regeneration runs fortnightly** (1st and 15th) rather than daily. Nothing upstream
  moves faster than weekly in a way that matters. Note GitHub schedules cron
  best-effort under load — the day is reliable, the hour is not.
- **`hvantk drift --ledger` rejects** being combined with `--all`, a dataset argument,
  `--regenerate` or `--json`, rather than silently ignoring them.

### Fixed

- **A failed mid-batch regeneration no longer contaminates an unrelated PR.** Earlier
  datasets' regenerated files were left in the working tree; the next handler's
  `git checkout -B` carried them onto its branch and `git add hvantk/skills` committed
  them — invisible in that PR's diff, body and ledger, while stderr claimed no PR had
  been opened for them.
- **The promotion gate can fire at all.** It excluded only `drift_fingerprint.json`,
  but every drift commit also writes `drift_ledger.json`, so from the first merged
  drift PR onward it would have skipped permanently, green and silent.
- **Stale-PR escalation is idempotent.** It was a pure function of the PR's creation
  time with no memory, so a PR open six months would have collected ~11 identical
  comments.
- **`load_ledger` no longer raises on truthy non-dict JSON**, honouring its documented
  contract; and ledger staleness compares timestamps as instants rather than strings,
  so a `-05:00` offset no longer reads as older than a `+00:00` one.

## 0.3.0 — 2026-08-05

### Added

- **`--n-partitions` now controls the VDS → MatrixTable partitioning** on both
  `hvantk hgc vds2mt` and `hvantk hgc pipeline`, coalescing the dense MatrixTable before
  the write. A VDS's on-disk layout is derived from its *reference-block* count, which is
  a property of the genome and saturates (on the order of 10^8 on chr1 by a few hundred
  samples) while the dense matrix keeps growing with N×M(N). Past that point the partition
  count stops tracking the size of the data it partitions and work-per-task collapses —
  measured on a ~1,000-sample chr1 cohort, where partitions shrank to under a MiB each and
  densify and QC sped up less than 2x going from 16 to 128 cores (8x) while well-sized
  stages scaled nearly linearly (#207).
  Implemented with `naive_coalesce`, which merges adjacent partitions without a shuffle so
  the densify for a merged group runs inside one task. Reduces only; the default is
  unchanged. **Not** implemented at the read: `hl.vds.read_vds(n_partitions=…)` looks
  tidier but derives intervals from the reference data via `_calculate_new_partitions`,
  whose count saturates independently of the request — measured on the 2,586-partition
  test VDS it returned 2 intervals for every request from 2 to 100, and `to_dense_mt` then
  failed a Scala `require` on the reference/variant mismatch for requests of 2, 4 and 16,
  where coalescing returned exactly 2, 4 and 16.

### Changed

- **Datasets that share a `drift_fingerprint` baseline are now treated as sharing one
  drift signal**, rather than as N independent ones. `hvantk drift` probes such a group
  once and fans the result out — every dataset still gets its own report entry, and each
  now carries `fingerprint_path` — and the drift workflow opens a single PR per signal.
  `ucsc-cellbrowser` is the case that forced it: `default`, `adult-ctx` and `dev-ctx` are
  genuinely distinct *schema* variants (their obs cell-type column is `celltype`, `Class`
  and `Type_v2`, which is why each earns its own snapshot), but `fetch_fingerprint()`
  takes no arguments and fingerprints the provider-wide catalog at
  `cells.ucsc.edu/dataset.json`. One upstream event therefore produced three identical
  PRs whose branches all wrote the same file, so merging any one made the other two
  conflict — #241 merged, #242 and #243 were closed as superseded. Grouping is keyed on
  *(baseline path, probe callable)*, not the path alone: two datasets sharing a baseline
  while declaring different probes would each overwrite the other's, so they deliberately
  do not group, and `hvantk plugins validate` now rejects that declaration outright. A
  lone dataset keeps its historical `drift/<provider>-<dataset>` branch name exactly, so
  existing open PRs are still matched; a group uses `drift/<provider>`.

### Fixed

- **`PipelineConfig.n_partitions` reached no pipeline stage.** It was accepted from the
  CLI and echoed back in the run plan while being read by nothing, so the run plan
  affirmatively told the user a setting had taken effect when it had not — the worst
  failure mode for a dead flag, and the first knob a user reaches for when they hit #207.
  It is now forwarded to `convert_vds_to_mt`; the run plan line names the stage it governs
  (#208).
- **Every open drift PR was force-pushed and its body re-edited once a day, forever.**
  Nine PRs churned daily for a week — roughly 63 notification events, none carrying new
  information. `hvantk drift --regenerate` rewrites `fetched_at` on every run, and the
  existing emptiness check compared against the *base branch*, so for a dataset that was
  still drifted it could never fire: the timestamp alone guaranteed a non-empty diff.
  Nothing compared the freshly regenerated fingerprint against what the branch already
  proposed. `drift_to_pr.py` now skips the push and the PR edit when the branch already
  carries a materially identical fingerprint — "materially" meaning equal once
  `fetched_at` and `probe_version` are dropped, the same keys the drift detector ignores.
  The skip additionally requires an **open PR** to still exist: a branch outlives its PR
  when one is closed, and a run whose push succeeded while `gh pr create` failed leaves a
  branch with no PR at all — in both cases the branch content matches, so a content-only
  check would suppress that dataset's drift forever. Every other outcome (no branch yet,
  a file the branch lacks, an unreadable blob, a git failure) still pushes: suppressing a
  real drift PR is far worse than one redundant force-push. A skipped dataset also
  restores the index and working tree before returning — `git checkout -B` does not clear
  the index, so a leftover staged fingerprint would be committed onto the *next*
  dataset's branch — and still reports itself in the job's step summary. Note this does **not** reduce how often drift is *detected* — a content
  revision, such as ClinGen's `content_length` moving while the checksum holds, is still
  a genuine change and still opens a PR.

- **The `cptac:expression` and `cptac:phospho` drift probes had never once succeeded.**
  The drift workflow installed with a bare `pip install -e .`, but the cptac probe
  fingerprints the *installed* `cptac` version (via `importlib.metadata`) against
  PayneLab's latest GitHub release — and `cptac` is declared in the `ptm` extra. Every
  scheduled run reported `The 'cptac' Python package is not installed; cannot
  fingerprint`, so upstream CPTAC drift has never been detectable. The workflow now
  installs `.[ptm]`. Of the 22 drift probes these two are the only ones needing an
  extra; the other 20 use `requests` or the stdlib alone.
- **A single rate-limited response could fail the whole scheduled drift run.** On
  2026-08-04 GenCC answered the regeneration request with HTTP 429; the probe had no
  retry, so it raised, no PR could be opened for the drifted dataset, and the run exited
  non-zero with nothing actually wrong. New `hvantk/core/utils/http.py` provides
  `request_with_retry`, which retries transient statuses (429 and the 5xx family) and
  connection errors with exponential backoff. It honours `Retry-After` but **clamps**
  it: `urllib3.util.Retry` sleeps for the header's full value with no upper bound
  (`backoff_max` caps only the exponential path), so a host answering `Retry-After: 3600`
  would park CI for an hour. The helper deliberately does not call `raise_for_status`,
  so callers keep their existing error handling and only the transient case changes.
  Wired into the GenCC probe; the other 12 HTTP probes can adopt it as needed.

## 0.2.0 — 2026-08-04

First tagged release. Everything below had accumulated under `Unreleased` since `0.1.0`,
which sat on `main` unchanged from 2025-05-04 across 61 merges.

### Added

- Declarative feature selection for `hvantk rerank` (Python API: `Config.selection`). Filters run within each axis — univariate AUC with within-axis BH-FDR, then Spearman redundancy — re-fitted inside every cross-validation fold on the training slice only, so the reported ΔAUC is not inflated by selection that has seen the held-out labels. A third RFECV step is available but **off by default** (`SelectionPolicy(wrapper="rfecv")`): it eliminated columns mostly where positives were fewest and can prune the ablation baseline axis, so it needs an out-of-fold outcome comparison before it can be trusted by default. `Config.selection = None` (the default) reproduces the previous code path exactly, and the CLI is unchanged.
- `rerank_arms(config)` runs each analysis as two arms, `clean` and `all`, over identical folds. `clean` (columns with no provenance conflict against the label source) is the headline; `all` adds conflicted and undeclared columns so the circularity channel is a measured number rather than an assumption. `RerankResult.selection` carries the per-fold selection frequency, the global-pass feature list, and both nested and global AUCs.
- Plugin manifests may declare per-predictor training provenance: an optional `scores: {<column>: {trained_on: [...]}}` block per dataset. `hvantk/skills/dbnsfp/plugin.yaml` declares it for 55 of its 57 rankscore predictors. An omitted score means unknown and is never treated as clean.
- Plugin system for data-provider adapters. Each provider now lives in a single folder under `hvantk/skills/<provider>/` with a `plugin.yaml` manifest, builder code, drift probe, downloader CLI, and tests. The loader auto-discovers plugins from the in-tree filesystem and Python entry points.
- `hvantk plugins {list,describe,errors,validate}` commands for inspecting the registry.
- `hvantk drift <provider:dataset>` for upstream-drift detection against committed expected fingerprints.
- `hvantk reprocess <provider:dataset>` for chaining download -> parse -> build -> drift-check from a single command.
- 13 migrated provider plugins: clingen (gene-disease), clinvar, cptac (expression + phospho), expression-atlas, gencc (submissions), gtex-eqtl, gwas-catalog, hgnc, insider, msigdb, peptideatlas (phospho), ucsc-cellbrowser (default / adult-ctx / dev-ctx), uniprot-ptm (sites).
- Scheduled CI workflow (`.github/workflows/drift.yml`) that runs `hvantk drift --all --json` daily and opens a draft PR per drifted plugin with the regenerated fingerprint pre-committed.

### Changed

- **Version bumped to `0.2.0`, and `pyproject.toml` migrated to PEP 621 `[project]`.** The
  version had been `0.1.0` since 2025-05-04, across 61 merges into `main` — so no release in
  fifteen months was distinguishable from any other by version. Separately, `name`,
  `version`, `description`, `authors`, `license`, `readme`,
  `keywords`, `urls`, `plugins`, `extras` and `scripts` all used the deprecated
  `[tool.poetry.*]` spelling — 11 warnings on every `poetry check`. They now live under
  `[project]`, `[project.optional-dependencies]`, `[project.entry-points]`,
  `[project.scripts]` and `[project.urls]`; only genuinely Poetry-specific keys
  (`include`/`exclude`, dependency groups) remain under `[tool.poetry]`.
  **The migration is resolution-neutral**: the lock resolves to the same 187 packages,
  name-for-name and version-for-version, before and after. Poetry's caret shorthand is
  spelled out as the PEP 508 equivalent it always meant (`^8.1.3` → `>=8.1.3,<9.0.0`), not
  re-pinned. With nothing deprecated left, the defensive `poetry>=2.0,<3.0` pin in the
  `poetry.lock in sync` CI job is unpinned again.
  One consequence is new: PEP 621 puts the full specifier in each extra, so `scipy>=1.8`
  is written six times and `scikit-learn>=1.4,<2.0` three times, and one could be re-pinned
  with the others left behind — resolving differently depending on which extra a user
  installs. `test_pyproject_extras.py` now asserts every extra spells a shared package
  identically, alongside a check that no extra re-declares a base dependency.
- **`gnomad` is no longer a dependency.** hvantk used exactly one function from it,
  `annotate_adj`, which is ~15 lines of Hail expression with no gnomAD data behind it.
  It is now ported into `hvantk/algorithms/hgc/adj.py` (gnomad_methods is MIT; the port
  keeps the logic and thresholds verbatim and carries the attribution), so `adj` means
  exactly what it means in a gnomAD callset. Dropping the dependency removes **35
  packages** from the lock — `hgvs`, `ga4gh-vrs`, `onnx`, `onnxruntime`, `skl2onnx`,
  `psycopg2`, `protobuf`, `sympy`, `slackclient` and more — and, critically, removes the
  transitive `jsonschema<4` pin that conflicted with hvantk's own declared
  `jsonschema>=4.0`. That conflict is what had made `poetry.lock` impossible to
  regenerate in place. `adjust_genotypes=True` no longer requires an optional install,
  so the `hgc` extra is now just `["matplotlib", "seaborn"]`.
- **`poetry.lock` regenerated and now consistent with `pyproject.toml`.** It had drifted
  across ~28 commits — pinning `jsonschema` 3.2.0 against a declared `>=4.0`, and missing
  the `ancestry`/`ml`/`constraint`/`expression` extras entirely — so `poetry install`
  failed on a clean checkout. 208 → 181 packages; the only version change besides the
  removals is `jsonschema` 3.2.0 → 4.26.0. `scanpy`'s move behind the `expression` extra
  is now actually in effect rather than merely declared.
- **Breaking (rerank).** `ArmAssignment.unknown` is renamed `undeclared`, and an undeclared
  predictor is now treated as *conflicted* rather than getting a bucket of its own. Arm
  membership is otherwise unchanged. Code reading `ArmAssignment.unknown` must be updated.
- **Rebuild your dbNSFP artifact.** `dbnsfp:variants` now parses the ~57 `*_rankscore`
  columns to `float64` with proper missingness, instead of leaving them as raw strings
  (`"."` for missing). `schema_id` stays `dbnsfp-v1` — the column *set* is unchanged and
  string rankscores were always a parsing bug rather than an intended schema — so nothing
  will warn you: an artifact built before this release carries strings where a fresh build
  carries floats. Re-run `hvantk reprocess dbnsfp:variants` before relying on those columns.
- Specificity features in the annotation matrix now emit a per-group **vector** by default
  (one column per surviving atlas group, named `{atlas}_{sanitized_group}`) instead of a
  single rolled-up scalar. A named roll-up is *additive* when `specificity.targets` is
  given; `specificity.emit: rollup` restores the previous single-column output. Two matrix
  axes may no longer share an `atlas` label, since their vector columns would collide.
- `hvantk drift` comparators can now actually detect an upstream change. Previously the
  comparison could pass regardless of source content, so drift went unreported.
- `scikit-learn` floor raised to `>=1.4` (NaN-tolerant tree estimators, needed by rerank's
  optional RFECV wrapper). `scipy` is now declared explicitly, as an optional dependency in
  the `ml` / `ancestry` / `psroc` extras.
- **Every command that imports `scipy` at module scope now has an extra that installs it.**
  Three modules do, and they sit on three *different* commands — a mapping the previous
  known-gaps note got wrong:
  `algorithms/ptm/constraint.py` → `hvantk ptm constraint` (`constraint` extra);
  `algorithms/enrichex/overlap.py` → `hvantk enrichex overlap` (new `enrichex` extra);
  `algorithms/burden/fet.py` → `hvantk cohort burden` (new `cohort` extra).
  `fet.py` is *not* reached by `hvantk enrichex burden`: its only importer is
  `algorithms/burden/pipeline.py`, imported by `tools/cohort/cohort_cli.py` alone.
  Previously `constraint` omitted `scipy`, so `pip install hvantk[constraint]` yielded a
  documented extra whose own command still raised `ModuleNotFoundError`, and neither
  `hvantk enrichex overlap` nor `hvantk cohort burden` had any extra to install.
  `enrichex` also carries matplotlib/seaborn, because `enrichex/__init__` imports
  `plot.py`/`report.py` unconditionally and a scipy-only extra would break on import.
  `scipy` was added to `expression` too — a resolution no-op, since scanpy already depends
  on it, but `visualization/expression/anndata.py` imports `scipy.sparse` directly and this
  project declares what it imports rather than inheriting it from a transitive edge that can
  move. A base install still raises a bare `ModuleNotFoundError` rather than a message
  naming the extra; a `require_scipy()` guard (cf. `require_scanpy`) would fix the *text*,
  and is tracked separately because it changes no extra's contents.
- **`hvantk` was unusable from a `pip install`.** The console script imported
  `hvantk.tools.enrichex` at module scope, which ran `algorithms/enrichex/__init__.py`,
  which eagerly imported `enrichex/plot.py`, `enrichex/report.py` and
  `visualization/base.py` — all three import `matplotlib` at module scope. matplotlib is
  optional, so on a base install **every** command including `hvantk --help` raised
  `ModuleNotFoundError`. Those three imports are now resolved on attribute access (PEP 562
  `__getattr__`), so the package imports without matplotlib and the CLI runs. The eight
  plotting/reporting names stay in `__all__` and stay importable; touching one without
  matplotlib now raises an `ImportError` naming the `enrichex` extra, matching
  `require_scanpy` and `_require_matplotlib`. No plotting behaviour changed — the enrichex
  CLIs already imported `generate_report` inside the functions that use it.
- **The wheel shipped 43.3 MB of test data.** `hvantk/tests/**` was absent from `exclude`
  (190 files, including a 14.7 MB VDS zip and an 11.4 MB expression-atlas fixture); the
  skills excludes were overridden by `include = "hvantk/skills/**/*.py"`, since a path named
  by `include` wins; and the excludes named `tests/data/**` where the skills actually use
  `tests/testdata/**`. Fixed all three: the wheel goes from **46.2 MB to 2.86 MB**
  uncompressed (28 MB to 897 KB on disk) with every manifest, skills module, catalog and
  drift fingerprint intact.
- **CI now installs the package.** New `packaging-smoke` job builds the wheel, checks its
  contents with `.github/scripts/check_wheel.py`, installs it into a clean environment with
  no extras, and runs `hvantk --help` / `hvantk plugins list` from a directory where the
  checkout is not importable — so the console script, the entry-point registrations and the
  packaging globs are exercised against the installed copy. It also asserts the provider
  count matches the tree, since a dropped manifest would otherwise still exit 0. Both bugs
  above were found by writing this job.
- **`hvantk/tests/hgc/` now runs in CI** as a new `hgc-hail` job — separate from
  `Plugin contract (hail)` rather than appended to it, so the contract signal is not delayed
  behind ~6 min of unrelated HGC work. This is the first automatic run of
  `test_convert_vds_to_mt`, the end-to-end exercise of the `adj` code ported in #252.
- **The Python version matrix now runs on `dev` PRs**, not only `main`, so an incompatibility
  is caught on one commit instead of at the release gate with a whole release to bisect.
  `actions/setup-python` moved v3 → v5 and both workflows now declare
  `permissions: contents: read` (both raised in review on #249).
- The extras table is now guarded by a test. It is duplicated in three places — the
  `[tool.poetry.extras]` block, `README.md` and `docs_site/getting-started/installation.md`
  — and only the first is executable, so the prose copies had drifted eight cells
  (`psroc`/`ancestry`/`ml` missing `scipy`, `ptm` missing `sorted-nearest`) across two
  releases. `hvantk/tests/test_pyproject_extras.py` now parses both markdown tables and
  fails if either disagrees with `pyproject.toml`, and also fails if an extra names a
  package that is not declared `optional = true`. Non-Hail, so it runs in the default suite.
- `scanpy` moved out of the base install into a new `expression` extra. It is required by
  `hvantk expression summarize`, `hvantk expression markers`, and `hvantk ptm constraint
  --expression-metric mean`; those now fail with an actionable message naming the extra
  rather than a bare `ModuleNotFoundError`. The extra cannot be installed on Intel macOS
  (scanpy → numba → llvmlite ships no x86_64 macOS wheel from 0.47).
- Package restructured into 4 purpose-driven roofs: `core/` (platform models, utilities, plugin/tool runtime, streamers, transient builders), `algorithms/` (analytical computation: ptm, psroc, qtlcascade, enrichex, hgc, ancestry, annotation, visualization, expression, statistics, training_sets), `skills/` (data ingestion plugins), `tools/` (CLI surface). Inside `core/` there are now sub-packages `models/`, `utils/`, `streamers/`, `plugin/`, `tool/`, `builders/` so adding a new format helper has one obvious home. One-way dependency rule (`skills/`, `tools/` → `algorithms/` → `core/`) is enforced by `hvantk/tests/test_dependency_directions.py`. `hvantk/data/`, `hvantk/utils/`, `hvantk/tables/`, and 8 top-level algorithm dirs (`hvantk/{ptm,psroc,qtlcascade,enrichex,hgc,ancestry,annotation,visualization}/`) are gone. `ClinVarStreamer` no longer imports from `hvantk.skills.clinvar.builder` — it accepts a pre-built Hail Table via its constructor.
- Registry keys for migrated providers use compound `provider:dataset` form. Recipe JSONs and any custom callers should update from bare names (e.g., `clinvar`) to compound (`clinvar:variants`). The legacy `hvantk mktable` / `hvantk mkmatrix` CLI surfaces have been retired; data builds now go through `hvantk reprocess <provider>:<dataset>` with `--plugin-arg key=value` for builder kwargs.
- Plugin manifests gain an optional `catalog: <path>` field pointing at a per-plugin `catalog/datasets.json`. `unified_registry.HvantkRegistry` now aggregates per-plugin catalogs from the plugin loader in addition to the legacy `resources/registry/genomics/datasets.json`.
- Per-domain catalogs `resources/registry/{transcriptomics,proteomics,epigenomics}/datasets.json` are removed; their entries now live inside each owning plugin's `catalog/datasets.json` (expression-atlas, ucsc-cellbrowser). `registry/genomics/datasets.json` is intentionally retained until orphan entries (dbNSFP, gnomad-metrics, ensembl-gene, gevir, cosmic-cgc) gain owning plugins.
- `hvantk catalog` CLI rewritten to read per-plugin catalogs via `HvantkRegistry`. New subcommands: `list` (with `--omics-type` / `--data-source` / `--organism` filters), `show`, `stats`, `search`. The legacy `catalog build` subcommand is removed; use `hvantk reprocess <provider:dataset>` instead.

### Removed

- Per-provider downloader modules under `hvantk/commands/*_downloader.py` for migrated providers (moved into their plugin folder's `cli.py`).
- Per-provider dataset classes under `hvantk/datasets/*_datasets.py` for migrated providers (moved into `hvantk/skills/<provider>/shared/`).
- Per-provider builder functions in `hvantk/tables/table_builders.py` and `matrix_builders.py` for migrated providers (moved into `hvantk/skills/<provider>/[<dataset>/]builder.py`).
- `hvantk/resources/generate_catalog.py` (regenerated the now-removed per-domain `datasets.json` files). Catalog regeneration is now a per-plugin concern; if a maintainer needs a packaged regenerator in the future it should live alongside each plugin's `catalog/datasets.json`.
- `hvantk/resources/catalog.yaml` (auto-generated summary file pointing at deleted per-domain JSON files). Equivalent information is available on demand via `hvantk catalog stats`.
- `openai`, `anthropic`, `google-genai` and `RestrictedPython` dropped from
  `requirements.txt` and `environment.yml`. None is imported anywhere in the tree, and
  none was ever declared in `pyproject.toml` — CI had been installing four packages the
  library does not use.

### Known gaps before first stable release

- The following plugins reference snapshot files (`schema.json`, `sample_rows.json`)
  in their `plugin.yaml` manifests that have not yet been seeded on disk:
  `clingen`, `gencc`, `hgnc`, `uniprot-ptm`, `expression-atlas`, `peptideatlas:phospho`,
  `cptac:expression`, and `cptac:phospho`. The first hail-enabled CI run with
  `--regenerate-snapshots` will bootstrap them. All `ucsc-cellbrowser` variants
  (`default`, `adult-ctx`, `dev-ctx`) already have populated snapshot dirs.
- (Resolved: version bumped to 0.2.0 — see Changed.) Releases are still not git-tagged, so
  a release is identifiable by version but not by a tag.
  (The three CI gaps previously listed here — no install job, the version matrix running
  only on `main`, and `hvantk/tests/hgc/` running in no job — are resolved; see the
  packaging and CI entries under Changed.)
