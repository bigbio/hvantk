# Changelog

## 0.2.0 — 2026-08-04

First tagged release. Everything below had accumulated under `Unreleased` since `0.1.0`,
which sat on `main` unchanged from 2025-05-04 across 61 merges.

### Added

- Declarative feature selection for `hvantk rerank` (Python API: `Config.selection`). Filters run within each axis — univariate AUC with within-axis BH-FDR, then Spearman redundancy — re-fitted inside every cross-validation fold on the training slice only, so the reported ΔAUC is not inflated by selection that has seen the held-out labels. A third RFECV step is available but **off by default** (`SelectionPolicy(wrapper="rfecv")`): across four real cohorts it eliminated columns almost exclusively in the one with the fewest positives, and pruned the ablation baseline axis, so it needs an out-of-fold outcome comparison before it can be trusted by default. `Config.selection = None` (the default) reproduces the previous code path exactly, and the CLI is unchanged.
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
