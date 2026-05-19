# hvantk Package Restructure — Design

- **Date:** 2026-05-19
- **Branch:** continues on `feat/data-handlers-refactoring`
- **Status:** Approved for implementation planning
- **Owner:** Yasset Perez-Riverol
- **Builds on:** the plugin/tool restructure already on this PR (~50 commits ending at `e04ba61`)

## 1. Problem

After the plugin and tool restructures landed, `hvantk/` has 14 top-level subpackages with no clear architectural roles:

```
ancestry/  annotation/  core/  data/  enrichex/  hgc/  psroc/
ptm/  qtlcascade/  resources/  skills/  tables/  tests/
tools/  utils/  visualization/  + hvantk.py
```

Three concrete problems:

1. **`hvantk/data/` is a junk drawer.** It mixes pure utilities (`file_utils.py`, `gene_mapper.py`) with provider-specific streamers (`clinvar_streamer.py`, `gencc_streamer.py`, …). The line is fuzzy and contributors don't know where new code belongs.
2. **`hvantk/core/` mixes contracts with workflow utilities.** Stable plugin/tool contracts (`plugin_api.py`, `tool_api.py`, manifest schemas) sit alongside format helpers (`anndata_utils.py`, `bgzf.py`, `converters.py`). Adding a new format helper looks the same as evolving a contract.
3. **Algorithms have no roof.** `ptm/`, `psroc/`, `qtlcascade/`, `enrichex/`, `hgc/`, `ancestry/`, `annotation/`, `visualization/` are all "analytical computation," but they sit at the top level alongside infrastructure packages, blurring the user's mental model.

## 2. Goals

1. **Three packages by purpose**, plus tools and resources: `core/` (platform models + utilities + plugin runtime), `algorithms/` (analytical computation), `skills/` (data ingestion plugins). `tools/` (CLI surface) and `resources/` (catalog data) sit alongside.
2. **One-way dependency flow** so the packages don't tangle. Enforced by a CI-runnable test, not just convention.
3. **Internal organization inside `core/`** via sub-packages (`core/models/`, `core/utils/`, `core/streamers/`, `core/plugin/`, `core/tool/`, `core/builders/`) so a new format helper has one obvious home.
4. **Preserve blame** for every moved file via `git mv`.
5. **No public-API regressions** for users who already imported from the old paths during the transition. Deprecation shims bridge old import paths until Phase 7 removes them.

## 3. Non-goals

- Authoring per-algorithm manifests (`algorithm.yaml`) parallel to plugin.yaml / tool.yaml. Worth doing later so `hvantk algorithms list` becomes a real command, but explicitly out of scope here.
- Plugin-izing the 5 residual builders that land in `core/builders/` (dbNSFP, gnomad-metrics, gevir, cosmic-cgc, ensembl-gene). Each is its own follow-up PR.
- Moving algorithm tests inside their algorithm folders. They stay under top-level `hvantk/tests/<algo>/` for now to keep diff size manageable.
- Renaming the test-harness directories (`tests/testdata/`, `tests/snapshots/`).

## 4. Decisions

| Question | Decision |
|---|---|
| Top-level layout | `core/`, `algorithms/`, `skills/`, `tools/`, `resources/`, `tests/`, `hvantk.py` |
| `core/` internal structure | Sub-packages: `models/`, `utils/`, `streamers/`, `plugin/`, `tool/`, `builders/`, plus `protocols.py`, `constants.py`, `config.py` at the top |
| Provider-specific streamers (clinvar, clingen, gencc, cosmic_cgc, alphagenome) | `core/streamers/` — they're consumed by multiple downstream callers (annotation, hgc, utils), not by their owning plugin |
| Streamer → plugin import (today's `core → skills` violation) | Fixed via dependency injection: streamers accept the built table as a parameter; callers provide it |
| Algorithm tests | Stay at `hvantk/tests/<algo>/` for this PR; algorithm-folder co-location deferred |
| Dependency-direction enforcement | New `hvantk/tests/test_dependency_directions.py` asserts no upward imports |
| Migration scope | All 14 top-level subpackages re-homed in one PR (12 phased commits) |
| Restructure on top of | The existing `feat/data-handlers-refactoring` branch (already ~50 commits) |

## 5. Dependency direction

The whole structure works only if dependencies flow one direction:

```
skills/      ──┐
               ├──→  algorithms/  ──→  core/
tools/       ──┘
```

| From | May import from |
|---|---|
| `core/<subpkg>` | other `core/` sub-packages (with caveats below) |
| `algorithms/<algo>` | `core/**` |
| `skills/<provider>` | `core/**` |
| `tools/<domain>/<tool>` | `core/**`, `algorithms/**`, `skills/**` |

**Forbidden:**

- `core/**` may NOT import from `algorithms/`, `skills/`, or `tools/`.
- `algorithms/**` may NOT import from `skills/` or `tools/`.
- `skills/**` may NOT import from `algorithms/`, `tools/`, or sibling `skills/`.
- `tools/**` is the only layer allowed to depend on everything, because the CLI is where workflows compose.

### Within `core/` sub-packages

```
core/protocols  ──┐
core/constants  ──┤
core/config     ──┤
                  ├──→  core/models   ──→  core/utils   ──→  core/streamers
core/models     ──┤                                          ↓
core/utils      ──┤                                       core/plugin
core/streamers  ──┤                                       core/tool
core/builders   ──┘
```

`core/utils/` only imports from `core/{protocols,constants,config}`. `core/streamers/` may import from `core/{models,utils}`. `core/plugin/` and `core/tool/` may import from any other `core/` sub-package (the loaders resolve callables across them).

### Enforcement

A new `hvantk/tests/test_dependency_directions.py` walks every `.py` file under each layer and asserts no forbidden imports. Caught at CI time, not runtime. Initially has `pytest.mark.xfail` markers because the rules are violated everywhere; each migration phase removes its xfails as the violations get fixed.

### The streamer-decoupling fix

Today `hvantk/data/clinvar_streamer.py` imports from `hvantk.skills.clinvar.builder`. After this restructure, `clinvar_streamer` lives in `core/streamers/` (because annotation, hgc, utils all consume it). That import — `core → skills` — violates the rule.

The streamer is refactored to accept the built Hail Table as a parameter instead of building it internally. Callers (annotation/clinvar-builder/etc.) provide the table at construction time. This is the cleanest decoupling and matches how `algorithms/annotation/` will compose builders + streamers anyway.

## 6. Package layout

### `core/` (reorganized in place)

```
hvantk/core/
  protocols.py
  constants.py
  config.py

  models/                   # platform data-model definitions
    metadata.py             ← from core/metadata.py
    dataset.py              ← from data/dataset.py
    anndata_utils.py        ← from core/anndata_utils.py
    backends.py             ← from core/backends.py

  utils/                    # format-agnostic / cross-cutting helpers
    file_utils.py           ← from data/file_utils.py
    bgzf.py                 ← from core/bgzf.py
    gene_mapper.py          ← from data/gene_mapper.py
    gene_aliases.py         ← from utils/gene_aliases.py
    gene_sets.py            ← from utils/gene_sets.py
    geneset_io.py           ← from utils/geneset_io.py
    table_utils.py          ← from utils/table_utils.py
    expressions.py          ← from utils/expressions.py
    genome.py               ← from utils/genome.py
    obo_parser.py           ← from utils/obo_parser.py
    mondo_parser.py         ← from utils/mondo_parser.py
    catalog.py              ← from utils/catalog.py
    hail_context.py         ← from core/hail_context.py
    converters.py           ← from core/converters.py
    readers.py              ← from core/readers.py
    writers.py              ← from core/writers.py
    router.py               ← from core/router.py

  streamers/                # cross-cutting streamers (consumed by algorithms/)
    base.py                 ← from data/data_streamer.py
    gene_disease.py         ← from data/gene_disease_streamer.py
    clinvar.py              ← from data/clinvar_streamer.py
    clingen.py              ← from data/clingen_streamer.py
    gencc.py                ← from data/gencc_streamer.py
    cosmic_cgc.py           ← from data/cosmic_cgc_streamer.py
    alphagenome.py          ← from data/alphagenome_streamer.py

  plugin/                   # plugin runtime
    api.py                  ← from core/plugin_api.py
    loader.py               ← from core/plugin_loader.py
    manifest.schema.json    ← from core/plugin_manifest.schema.json
    drift_runner.py         ← from core/drift_runner.py
    registry.py             ← from tables/registry.py

  tool/                     # tool runtime
    api.py                  ← from core/tool_api.py
    loader.py               ← from core/tool_loader.py
    manifest.schema.json    ← from core/tool_manifest.schema.json

  builders/                 # TRANSIENT — unmigrated builders, become plugins next
    table.py                ← from tables/table_builders.py
    matrix.py               ← from tables/matrix_builders.py
    genome.py               ← from tables/genome_builders.py
```

### `algorithms/` (new top-level package)

```
hvantk/algorithms/
  ptm/             ← from hvantk/ptm/
  psroc/           ← from hvantk/psroc/
  qtlcascade/      ← from hvantk/qtlcascade/
  enrichex/        ← from hvantk/enrichex/
  hgc/             ← from hvantk/hgc/
  ancestry/        ← from hvantk/ancestry/
  annotation/      ← from hvantk/annotation/
  visualization/   ← from hvantk/visualization/
  expression/      ← matrix_utils.py, tissue_specificity.py from utils/
  statistics/      ← wilcoxon.py, correction.py from utils/
  training_sets/   ← generate_training_set_streamed.py,
                     generate_enhanced_training_set.py from utils/
```

### Untouched

- `hvantk/skills/` — plugin layout already final (13 providers, 16 datasets)
- `hvantk/tools/` — tool layout already final (24 tools across 11 domains)
- `hvantk/resources/` — catalog + JSON schemas
- `hvantk/tests/` — cross-cutting + algorithm tests stay here
- `hvantk/hvantk.py` — Click entry point (imports get rewritten in Phase 6)

### Deleted (after content moves)

```
hvantk/data/              ← removed
hvantk/utils/             ← removed
hvantk/tables/            ← removed
hvantk/ancestry/          ← moved
hvantk/annotation/        ← moved
hvantk/enrichex/          ← moved
hvantk/hgc/               ← moved
hvantk/psroc/             ← moved
hvantk/ptm/               ← moved
hvantk/qtlcascade/        ← moved
hvantk/visualization/     ← moved
```

### Per-file judgment calls during implementation

Two files have shape that requires reading their contents to place correctly:

- `hvantk/utils/catalog.py` and `hvantk/utils/genome.py`: default placement is `core/utils/`. If either turns out to be algorithm-shaped (operating on loaded data, not just helping with IO), the implementer relocates it to `algorithms/` during the migration commit.
- `hvantk/tables/genome_builders.py` (`build_1k_genome`): default placement is `core/builders/`. If the team prefers treating it as an algorithm (it processes user-supplied raw data), it can move to `algorithms/build/`. Implementer flags this during Phase 3.

## 7. Migration plan (12 commits, 8 phases)

### Phase 0 — Scaffolding (1 commit)

- Create empty sub-package skeletons: `core/{models,utils,streamers,plugin,tool,builders}/__init__.py` and `algorithms/{ptm,psroc,qtlcascade,enrichex,hgc,ancestry,annotation,visualization,expression,statistics,training_sets}/__init__.py`.
- Add `hvantk/tests/test_dependency_directions.py` with `pytest.mark.xfail` markers (the rules are violated everywhere; phases remove xfails as they fix the violations).
- Commit: `chore(restructure): scaffold core sub-packages + algorithms package`

### Phase 1 — `core/` internal reorganization (1 commit)

- `git mv` all 17 files currently in `hvantk/core/*.py` into their sub-package homes.
- Add re-export shims to `hvantk/core/__init__.py` for the prior import paths (`from hvantk.core import plugin_loader`, etc.) with `DeprecationWarning`.
- Commit: `refactor(core): reorganize into models/utils/plugin/tool sub-packages`

### Phase 2 — `data/` → `core/{models,utils,streamers}/` (1 commit)

- `git mv` the 11 files.
- Apply the streamer-decoupling fix: `clinvar_streamer` (and any other streamer that imports from `skills/`) accepts the built table as a parameter.
- Update the ~25 callers (annotation, hgc, plugins, tests).
- Remove `hvantk/data/__init__.py` re-exports.
- Commit: `refactor(streamers): move data/ into core/{models,utils,streamers}/, decouple from skills/`

### Phase 3 — `tables/` → `core/plugin/registry.py` + `core/builders/` (1 commit)

- `tables/registry.py` → `core/plugin/registry.py` (plugin runtime).
- `tables/table_builders.py`, `matrix_builders.py`, `genome_builders.py` → `core/builders/`.
- Update any callers (mostly tool wrappers in `tools/build/`).
- Commit: `refactor(builders): move tables/ into core/plugin/registry + core/builders/`

### Phase 4 — `utils/` split (1 commit)

- Generic utilities (`gene_aliases`, `gene_sets`, `geneset_io`, `table_utils`, `expressions`, `genome`, `obo_parser`, `mondo_parser`, `catalog`) → `core/utils/`.
- Computation-flavored (`matrix_utils`, `tissue_specificity`, `wilcoxon`, `correction`, training-set generators) → respective `algorithms/<sub>/`.
- Commit: `refactor(utils): split into core/utils/ + algorithms/{expression,statistics,training_sets}/`

### Phase 5 — `algorithms/` moves (8 sub-commits, one per algorithm)

- `hvantk/<algo>/` → `algorithms/<algo>/` for ptm, psroc, qtlcascade, enrichex, hgc, ancestry, annotation, visualization.
- Each commit updates the algorithm's own imports + any callers (mostly `tools/<domain>/`).
- Commits: `refactor(algorithms/ptm): move from hvantk/ptm/` (and 7 more)

### Phase 6 — Tool import-path updates (1 commit)

- `hvantk/tools/<domain>/*.py` and top-level `hvantk/hvantk.py` get imports rewritten.
- Tool manifests (`<basename>.tool.yaml`) get `cli.module:` fields updated.
- Commit: `refactor(tools): rewire imports for the new core/ + algorithms/ layout`

### Phase 7 — Remove deprecation shims + enforce direction (1 commit)

- Drop the `hvantk/core/__init__.py` re-export shims added in Phase 1.
- Remove `xfail` from `test_dependency_directions.py`.
- Delete any now-empty package directories.
- Commit: `refactor(restructure): remove deprecation shims, enforce dependency direction`

### Phase 8 — Documentation (1 commit)

- Update `hvantk/skills/_conventions/SKILL.md` and `docs_site/architecture.md` with the new layout.
- CHANGELOG entry describes the four-package shape and the dependency-direction rule.
- Commit: `docs(restructure): document core/algorithms/skills/tools architecture`

### Total churn

- ~50 file moves (`git mv` for blame preservation)
- ~400 import-path updates (mechanical, batched per phase)
- 1 new test file (`test_dependency_directions.py`)
- 12 commits on the same branch

## 8. Risks

1. **Circular imports during streamer decoupling (Phase 2).** If any caller relies on the streamer to build its own table, the decoupling requires refactoring the caller. Mitigation: spot-check the 5 callers (annotation_streamer, hgc, generate_training_set_streamed, PSROC, util) BEFORE Phase 2 lands; refactor any caller that requires the streamer to build its own table in a precursor commit.

2. **Test import path churn.** Tests reach into every layer. Mechanical sed updates work, but `unittest.mock.patch("hvantk.data.clinvar_streamer.ClinVarStreamer")` strings get easily missed. Mitigation: after each phase, run `pytest --collect-only` to surface import-resolution failures before actual test execution.

3. **Deprecation shims (Phase 1) accidentally become permanent.** Easy to forget the Phase 7 cleanup. Mitigation: each shim's docstring includes `# REMOVE IN PHASE 7` as a grep target.

4. **`build_1k_genome` reclassification.** Currently in `hvantk/tables/genome_builders.py`. Default placement is `core/builders/`. If the team prefers `algorithms/build/`, the implementer flags during Phase 3.

5. **Algorithm tests don't move.** They stay under `hvantk/tests/<algo>/`. This is deliberate but creates a layout asymmetry with plugins (which co-locate tests). Future PR can move them if the asymmetry is bothersome.

## 9. Testing strategy

- **Per-phase baseline check**: before each phase, capture `pytest hvantk/ -m "" --no-header | tail -3`. After each phase, re-run; pass count must match or increase.
- **`test_dependency_directions.py`** is the load-bearing acceptance check. Once it's green (Phase 7), the architecture is enforced for the future.
- **CLI smoke after Phase 7**: `hvantk plugins list`, `hvantk tools list`, `hvantk catalog stats` must produce the same output as today.

## 10. Out of scope for this branch

- Per-algorithm manifests (`algorithm.yaml`) parallel to plugin.yaml / tool.yaml. Worth a follow-up so `hvantk algorithms list` becomes real.
- Plugin-izing the 5 residual builders in `core/builders/` (dbNSFP, gnomad-metrics, gevir, cosmic-cgc, ensembl-gene). Each is its own follow-up PR.
- Moving algorithm tests inside their algorithm folders.
- Renaming test-harness directories.
- A new top-level `model/` (singular) package; the user's mental model used "core/models" which is captured as `core/models/`.
