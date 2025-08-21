# HVANTK Minimal Refactor Plan (Short Sprint)

## Status checkpoint (as of 2025-08-20)
- Phase 0–2: Done. Core data and annotation modules moved under hvantk.data and hvantk.annotation; imports updated.
- Phase 3: Implemented. New hvantk.tables package created; ucsc.py, expression_atlas.py, cptac.py co-located there. utils/make_tables.py moved to hvantk.tables.creators.
- Phase 4: Mostly done. hvantk/core/constants.py and hvantk/core/config.py in place; commands and internal imports updated. CLI still loads.
- Temporary deprecation shims added for stability:
  - hvantk/htables/{ucsc,expression_atlas,cptac}.py → re-export hvantk.tables.* with warnings
  - hvantk/utils/make_tables.py → re-exports hvantk.tables.creators with warning
  - hvantk/utils/constants.py, hvantk/settings.py → re-export hvantk.core.{constants,config} with warnings
- Tests updated to import hvantk.tables.* and hvantk.tables.creators where relevant; targeted test suites are passing.

Open actions (short list)
- Sweep and update any lingering imports (tests, examples, docs) to hvantk.tables.* and hvantk.core.*; avoid old htables/ and utils.make_tables paths.
- Trim shims remain minimal re-exports only (done for listed files); keep deprecation warnings.
- Define public APIs in hvantk.tables modules via __all__ (done for creators, ucsc, expression_atlas, cptac).
- Document deprecation/removal timeline for shims (see Deprecation plan below).
- Decide final location for matrix utilities and update imports (see Handling utils/matrix_utils).

## Goals and Constraints
- Timebox: a few hours, not weeks.
- Minimize blast radius: small, mechanical moves; no large redesigns.
- Temporary deprecation shims allowed to reduce churn; document removal timeline.
- Improve clarity by reducing utils/ sprawl and aligning modules with domains used today.

## Out-of-Scope (defer)
- Deep CLI regrouping or new command groups.
- AI provider modularization and plugin system.
- Complex “core” abstractions or interface layers.

## Target Shape (lightweight)
```
hvantk/
├── core/                 # small, shared basics
│   ├── constants.py
│   └── config.py         # moved from settings.py
├── data/
│   ├── data_streamer.py  # moved from utils/
│   ├── dataset.py        # moved from utils/
│   └── file_utils.py     # moved from utils/
├── annotation/
│   ├── annotate.py               # moved from utils/
│   ├── flexible_annotation.py    # moved from utils/
│   └── annotation_streamer.py    # moved from utils/
├── tables/
│   ├── creators.py              # moved from utils/make_tables.py
│   ├── expression_atlas.py      # moved from htables/
│   ├── ucsc.py                  # moved from htables/
│   ├── cptac.py                 # moved from htables/
│   └── matrix_utils.py          # moved from utils/ (see below)
├── commands/              # keep as-is for now
└── tests/                 # update imports only
```

Notes:
- Keep existing CLI commands and entry point. Only adjust imports inside them.
- Leave visualization/, resources/, examples/ untouched except for import path updates.
- ai/ stays as-is (utils/llm_interface.py) for now.

## Step-by-Step Plan (feasible in hours)

Phase 0 – Prep (30–45 min)
- Create new packages: hvantk/core, hvantk/data, hvantk/annotation, hvantk/tables, each with __init__.py.
- Decide import policy: prefer top-level from hvantk.<pkg> imports in code/tests.
- Grep for references to moved files to estimate edits.

Phase 1 – Data cleanup (45–60 min)
- Move: utils/data_streamer.py → hvantk/data/data_streamer.py
- Move: utils/dataset.py → hvantk/data/dataset.py
- Move: utils/file_utils.py → hvantk/data/file_utils.py
- Update all imports in hvantk/ and tests/ accordingly.
- Run unit tests related to data and fix import errors.

Phase 2 – Annotation consolidation (45–60 min)
- Move: utils/annotate.py → hvantk/annotation/annotate.py
- Move: utils/flexible_annotation.py → hvantk/annotation/flexible_annotation.py
- Move: utils/annotation_streamer.py → hvantk/annotation/annotation_streamer.py
- Update imports in code/tests/examples.
- Re-run tests for annotation and fix import errors.

Phase 3 – Tables co-location (45–60 min)
- Move: htables/expression_atlas.py → hvantk/tables/expression_atlas.py
- Move: htables/ucsc.py → hvantk/tables/ucsc.py
- Move: htables/cptac.py → hvantk/tables/cptac.py
- Move: utils/make_tables.py → hvantk/tables/creators.py
- Update imports in commands and tests that use table modules.
- Add temporary deprecation shims for htables/* and utils/make_tables to re-export hvantk.tables.* with warnings (remove per timeline below).

Phase 3.5 – Matrix utilities (15–30 min)
- Move: utils/matrix_utils.py → hvantk/tables/matrix_utils.py.
- Update imports in tests and any code to hvantk.tables.matrix_utils.
- Add a minimal deprecation shim at hvantk/utils/matrix_utils.py that re-exports hvantk.tables.matrix_utils and warns.
- Optional: define __all__ in hvantk/tables/matrix_utils.py.

Phase 4 – Small core cleanup (30–45 min)
- Move: utils/constants.py → hvantk/core/constants.py
- Move: hvantk/settings.py → hvantk/core/config.py
- Update hvantk/hvantk.py to import CONTEXT_SETTINGS from hvantk.core.config.
- Verify CLI still runs; adjust logging initialization if needed.

Phase 5 – Finish line (30–45 min)
- Run full test suite; fix any remaining import paths.
- Update minimal docs: README, docs/data_streamer_architecture.md references to new paths.
- Document public APIs via __all__ in hvantk.tables modules.
- Commit with a migration summary in the PR description.

## Handling utils/matrix_utils
- New location: hvantk/tables/matrix_utils.py (it operates on Hail MatrixTables and fits the tables domain).
- Steps:
  - Move the file and update imports to from hvantk.tables.matrix_utils import ... across tests and code.
  - Add hvantk/utils/matrix_utils.py shim:
    - Re-export hvantk.tables.matrix_utils (from hvantk.tables.matrix_utils import *).
    - Emit DeprecationWarning guiding users to import from hvantk.tables.matrix_utils.
  - Add __all__ to hvantk/tables/matrix_utils.py to define its public API explicitly.
  - Update docs/examples mentioning matrix utils to the new path.
- Removal timeline: remove the utils shim per the deprecation plan below.

## File Migration Map (this sprint)
- utils/data_streamer.py → hvantk/data/data_streamer.py
- utils/dataset.py → hvantk/data/dataset.py
- utils/file_utils.py → hvantk/data/file_utils.py
- utils/annotate.py → hvantk/annotation/annotate.py
- utils/flexible_annotation.py → hvantk/annotation/flexible_annotation.py
- utils/annotation_streamer.py → hvantk/annotation/annotation_streamer.py
- htables/expression_atlas.py → hvantk/tables/expression_atlas.py
- htables/ucsc.py → hvantk/tables/ucsc.py
- htables/cptac.py → hvantk/tables/cptac.py
- utils/make_tables.py → hvantk/tables/creators.py
- utils/matrix_utils.py → hvantk/tables/matrix_utils.py
- utils/constants.py → hvantk/core/constants.py
- hvantk/settings.py → hvantk/core/config.py

## Acceptance Criteria
- All imports in repo updated; no references to moved paths remain, except for temporary deprecation shims that re-export with warnings.
- hvantk CLI loads and basic commands still execute.
- Tests run green (or same status as before) after import updates.
- Documentation mentions new locations for moved modules.
- Deprecation plan is documented.

## Deprecation plan (temporary shims)
- Shims: hvantk/htables/*, hvantk/utils/make_tables.py, hvantk/utils/constants.py, hvantk/settings.py, and hvantk/utils/matrix_utils.py (once added).
- Warnings: each shim emits DeprecationWarning pointing to the new hvantk.tables.* or hvantk.core.* path.
- Timeline: keep shims for 1–2 minor releases (or ~60–90 days), then remove.
- Tracking: add a CHANGELOG entry and a short section in README about the migration paths.

## Quick Risks and Mitigations
- Hail import heaviness: keep tables code isolated under hvantk/tables/ and avoid importing it at package import time.
- Hidden cross-dependencies: if any circular import appears, prefer late imports inside functions.
- Examples may lag: prioritize updating examples that are part of CI/tests, defer the rest to a follow-up.

## Next Small Follow-ups (optional, not in this sprint)
- Group commands under subpackages (data/, annotation/, tables/) without changing UX.
- Introduce optional extras for AI providers.
- Add type checking and pre-commit hooks.
