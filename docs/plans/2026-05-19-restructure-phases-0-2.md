# Package Restructure — Phases 0–2 Implementation Plan

> **For agentic workers:** Execute this plan task-by-task — write each failing test, run it red, implement, run it green, commit. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the foundation of the package restructure spec: scaffold the new `core/` sub-packages and `algorithms/` package, reorganize `core/`'s internal layout, dissolve `hvantk/data/` into `core/{models,utils,streamers}/`, and apply the streamer-decoupling fix so streamers no longer import from `skills/`.

**Architecture:** Four-package layout (`core/` + `algorithms/` + `skills/` + `tools/`) with one-way dependency flow enforced by a CI test. `core/` is sub-divided internally (`models/`, `utils/`, `streamers/`, `plugin/`, `tool/`, `builders/`). After this plan, `data/` is gone and `core/` has its final shape; `tables/`, `utils/`, and the algorithm packages still need migration (covered by follow-up plans).

**Tech Stack:** Python 3.11+, pytest, `git mv` for blame preservation, Poetry. No new dependencies.

**Spec reference:** [docs/specs/2026-05-19-package-restructure-design.md](docs/specs/2026-05-19-package-restructure-design.md). This plan implements §§5–7 (Phases 0–2). Phases 3–8 are follow-up plans.

---

## File Structure (created or modified in this plan)

| Path | Responsibility |
|---|---|
| `hvantk/core/models/__init__.py` | New sub-package marker |
| `hvantk/core/utils/__init__.py` | New sub-package marker |
| `hvantk/core/streamers/__init__.py` | New sub-package marker |
| `hvantk/core/plugin/__init__.py` | New sub-package marker |
| `hvantk/core/tool/__init__.py` | New sub-package marker |
| `hvantk/core/builders/__init__.py` | New sub-package marker (placeholder; populated in Phase 3) |
| `hvantk/algorithms/__init__.py` | New top-level package marker |
| `hvantk/algorithms/{ptm,psroc,qtlcascade,enrichex,hgc,ancestry,annotation,visualization,expression,statistics,training_sets}/__init__.py` | Empty package markers for follow-up plans |
| `hvantk/tests/test_dependency_directions.py` | NEW — asserts no upward imports; uses xfail for layers not yet fixed |
| `hvantk/core/__init__.py` | Re-export shims for backwards compatibility (Phase 1) |
| `hvantk/core/{plugin_api,plugin_loader,plugin_manifest.schema,drift_runner}.py` → `hvantk/core/plugin/` | Moved via `git mv` |
| `hvantk/core/{tool_api,tool_loader,tool_manifest.schema}.py` → `hvantk/core/tool/` | Moved |
| `hvantk/core/{anndata_utils,backends,metadata}.py` → `hvantk/core/models/` | Moved |
| `hvantk/core/{bgzf,converters,readers,writers,router,hail_context}.py` → `hvantk/core/utils/` | Moved |
| `hvantk/data/file_utils.py` → `hvantk/core/utils/file_utils.py` | Moved (Phase 2) |
| `hvantk/data/gene_mapper.py` → `hvantk/core/utils/gene_mapper.py` | Moved |
| `hvantk/data/dataset.py` → `hvantk/core/models/dataset.py` | Moved |
| `hvantk/data/data_streamer.py` → `hvantk/core/streamers/base.py` | Moved |
| `hvantk/data/gene_disease_streamer.py` → `hvantk/core/streamers/gene_disease.py` | Moved |
| `hvantk/data/clinvar_streamer.py` → `hvantk/core/streamers/clinvar.py` | Moved + decoupled from skills/ |
| `hvantk/data/{clingen,gencc,cosmic_cgc,alphagenome}_streamer.py` → `hvantk/core/streamers/<name>.py` | Moved |
| `hvantk/data/__init__.py` | Deleted at the end of Phase 2 |
| ~25 caller files | Imports updated (concrete paths listed inside each task) |

---

## Phase 0 — Scaffolding

### Task 1: Add the dependency-direction test (xfail-gated)

**Files:**
- Create: `hvantk/tests/test_dependency_directions.py`

- [ ] **Step 1: Write the test**

```python
# hvantk/tests/test_dependency_directions.py
"""Assert that hvantk's four-package layout (core / algorithms / skills / tools)
honors the one-way dependency rule documented in the design spec:

    skills/      ──┐
                   ├──→  algorithms/  ──→  core/
    tools/       ──┘

Implementation note: each layer-pair is checked independently and gated with
xfail until the corresponding migration phase fixes the violations. Phases
remove xfail markers as they land. Phase 7 deletes every xfail marker —
after that, this test enforces the contract for the future.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _imports_in(layer: str) -> list[tuple[Path, str]]:
    """Return (file, dotted-import) pairs for every import statement in <layer>."""
    out: list[tuple[Path, str]] = []
    root = REPO_ROOT / layer
    if not root.is_dir():
        return out
    for py in root.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        try:
            tree = ast.parse(py.read_text())
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    out.append((py, alias.name))
            elif isinstance(node, ast.ImportFrom) and node.module:
                out.append((py, node.module))
    return out


def _forbidden_matches(layer: str, forbidden_prefixes: list[str]) -> list[tuple[Path, str]]:
    bad: list[tuple[Path, str]] = []
    for file, dotted in _imports_in(layer):
        for prefix in forbidden_prefixes:
            if dotted == prefix or dotted.startswith(prefix + "."):
                bad.append((file, dotted))
    return bad


@pytest.mark.xfail(reason="enforced from Phase 7 onward", strict=False)
def test_core_does_not_import_upward():
    bad = _forbidden_matches(
        "core", ["hvantk.algorithms", "hvantk.skills", "hvantk.tools"]
    )
    assert not bad, (
        "core/ must not import from algorithms/skills/tools. Offenders:\n"
        + "\n".join(f"  {p.relative_to(REPO_ROOT)} -> {d}" for p, d in bad)
    )


@pytest.mark.xfail(reason="enforced once algorithms/ migration lands (Phase 5)", strict=False)
def test_algorithms_does_not_import_skills_or_tools():
    bad = _forbidden_matches("algorithms", ["hvantk.skills", "hvantk.tools"])
    assert not bad, (
        "algorithms/ must not import from skills/ or tools/. Offenders:\n"
        + "\n".join(f"  {p.relative_to(REPO_ROOT)} -> {d}" for p, d in bad)
    )


def test_skills_does_not_import_algorithms_or_tools():
    """Skills are siblings — they meet only through core/. Already true today
    (Phase 1 of the original plugin migration enforced this for the 13
    migrated plugins). Should stay green; no xfail."""
    bad = _forbidden_matches("skills", ["hvantk.algorithms", "hvantk.tools"])
    sibling_bad = [
        (file, dotted)
        for file, dotted in _imports_in("skills")
        if dotted.startswith("hvantk.skills.")
        and dotted.split(".")[2] != file.relative_to(REPO_ROOT / "skills").parts[0]
    ]
    assert not bad and not sibling_bad, (
        "skills/ must not import from algorithms/, tools/, or sibling skills/. "
        f"Offenders: {bad + sibling_bad}"
    )
```

- [ ] **Step 2: Run test to verify the xfail behavior**

Run: `/usr/local/bin/pytest hvantk/tests/test_dependency_directions.py -v -m ""`
Expected:
- `test_core_does_not_import_upward` → XFAIL (today's `core/data/` streamers import from `skills/`)
- `test_algorithms_does_not_import_skills_or_tools` → XFAIL (no `algorithms/` exists yet)
- `test_skills_does_not_import_algorithms_or_tools` → PASS (already enforced)

If `test_skills_...` fails: that's a real regression from the prior migration. STOP and report.

- [ ] **Step 3: Commit**

```bash
git add hvantk/tests/test_dependency_directions.py
git commit -m "test(deps): add dependency-direction test (xfail until phase 7)"
```

---

### Task 2: Scaffold `algorithms/` package skeleton

**Files:**
- Create: `hvantk/algorithms/__init__.py`
- Create: `hvantk/algorithms/{ancestry,annotation,enrichex,expression,hgc,ptm,psroc,qtlcascade,statistics,training_sets,visualization}/__init__.py`

- [ ] **Step 1: Create the directories and `__init__.py` files**

```bash
mkdir -p hvantk/algorithms
touch hvantk/algorithms/__init__.py
for sub in ancestry annotation enrichex expression hgc ptm psroc qtlcascade statistics training_sets visualization; do
  mkdir -p "hvantk/algorithms/$sub"
  touch "hvantk/algorithms/$sub/__init__.py"
done
```

Each `__init__.py` should be empty (the package marker only). Verify with:

```bash
find hvantk/algorithms -type f -name '__init__.py' | sort
```

Expected: 12 lines, one per directory (root + 11 sub-packages).

- [ ] **Step 2: Verify nothing else broke**

Run: `/usr/local/bin/pytest hvantk/tests/test_dependency_directions.py -v -m ""`
Expected: same outputs as Task 1 step 2 (the empty `algorithms/` directory has no imports, so the algorithms xfail is still XFAIL because it asserts "no algorithms/ imports from skills/tools" — empty case is vacuously true, will now PASS, not XFAIL).

If the second test (`test_algorithms_does_not_import_skills_or_tools`) flips from XFAIL to XPASS, that's fine — it'll stay green from now on. (Phase 5 fills algorithms/ with real code; the test continues to gate against violations.)

- [ ] **Step 3: Commit**

```bash
git add hvantk/algorithms/
git commit -m "chore(restructure): scaffold algorithms/ package + 11 sub-packages"
```

---

## Phase 1 — `core/` internal reorganization

### Task 3: Move plugin-runtime files into `core/plugin/`

**Files moved:**
- `hvantk/core/plugin_api.py` → `hvantk/core/plugin/api.py`
- `hvantk/core/plugin_loader.py` → `hvantk/core/plugin/loader.py`
- `hvantk/core/plugin_manifest.schema.json` → `hvantk/core/plugin/manifest.schema.json`
- `hvantk/core/drift_runner.py` → `hvantk/core/plugin/drift_runner.py`

**Files modified** (imports + schema-load paths):
- `hvantk/core/plugin/loader.py` (after move): its `_SCHEMA_PATH = Path(__file__).parent / "plugin_manifest.schema.json"` becomes `Path(__file__).parent / "manifest.schema.json"`.
- `hvantk/tools/plugins/{plugins_cli,drift_cli,reprocess_cli,tools_cli}.py` and `hvantk/resources/unified_registry.py`: any `from hvantk.core import plugin_loader` stays valid via the shim added in Task 7, so no immediate update needed.

- [ ] **Step 1: Create the `core/plugin/` package**

```bash
mkdir -p hvantk/core/plugin
touch hvantk/core/plugin/__init__.py
```

- [ ] **Step 2: Move the 4 files via `git mv`**

```bash
git mv hvantk/core/plugin_api.py hvantk/core/plugin/api.py
git mv hvantk/core/plugin_loader.py hvantk/core/plugin/loader.py
git mv hvantk/core/plugin_manifest.schema.json hvantk/core/plugin/manifest.schema.json
git mv hvantk/core/drift_runner.py hvantk/core/plugin/drift_runner.py
```

- [ ] **Step 3: Update the schema path inside `core/plugin/loader.py`**

The file currently has:
```python
_SCHEMA_PATH = Path(__file__).parent / "plugin_manifest.schema.json"
```

Since `__file__` is now `hvantk/core/plugin/loader.py`, `Path(__file__).parent` already points at `hvantk/core/plugin/` which contains `manifest.schema.json` (after the rename above). Change the filename in the path:

```python
_SCHEMA_PATH = Path(__file__).parent / "manifest.schema.json"
```

- [ ] **Step 4: Update internal imports inside the moved files**

`core/plugin/loader.py` originally had `from .plugin_api import ...`. The `plugin_api` module is now `core/plugin/api.py` — change to `from .api import ...`:

```python
from .api import (
    DatasetSpec,
    PluginLoadError,
    PluginNameCollision,
    Provider,
    TestPaths,
)
```

`core/plugin/drift_runner.py` originally had `from .plugin_api import DatasetSpec, DriftProbeError`. Change:

```python
from .api import DatasetSpec, DriftProbeError
```

It also had `from . import plugin_loader` inside `run_drift_check`. Change:

```python
from . import loader as plugin_loader
```

(Or just `from . import loader` and reference `loader.get_registry()`. Pick whichever reads cleaner; the import is local to a function.)

- [ ] **Step 5: Run plugin-system tests to confirm no regression**

Many call-sites still use `from hvantk.core import plugin_loader` — that import path is broken until Task 7 adds the shim. So at this point, expect failures. That's OK; we're in the middle of the phase. Run anyway to confirm the failures are import-path-only:

Run: `/usr/local/bin/pytest hvantk/tests/test_plugin_api.py -v -m "" 2>&1 | tail -10`

Expected: collection or import errors mentioning `cannot import name 'plugin_api'` or similar. Do NOT commit yet — the next tasks fix this.

- [ ] **Step 6: Do NOT commit yet**

Task 7 adds the back-compat shim. Single commit covers Tasks 3–7.

---

### Task 4: Move tool-runtime files into `core/tool/`

**Files moved:**
- `hvantk/core/tool_api.py` → `hvantk/core/tool/api.py`
- `hvantk/core/tool_loader.py` → `hvantk/core/tool/loader.py`
- `hvantk/core/tool_manifest.schema.json` → `hvantk/core/tool/manifest.schema.json`

- [ ] **Step 1: Create the package + move files**

```bash
mkdir -p hvantk/core/tool
touch hvantk/core/tool/__init__.py
git mv hvantk/core/tool_api.py hvantk/core/tool/api.py
git mv hvantk/core/tool_loader.py hvantk/core/tool/loader.py
git mv hvantk/core/tool_manifest.schema.json hvantk/core/tool/manifest.schema.json
```

- [ ] **Step 2: Update the schema path inside `core/tool/loader.py`**

Change:
```python
_SCHEMA_PATH = Path(__file__).parent / "tool_manifest.schema.json"
```
to:
```python
_SCHEMA_PATH = Path(__file__).parent / "manifest.schema.json"
```

- [ ] **Step 3: Update internal imports inside the moved file**

`core/tool/loader.py` originally had `from .tool_api import ...`. Change:

```python
from .api import (
    Subcommand,
    ToolLoadError,
    ToolRequirements,
    ToolSpec,
)
```

- [ ] **Step 4: Do NOT commit yet — Task 7 commits the whole `core/` reorg as one unit**

---

### Task 5: Move model files into `core/models/`

**Files moved:**
- `hvantk/core/anndata_utils.py` → `hvantk/core/models/anndata_utils.py`
- `hvantk/core/backends.py` → `hvantk/core/models/backends.py`
- `hvantk/core/metadata.py` → `hvantk/core/models/metadata.py`

- [ ] **Step 1: Create + move**

```bash
mkdir -p hvantk/core/models
touch hvantk/core/models/__init__.py
git mv hvantk/core/anndata_utils.py hvantk/core/models/anndata_utils.py
git mv hvantk/core/backends.py hvantk/core/models/backends.py
git mv hvantk/core/metadata.py hvantk/core/models/metadata.py
```

- [ ] **Step 2: Check for intra-`core/` imports between these files**

Run: `grep -n 'from hvantk.core\|from \.' hvantk/core/models/*.py`

If `metadata.py` imports from `anndata_utils.py` via `from .anndata_utils import ...` — that relative import still works (same package). No change needed.

If any of these files import from `hvantk.core.backends` directly (absolute), change to `from hvantk.core.models.backends`. Same for the other two file names.

- [ ] **Step 3: Do NOT commit yet**

---

### Task 6: Move utility files into `core/utils/`

**Files moved:**
- `hvantk/core/bgzf.py` → `hvantk/core/utils/bgzf.py`
- `hvantk/core/converters.py` → `hvantk/core/utils/converters.py`
- `hvantk/core/readers.py` → `hvantk/core/utils/readers.py`
- `hvantk/core/writers.py` → `hvantk/core/utils/writers.py`
- `hvantk/core/router.py` → `hvantk/core/utils/router.py`
- `hvantk/core/hail_context.py` → `hvantk/core/utils/hail_context.py`

**Files NOT moved (stay at `core/` top level):**
- `hvantk/core/__init__.py`
- `hvantk/core/config.py`
- `hvantk/core/constants.py`
- `hvantk/core/protocols.py`

- [ ] **Step 1: Create + move**

```bash
mkdir -p hvantk/core/utils
touch hvantk/core/utils/__init__.py
git mv hvantk/core/bgzf.py hvantk/core/utils/bgzf.py
git mv hvantk/core/converters.py hvantk/core/utils/converters.py
git mv hvantk/core/readers.py hvantk/core/utils/readers.py
git mv hvantk/core/writers.py hvantk/core/utils/writers.py
git mv hvantk/core/router.py hvantk/core/utils/router.py
git mv hvantk/core/hail_context.py hvantk/core/utils/hail_context.py
```

- [ ] **Step 2: Check for intra-`core/` imports**

Run: `grep -rn 'from hvantk.core\|from \.\.\?' hvantk/core/utils/`

For each match:
- Relative imports between sibling utils files (`from .bgzf import ...`) still work — no change.
- Absolute imports like `from hvantk.core.bgzf import ...` need to become `from hvantk.core.utils.bgzf import ...`.

Common patterns to fix:
- `converters.py` may import from `bgzf.py` and `hail_context.py` — relative imports work.
- `readers.py` / `writers.py` may import from `anndata_utils` (now in `models/`) — change to `from hvantk.core.models.anndata_utils import ...`.

After editing, verify each utility module can be imported standalone (will hail-fail in this env, but the syntax + import resolution should be clean):

```bash
python3 -c "import ast; ast.parse(open('hvantk/core/utils/converters.py').read()); print('ok')"
```

Repeat for the other 5 files.

- [ ] **Step 3: Do NOT commit yet — Task 7 commits the whole `core/` reorg**

---

### Task 7: Add backwards-compat shims to `hvantk/core/__init__.py`, then commit

**Files modified:**
- `hvantk/core/__init__.py` — add lazy re-exports for the prior import paths

- [ ] **Step 1: Read the current `__init__.py`**

Check what's already there:
```bash
cat hvantk/core/__init__.py
```

Whatever is there, preserve it. Append the shims below.

- [ ] **Step 2: Add re-export shims with DeprecationWarning**

Append to `hvantk/core/__init__.py`:

```python
# ---------------------------------------------------------------------------
# Backwards-compat shims for callers using the pre-Phase-1 import paths.
# REMOVE IN PHASE 7 of the restructure plan.
# ---------------------------------------------------------------------------
import warnings as _warnings


def __getattr__(name: str):
    """Forward old `from hvantk.core import plugin_loader` style imports to
    their new sub-package locations, with a DeprecationWarning so callers
    can find the import to update.
    """
    _shim_map = {
        # plugin runtime
        "plugin_api":              "hvantk.core.plugin.api",
        "plugin_loader":           "hvantk.core.plugin.loader",
        "drift_runner":            "hvantk.core.plugin.drift_runner",
        # tool runtime
        "tool_api":                "hvantk.core.tool.api",
        "tool_loader":             "hvantk.core.tool.loader",
        # models
        "anndata_utils":           "hvantk.core.models.anndata_utils",
        "backends":                "hvantk.core.models.backends",
        "metadata":                "hvantk.core.models.metadata",
        # utils
        "bgzf":                    "hvantk.core.utils.bgzf",
        "converters":              "hvantk.core.utils.converters",
        "readers":                 "hvantk.core.utils.readers",
        "writers":                 "hvantk.core.utils.writers",
        "router":                  "hvantk.core.utils.router",
        "hail_context":            "hvantk.core.utils.hail_context",
    }
    if name in _shim_map:
        _warnings.warn(
            f"hvantk.core.{name} moved to {_shim_map[name]}; "
            "update the import (this shim is removed in restructure Phase 7)",
            DeprecationWarning,
            stacklevel=2,
        )
        import importlib
        return importlib.import_module(_shim_map[name])
    raise AttributeError(f"module 'hvantk.core' has no attribute '{name}'")
```

- [ ] **Step 3: Verify the shim works**

```bash
python3 -W ignore::DeprecationWarning -c "
from hvantk.core import plugin_loader, tool_loader, anndata_utils, bgzf
print('shim ok:', plugin_loader.__name__, tool_loader.__name__, anndata_utils.__name__, bgzf.__name__)
"
```

Expected output: `shim ok: hvantk.core.plugin.loader hvantk.core.tool.loader hvantk.core.models.anndata_utils hvantk.core.utils.bgzf`

- [ ] **Step 4: Run the plugin-system test sweep**

```bash
/usr/local/bin/pytest hvantk/tests/test_plugin_api.py hvantk/tests/test_plugin_loader.py hvantk/tests/test_plugin_manifest_schema.py hvantk/tests/test_drift_runner.py hvantk/tests/test_plugins_cli.py hvantk/tests/test_drift_cli.py hvantk/tests/test_reprocess_cli.py hvantk/tests/test_plugin_loader_lifecycle.py hvantk/tests/test_plugin_loader_registry_integration.py hvantk/tests/test_plugin_packaging.py hvantk/tests/test_tool_loader.py hvantk/tests/test_tool_manifest_schema.py hvantk/tests/test_tools_cli.py hvantk/tests/test_tool_smoke.py hvantk/tests/test_dependency_directions.py -v -m "" 2>&1 | tail -10
```

Expected: every test that was passing before still passes (plugin_api, manifest_schema, tool_loader, tool_manifest_schema, drift_runner, reprocess_cli, dependency-direction xfails still XFAIL).

A small number of test files reference the OLD import paths in `monkeypatch.setattr("hvantk.core.plugin_loader", ...)` calls — those should be updated, but in this task we rely on the shim. If a test fails with `AttributeError: module 'hvantk.core' has no attribute 'plugin_loader'`, that means `__getattr__` isn't being consulted (this can happen with `monkeypatch.setattr` calls that bypass `__getattr__`). For those test files, update the dotted path explicitly:

```python
monkeypatch.setattr("hvantk.core.plugin.loader.get_registry", lambda: reg)
```

Search and fix any monkeypatch calls touching the old paths:

```bash
grep -rn '"hvantk\.core\.\(plugin_loader\|plugin_api\|tool_loader\|tool_api\|drift_runner\|anndata_utils\|backends\|bgzf\|converters\|hail_context\|metadata\|readers\|writers\|router\)' hvantk/tests/
```

Update each match to the new dotted path under `hvantk.core.{plugin,tool,models,utils}.<basename>`.

- [ ] **Step 5: Commit the entire Phase 1 reorg**

```bash
git add hvantk/core/
git add hvantk/tests/  # any monkeypatch fixes
git commit -m "refactor(core): reorganize into models/utils/plugin/tool sub-packages

Move all 17 files in hvantk/core/ into purpose-driven sub-packages:
- core/plugin/ (api, loader, drift_runner, manifest schema)
- core/tool/ (api, loader, manifest schema)
- core/models/ (anndata_utils, backends, metadata)
- core/utils/ (bgzf, converters, readers, writers, router, hail_context)
- core/ root keeps protocols.py, constants.py, config.py

hvantk/core/__init__.py gains a __getattr__ shim that forwards the prior
flat-namespace imports (from hvantk.core import plugin_loader, etc.) to
the new locations with a DeprecationWarning. Shim removal is scheduled
for Phase 7 of the restructure."
```

---

## Phase 2 — `data/` → `core/{models,utils,streamers}/` + streamer decoupling

### Task 8: Move non-streamer files from `data/` to `core/`

**Files moved:**
- `hvantk/data/dataset.py` → `hvantk/core/models/dataset.py`
- `hvantk/data/file_utils.py` → `hvantk/core/utils/file_utils.py`
- `hvantk/data/gene_mapper.py` → `hvantk/core/utils/gene_mapper.py`

- [ ] **Step 1: Move the 3 files**

```bash
git mv hvantk/data/dataset.py hvantk/core/models/dataset.py
git mv hvantk/data/file_utils.py hvantk/core/utils/file_utils.py
git mv hvantk/data/gene_mapper.py hvantk/core/utils/gene_mapper.py
```

- [ ] **Step 2: Update intra-module imports inside the moved files**

For each moved file, audit absolute imports:

```bash
grep -n 'from hvantk\.\|import hvantk\.' hvantk/core/models/dataset.py hvantk/core/utils/file_utils.py hvantk/core/utils/gene_mapper.py
```

Likely findings:
- `gene_mapper.py` may import from `data/dataset.py` — change `from hvantk.data.dataset` → `from hvantk.core.models.dataset`.
- `file_utils.py` may import from itself or from other utils — usually clean, but verify.

For each absolute import like `from hvantk.data.X`, update to the new path (`hvantk.core.{models,utils}.X` per the move table above).

- [ ] **Step 3: Update CALLERS of these three files**

Run a project-wide search:

```bash
grep -rln 'from hvantk.data.file_utils\|from hvantk.data.gene_mapper\|from hvantk.data.dataset\|hvantk\.data\.file_utils\|hvantk\.data\.gene_mapper\|hvantk\.data\.dataset' \
    hvantk/ examples/ 2>/dev/null | grep -v __pycache__
```

For each match, sed-rewrite:
- `from hvantk.data.file_utils import` → `from hvantk.core.utils.file_utils import`
- `from hvantk.data.gene_mapper import` → `from hvantk.core.utils.gene_mapper import`
- `from hvantk.data.dataset import` → `from hvantk.core.models.dataset import`
- string-literal `"hvantk.data.file_utils"` (in `monkeypatch.setattr`) → `"hvantk.core.utils.file_utils"` (and similarly for the other two)

Expected caller count: about 8 files (7 plugins use `file_utils.download_file`; `gene_mapper` is used by `algorithms/annotation`, `psroc`, `hgnc/SKILL.md` references, etc.).

- [ ] **Step 4: Verify imports resolve**

```bash
python3 -c "
from hvantk.core.utils.file_utils import download_file
from hvantk.core.utils.gene_mapper import GeneMapper
from hvantk.core.models.dataset import Dataset
print('moves ok')
"
```

Expected: `moves ok`. If any `ModuleNotFoundError` other than `hail`, fix the source caller's import.

- [ ] **Step 5: Run the plugin-system tests**

```bash
/usr/local/bin/pytest hvantk/tests/test_plugin_api.py hvantk/tests/test_plugin_loader.py hvantk/tests/test_plugin_manifest_schema.py hvantk/tests/test_tool_loader.py hvantk/tests/test_unified_registry_per_plugin.py -v -m "" 2>&1 | tail -8
```

Expected: all pass.

- [ ] **Step 6: Do NOT commit yet — the streamer moves in Tasks 9–11 are the rest of Phase 2; commit at the end as one atomic unit**

---

### Task 9: Decouple `clinvar_streamer` from `skills/clinvar/`

**Goal:** today's `data/clinvar_streamer.py` calls `from hvantk.skills.clinvar.builder import create_clinvar_tb`. After the move it would be `core → skills`, which violates the dependency rule. Fix: streamer accepts the built table (or a callable that returns one) as a parameter.

**Files modified:**
- `hvantk/data/clinvar_streamer.py` (about to move) — change constructor signature
- Callers of `ClinVarStreamer` — pass the built table

- [ ] **Step 1: Find every caller of `ClinVarStreamer`**

```bash
grep -rn 'ClinVarStreamer\|clinvar_streamer\.ClinVar' hvantk/ examples/ 2>/dev/null | grep -v __pycache__
```

Capture the list. Typical callers (from prior session audit):
- `hvantk/annotation/annotation_streamer.py`
- `hvantk/utils/generate_training_set_streamed.py`
- `hvantk/tests/test_data_streamer.py`
- `hvantk/tests/test_flexible_annotation.py`
- `examples/clinvar/clinvar_streamer_example.py`
- `hvantk/data/__init__.py` (re-export — will be deleted in Task 12)

Read each caller to understand how it currently uses `ClinVarStreamer`. Common pattern: caller passes a path to a raw ClinVar VCF, and the streamer builds the Hail Table internally via `create_clinvar_tb`. After decoupling, the caller will build the table and pass the table itself.

- [ ] **Step 2: Read the current streamer**

```bash
cat hvantk/data/clinvar_streamer.py | head -60
```

Note the current `__init__` signature and where `create_clinvar_tb` is called.

- [ ] **Step 3: Refactor the streamer to accept a pre-built table**

Edit `hvantk/data/clinvar_streamer.py`. Replace the `from hvantk.skills.clinvar.builder import create_clinvar_tb` import with a removal, and change the constructor to take a Hail Table:

```python
# hvantk/data/clinvar_streamer.py  (will move to core/streamers/clinvar.py in Task 10)
"""Stream ClinVar records grouped by gene for downstream consumers.

The streamer operates on a PRE-BUILT Hail Table. The caller is responsible
for producing the table (typically via the clinvar plugin's create_clinvar_tb
builder) and passing it in. This keeps the streamer free of any plugin-layer
dependencies, honoring the core/ → skills/ direction rule documented in
docs/specs/2026-05-19-package-restructure-design.md §5.
"""

from __future__ import annotations

import hail as hl
import logging

logger = logging.getLogger(__name__)


class ClinVarStreamer:
    """Streams ClinVar records grouped by gene from a pre-built Hail Table.

    Parameters
    ----------
    table
        A Hail Table produced by `hvantk.skills.clinvar.builder.create_clinvar_tb`
        (or any equivalent caller-provided source) keyed by (locus, alleles).
    """

    def __init__(self, table: hl.Table):
        self._table = table

    # ... keep the rest of the original streaming logic, but use self._table
    # wherever the old code referenced the internally-built table.
```

Preserve every public method that callers rely on (most of `ClinVarStreamer`'s body stays the same — only the constructor changes).

- [ ] **Step 4: Update each caller to build the table and pass it**

For each caller found in Step 1, replace:

```python
streamer = ClinVarStreamer(vcf_path="/data/clinvar.vcf.bgz", ...)
```

with:

```python
from hvantk.skills.clinvar.builder import create_clinvar_tb

ht = create_clinvar_tb(input_path="/data/clinvar.vcf.bgz", output_path="/tmp/clinvar.ht")
streamer = ClinVarStreamer(table=ht)
```

For tests that use `monkeypatch` to stub the builder, switch to stubbing the table directly:

```python
fake_ht = MagicMock(spec=hl.Table)
streamer = ClinVarStreamer(table=fake_ht)
```

The clinvar example script (`examples/clinvar/clinvar_streamer_example.py`) needs the same change — build first, then stream.

- [ ] **Step 5: Verify nothing else imports from `skills/` inside the streamer**

```bash
grep -n 'hvantk.skills' hvantk/data/clinvar_streamer.py
```

Expected: NO matches. If any remain (e.g., a `from hvantk.skills.clinvar.shared.datasets import ClinVarDataset` for type hints), evaluate whether it's needed at runtime or can be moved behind a `TYPE_CHECKING` guard or removed entirely.

- [ ] **Step 6: Run the relevant tests**

```bash
/usr/local/bin/pytest hvantk/tests/test_data_streamer.py hvantk/tests/test_flexible_annotation.py -v -m "" 2>&1 | tail -10
```

Expected: pass. (These tests will hail-fail in this env, but the import resolution must be clean. If a test fails with a non-hail import error, fix it.)

- [ ] **Step 7: Do NOT commit yet — Task 11 commits all of Phase 2**

---

### Task 10: Move the 6 streamer files into `core/streamers/`

**Files moved:**
- `hvantk/data/data_streamer.py` → `hvantk/core/streamers/base.py`
- `hvantk/data/gene_disease_streamer.py` → `hvantk/core/streamers/gene_disease.py`
- `hvantk/data/clinvar_streamer.py` → `hvantk/core/streamers/clinvar.py`
- `hvantk/data/clingen_streamer.py` → `hvantk/core/streamers/clingen.py`
- `hvantk/data/gencc_streamer.py` → `hvantk/core/streamers/gencc.py`
- `hvantk/data/cosmic_cgc_streamer.py` → `hvantk/core/streamers/cosmic_cgc.py`
- `hvantk/data/alphagenome_streamer.py` → `hvantk/core/streamers/alphagenome.py`

- [ ] **Step 1: Create `core/streamers/` package + move files**

```bash
mkdir -p hvantk/core/streamers
touch hvantk/core/streamers/__init__.py
git mv hvantk/data/data_streamer.py hvantk/core/streamers/base.py
git mv hvantk/data/gene_disease_streamer.py hvantk/core/streamers/gene_disease.py
git mv hvantk/data/clinvar_streamer.py hvantk/core/streamers/clinvar.py
git mv hvantk/data/clingen_streamer.py hvantk/core/streamers/clingen.py
git mv hvantk/data/gencc_streamer.py hvantk/core/streamers/gencc.py
git mv hvantk/data/cosmic_cgc_streamer.py hvantk/core/streamers/cosmic_cgc.py
git mv hvantk/data/alphagenome_streamer.py hvantk/core/streamers/alphagenome.py
```

- [ ] **Step 2: Update intra-streamer imports**

`gene_disease_streamer.py` is the parent class for clingen/gencc/cosmic_cgc. The child classes import it as `from hvantk.data.gene_disease_streamer import GeneDiseaseValidityStreamer`. After move:

```python
# new import inside each child streamer:
from hvantk.core.streamers.gene_disease import GeneDiseaseValidityStreamer
```

Update each of `clingen.py`, `gencc.py`, `cosmic_cgc.py`. Also update any cross-imports between streamers (e.g., a streamer importing the base `DataStreamer` from `data_streamer.py`).

Also update `data_streamer.py` → now `base.py`. If any streamer imports `from hvantk.data.data_streamer import DataStreamer`, change to:

```python
from hvantk.core.streamers.base import DataStreamer
```

- [ ] **Step 3: Update CALLERS of every streamer**

```bash
grep -rln 'from hvantk\.data\.\(clinvar\|clingen\|gencc\|cosmic_cgc\|alphagenome\|gene_disease\|data\)_streamer' hvantk/ examples/ 2>/dev/null | grep -v __pycache__
```

For each match, sed-rewrite:
- `from hvantk.data.clinvar_streamer import` → `from hvantk.core.streamers.clinvar import`
- `from hvantk.data.clingen_streamer import` → `from hvantk.core.streamers.clingen import`
- `from hvantk.data.gencc_streamer import` → `from hvantk.core.streamers.gencc import`
- `from hvantk.data.cosmic_cgc_streamer import` → `from hvantk.core.streamers.cosmic_cgc import`
- `from hvantk.data.alphagenome_streamer import` → `from hvantk.core.streamers.alphagenome import`
- `from hvantk.data.gene_disease_streamer import` → `from hvantk.core.streamers.gene_disease import`
- `from hvantk.data.data_streamer import` → `from hvantk.core.streamers.base import`

Affected files (from prior audit):
- `hvantk/annotation/annotation_streamer.py`
- `hvantk/annotation/annotation_pipeline.py`
- `hvantk/annotation/annotate.py`
- `hvantk/hgc/__init__.py`
- `hvantk/utils/generate_training_set_streamed.py`
- `hvantk/tools/genesets/genesets_cli.py`
- `hvantk/skills/expression_atlas/shared/datasets.py` (uses `load_expression_atlas_datasets` from `unified_registry`; no streamer import — but verify)
- Test files: `test_data_streamer.py`, `test_clingen_streamer.py`, `test_gencc.py`, `test_flexible_annotation.py`
- Example files under `examples/`

For each: update the import paths.

- [ ] **Step 4: Verify imports resolve standalone**

```bash
python3 -c "
import ast
for p in [
    'hvantk/core/streamers/base.py',
    'hvantk/core/streamers/gene_disease.py',
    'hvantk/core/streamers/clinvar.py',
    'hvantk/core/streamers/clingen.py',
    'hvantk/core/streamers/gencc.py',
    'hvantk/core/streamers/cosmic_cgc.py',
    'hvantk/core/streamers/alphagenome.py',
]:
    ast.parse(open(p).read())
    print(f'parse ok: {p}')
"
```

Expected: 7 lines of `parse ok:`.

- [ ] **Step 5: Do NOT commit yet — Task 11 commits all of Phase 2**

---

### Task 11: Remove `hvantk/data/__init__.py` re-exports + delete `data/` + commit Phase 2

**Files modified:**
- Delete: `hvantk/data/__init__.py`
- Delete: `hvantk/data/` (whole directory should be empty after Tasks 8–10)

- [ ] **Step 1: Check `data/` is now empty**

```bash
ls hvantk/data/ 2>/dev/null
```

Expected: only `__init__.py` (and maybe a `__pycache__/`).

- [ ] **Step 2: Read what `data/__init__.py` re-exports**

```bash
cat hvantk/data/__init__.py
```

It typically has lines like:

```python
from hvantk.data.clinvar_streamer import ClinVarStreamer
from hvantk.data.clingen_streamer import ClinGenStreamer
# etc.
```

Any code that imported `from hvantk.data import ClinVarStreamer` (without the submodule) relied on this re-export. Check for such callers:

```bash
grep -rn 'from hvantk\.data import' hvantk/ examples/ 2>/dev/null | grep -v __pycache__
```

For each match, rewrite to the new full path. E.g.:
- `from hvantk.data import ClinVarStreamer` → `from hvantk.core.streamers.clinvar import ClinVarStreamer`

- [ ] **Step 3: Delete `data/`**

```bash
git rm hvantk/data/__init__.py
# remove __pycache__ left over
rm -rf hvantk/data/__pycache__/ 2>/dev/null
rmdir hvantk/data 2>/dev/null
```

Verify:

```bash
ls hvantk/data 2>&1
```

Expected: `ls: hvantk/data: No such file or directory`.

- [ ] **Step 4: Run the relevant test suites**

```bash
/usr/local/bin/pytest \
    hvantk/tests/test_plugin_api.py \
    hvantk/tests/test_plugin_loader.py \
    hvantk/tests/test_plugin_manifest_schema.py \
    hvantk/tests/test_drift_runner.py \
    hvantk/tests/test_plugins_cli.py \
    hvantk/tests/test_drift_cli.py \
    hvantk/tests/test_reprocess_cli.py \
    hvantk/tests/test_plugin_loader_lifecycle.py \
    hvantk/tests/test_plugin_loader_registry_integration.py \
    hvantk/tests/test_plugin_packaging.py \
    hvantk/tests/test_tool_loader.py \
    hvantk/tests/test_tool_manifest_schema.py \
    hvantk/tests/test_tools_cli.py \
    hvantk/tests/test_tool_smoke.py \
    hvantk/tests/test_unified_registry_per_plugin.py \
    hvantk/tests/test_cli_catalog.py \
    hvantk/tests/test_dependency_directions.py \
    -v -m "" 2>&1 | tail -10
```

Expected: every test that passed before still passes. The dependency-direction test's `test_core_does_not_import_upward` MAY now flip from XFAIL to XPASS — if it does, that's the desired outcome (decoupling worked). If it XPASS-es, remove the `@pytest.mark.xfail` decorator on that one function so it permanently enforces the rule.

Specifically: run:

```bash
/usr/local/bin/pytest hvantk/tests/test_dependency_directions.py -v -m "" 2>&1 | tail -10
```

If `test_core_does_not_import_upward` reports XPASS (i.e., the test passed despite the xfail marker), edit `hvantk/tests/test_dependency_directions.py` and DELETE the `@pytest.mark.xfail(...)` decorator on that specific test function. The other tests' xfails remain (algorithms/ migration is Phase 5; that test stays xfailed until then).

- [ ] **Step 5: Commit Phase 2**

```bash
git add hvantk/
git add examples/  # for clinvar_streamer_example.py and similar
git commit -m "refactor(streamers): dissolve hvantk/data/ into core/{models,utils,streamers}/

Phase 2 of the package restructure (docs/specs/2026-05-19-package-restructure-design.md).

Moves:
  data/file_utils.py    -> core/utils/file_utils.py
  data/gene_mapper.py   -> core/utils/gene_mapper.py
  data/dataset.py       -> core/models/dataset.py
  data/data_streamer.py -> core/streamers/base.py
  data/gene_disease_streamer.py -> core/streamers/gene_disease.py
  data/{clinvar,clingen,gencc,cosmic_cgc,alphagenome}_streamer.py
       -> core/streamers/{clinvar,clingen,gencc,cosmic_cgc,alphagenome}.py

Streamer decoupling: ClinVarStreamer (and any other streamer that previously
imported from hvantk.skills.<provider>.builder) now accepts a pre-built Hail
Table via its constructor. Callers build the table first, then construct the
streamer. This eliminates the core->skills import cycle and unblocks the
dependency-direction enforcement test.

If test_dependency_directions.py::test_core_does_not_import_upward now XPASSes,
removes its xfail marker (still kept on the algorithms test pending Phase 5)."
```

---

## Phase 0–2 Verification (run after Task 11)

- [ ] **Step 1: Full plugin/tool/catalog test sweep**

```bash
/usr/local/bin/pytest \
    hvantk/tests/test_plugin_api.py \
    hvantk/tests/test_plugin_loader.py \
    hvantk/tests/test_plugin_manifest_schema.py \
    hvantk/tests/test_plugin_loader_lifecycle.py \
    hvantk/tests/test_plugin_loader_registry_integration.py \
    hvantk/tests/test_plugin_packaging.py \
    hvantk/tests/test_drift_runner.py \
    hvantk/tests/test_drift_cli.py \
    hvantk/tests/test_plugins_cli.py \
    hvantk/tests/test_reprocess_cli.py \
    hvantk/tests/test_tool_loader.py \
    hvantk/tests/test_tool_manifest_schema.py \
    hvantk/tests/test_tools_cli.py \
    hvantk/tests/test_tool_smoke.py \
    hvantk/tests/test_unified_registry_per_plugin.py \
    hvantk/tests/test_cli_catalog.py \
    hvantk/tests/test_top_level_cli.py \
    hvantk/tests/test_dependency_directions.py \
    -v -m "" 2>&1 | tail -10
```

Expected: pass count unchanged from pre-restructure baseline (~75 passing). `test_dependency_directions.py::test_core_does_not_import_upward` either XFAIL or PASS depending on whether the decorator was removed in Task 11; the other two dependency-direction tests are still XFAIL (algorithms not yet migrated).

- [ ] **Step 2: Plugin loader end-to-end smoke**

```bash
python3 -c "
from hvantk.core.plugin import loader as plugin_loader
plugin_loader.reset_registry_for_tests()
reg = plugin_loader.get_registry()
print('providers:', len(reg.list_providers()))
print('datasets:', len(reg.list_datasets()))
print('errors:', len(reg.load_errors()))
"
```

Expected: 13 providers, ~5 datasets visible without hail, ~11 load errors (hail-import only). Same numbers as the pre-Phase-1 baseline — the reorganization is internally consistent.

- [ ] **Step 3: Verify the shim still works for back-compat callers**

```bash
python3 -W ignore::DeprecationWarning -c "
from hvantk.core import plugin_loader
print('shim ok:', plugin_loader.__name__)
"
```

Expected: `shim ok: hvantk.core.plugin.loader`.

- [ ] **Step 4: Confirm hvantk/data is gone**

```bash
ls hvantk/data 2>&1
```

Expected: `ls: hvantk/data: No such file or directory`.

---

## Out of scope (covered by follow-up plans)

- **Phase 3** — `tables/` → `core/plugin/registry.py` + `core/builders/`. Mechanical move; one commit.
- **Phase 4** — `utils/` split into `core/utils/` (generic) and `algorithms/{expression,statistics,training_sets}/` (computation). One commit, ~10 file moves with judgment calls per file documented in the spec.
- **Phase 5** — `algorithms/` package moves. 8 sub-commits (ptm, psroc, qtlcascade, enrichex, hgc, ancestry, annotation, visualization).
- **Phase 6** — Tool import-path updates. 1 commit, ~24 tool manifests' `cli.module:` fields touched.
- **Phase 7** — Remove the `hvantk/core/__init__.py` deprecation shim added in Task 7. Remove the remaining xfail markers from `test_dependency_directions.py`. Delete now-empty package dirs.
- **Phase 8** — Documentation: update `_conventions/SKILL.md`, `docs_site/architecture.md`, CHANGELOG.

Each of those phases follows a similar shape to this plan and can be handed to subagent execution once Phases 0–2 are merged-ready.
