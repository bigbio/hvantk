# hvantk Data Handlers Plugin System — Design

- **Date:** 2026-05-15
- **Branch:** `feat/data-handlers-refactoring` (from `dev`)
- **Status:** Approved for implementation planning
- **Owner:** Yasset Perez-Riverol

## 1. Problem

hvantk wraps ten external data providers (clinvar, gtex-eqtl, gwas-catalog, hgnc, insider, msigdb, ucsc-cellbrowser, cptac, peptideatlas, expression-atlas). For each provider, the adapter code is scattered across four to six locations in the repository:

- `hvantk/datasets/<x>_datasets.py` — dataset class, download wrapper
- `hvantk/tables/<x>.py` or a slice of `tables/table_builders.py` / `tables/matrix_builders.py` — Hail or AnnData builder
- `hvantk/commands/<x>_downloader.py` — CLI downloader
- `hvantk/tests/test_<x>_builder.py` — round-trip test
- `hvantk/tests/testdata/raw/<x>/` — fixture data
- `hvantk/tests/snapshots/<x>/` — schema + row snapshots
- `hvantk/skills/<x>/SKILL.md` — agent-facing skill (7 of 10 providers have one)
- `hvantk/tables/registry.py` — registration in `TABLE_BUILDERS` / `MATRIX_BUILDERS`
- `hvantk/resources/catalog.yaml`, `hvantk/resources/registry/<domain>/datasets.json` — metadata

Adding a provider, updating an adapter, or detecting upstream schema drift requires touching six-plus paths. When the upstream source changes its data structure (a new column, a renamed field, a deprecated endpoint), there is no automated signal — the change shows up as a downstream test failure or, worse, as silently wrong data.

## 2. Goals

1. **One package per provider.** Every file that an agent or human needs to touch when fixing, updating, or releasing an adapter lives under a single folder: `hvantk/skills/<provider>/`.
2. **A machine-readable contract per provider.** A `plugin.yaml` manifest declares the provider's name, version, datasets, builders, drift probe, tests, and CLI commands. Agents read the manifest before opening any code.
3. **Active drift detection.** Each plugin ships a lightweight `drift_probe` that fingerprints the live source's header/manifest. A scheduled runner diffs the observed fingerprint against a committed expected fingerprint and flags drift.
4. **Extensibility.** Third parties can publish provider plugins as separate pip packages and register them via Python entry points without forking hvantk.
5. **Backwards-compatible internal API.** Existing callers of `TABLE_BUILDERS` / `MATRIX_BUILDERS` continue to work. CLI commands keep dispatching through `make_table_cli.py` / `make_matrix_cli.py`.

## 3. Non-goals

- Replacing `hvantk/resources/catalog.yaml` or `hvantk/resources/registry/<domain>/datasets.json`. They remain the source of truth for provider metadata; plugins reference them via `source.catalog_ref` / `source.registry_refs`.
- Refactoring non-provider commands (`ptm_cli`, `enrichex_cli`, `make_table_cli`, `make_matrix_cli`, `download_cli`). These stay in `hvantk/commands/`.
- Shipping an out-of-tree example plugin package. The entry-point hook is added; a worked external example is a follow-up doc PR.
- Wiring the scheduled drift workflow into GitHub Actions. `hvantk drift --all --json` is available; the cron job and auto-PR mechanic are a separate PR.

## 4. Decisions

| Question | Decision |
|---|---|
| Co-location scope | Maximum — code + SKILL.md + tests + fixtures + snapshots + manifest all inside `hvantk/skills/<provider>/` |
| Discovery mechanism | Hybrid — filesystem scan of `hvantk/skills/*/plugin.yaml` plus Python entry points (`hvantk.providers` group) |
| Drift detection | Dedicated probe per plugin (header fingerprint + version probe), separate from snapshot tests |
| Granularity | One plugin folder per source; a source may host multiple dataset subfolders (one SKILL.md per dataset) |
| Versioning | Per-plugin semver in `plugin.yaml`; release tag `<provider>-X.Y.Z` |
| Migration scope | All 10 providers in this PR |
| Contract style | YAML manifest + Python `Provider` protocol (Approach A) |
| Test location | Inside the plugin folder: `hvantk/skills/<provider>/<dataset>/tests/` |
| Registry keying | Compound keys `provider:dataset` (e.g., `cptac:phospho`, `clinvar:variants`) |
| Backwards-compat aliases | None — every caller migrates to compound keys in this PR |

## 5. Package layout

Folder shape per provider:

```
hvantk/skills/<provider>/
  plugin.yaml                      # machine-readable manifest (REQUIRED)
  __init__.py                      # package marker; entry-point target (may be empty)
  README.md                        # optional, human-facing
  shared/                          # shared code reused across this source's datasets
    downloader.py
    parser.py
  <dataset>/                       # one folder per dataset shipped by the source
    SKILL.md                       # the 9-section convention (per _conventions/SKILL.md)
    builder.py                     # builder function(s)
    drift_probe.py                 # fingerprint probe
    cli.py                         # optional Click subcommand
    tests/
      test_builder.py
      testdata/raw/<dataset>/...
      snapshots/{schema,sample_rows}.json
      drift_fingerprint.json       # expected fingerprint, compared by drift runner
```

Single-dataset providers (hgnc, msigdb, clinvar, etc.) collapse: `SKILL.md`, `builder.py`, `drift_probe.py`, `tests/` live directly under `hvantk/skills/<provider>/`. The manifest's `datasets:` list has one entry.

## 6. `plugin.yaml` schema (api_version: 1)

```yaml
api_version: 1                     # manifest schema version
name: cptac                        # globally unique provider id
version: 1.2.0                     # plugin semver (drives release tag cptac-1.2.0)
status: stable                     # provisional | stable | deprecated
maintainers:
  - ypriverol@gmail.com
description: CPTAC phospho / acetyl proteomics adapters

source:
  catalog_ref: cptac                       # key in hvantk/resources/catalog.yaml
  registry_refs:                           # entries in resources/registry/<domain>/datasets.json
    - proteomics:CPTAC_phospho_brca

datasets:
  - name: phospho                          # composed with provider as cptac:phospho
    domain: proteomics                     # proteomics | transcriptomics | genomics | epigenomics
    backend: hail                          # hail | anndata | pandas
    builder:
      module: hvantk.skills.cptac.phospho.builder
      function: build_cptac_phospho_tb
    drift_probe:
      module: hvantk.skills.cptac.phospho.drift_probe
      function: fetch_fingerprint
    skill: phospho/SKILL.md                # relative to plugin folder
    tests:
      command: pytest hvantk/skills/cptac/phospho/tests -m hail
      fixture: phospho/tests/testdata/raw/cptac-phospho
      schema_snapshot: phospho/tests/snapshots/schema.json
      row_snapshot: phospho/tests/snapshots/sample_rows.json
      drift_fingerprint: phospho/tests/drift_fingerprint.json

cli:
  - command: cptac-phospho-download        # MUST be prefixed with provider name
    module: hvantk.skills.cptac.phospho.cli
    function: download_cmd

runtime_requirements:
  python:
    - "cptac>=1.5"
  system: []                               # e.g. ["java>=11"] for Hail-backed plugins
```

JSON schema for the manifest lives at `hvantk/core/plugin_manifest.schema.json`; CI validates every `plugin.yaml` against it. The loader rejects manifests that fail validation.

### Compound keys

The registry key is `<plugin.name>:<dataset.name>`. The plugin author never writes the composed string; the loader composes it. For single-dataset providers, the dataset name is descriptive (not `default`):

| Provider | Compound key |
|---|---|
| clinvar | `clinvar:variants` |
| hgnc | `hgnc:lookup` |
| msigdb | `msigdb:genesets` |
| insider | `insider:variants` (currently registered as `interactome`) |
| gtex-eqtl | `gtex-eqtl:eqtls` |
| gwas-catalog | `gwas-catalog:associations` |
| ucsc-cellbrowser | `ucsc-cellbrowser:<collection>` (one per collection) |
| expression-atlas | `expression-atlas:<accession>` (one per accession) |
| cptac | `cptac:phospho` (acetyl/glyco future-ready) |
| peptideatlas | `peptideatlas:phospho` |

## 7. Python `Provider` protocol

New module `hvantk/core/plugin_api.py` is the stable import path for in-tree and out-of-tree plugins.

```python
from typing import Callable, Any, Mapping
from dataclasses import dataclass

@dataclass(frozen=True)
class DatasetSpec:
    name: str                             # compound key, e.g. "cptac:phospho"
    domain: str                           # proteomics | transcriptomics | genomics | epigenomics
    backend: str                          # hail | anndata | pandas
    builder: Callable[..., Any]
    drift_probe: Callable[[], Mapping[str, Any]]
    skill_path: str                       # absolute path to SKILL.md
    test_paths: "TestPaths"

@dataclass(frozen=True)
class TestPaths:
    command: str
    fixture: str
    schema_snapshot: str
    row_snapshot: str
    drift_fingerprint: str

class PluginLoadError(Exception):
    """Raised when a plugin fails the protocol or its requirements."""

@dataclass(frozen=True)
class Provider:
    """Registry record for one provider.

    CONSTRUCTED by the loader from plugin.yaml + resolved callables — plugin
    authors do not subclass or instantiate this. The plugin's __init__.py
    only needs to be importable; the loader materialises this object."""
    name: str                             # matches plugin.yaml name
    version: str                          # matches plugin.yaml version
    datasets: tuple[DatasetSpec, ...]
```

### Builder signature contract

```python
# Hail Table builder (backend: hail)
def build_<dataset>_tb(
    input_path: str,
    output_path: str,
    *,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
    **kwargs,
) -> "hl.Table": ...

# AnnData builder (backend: anndata)
def build_<dataset>_ad(
    expression_matrix_path: str,
    metadata_path: str,
    output_path: str,
    *,
    overwrite: bool = False,
    **kwargs,
) -> "anndata.AnnData": ...
```

Builders MUST be idempotent under `overwrite=True` and MUST checkpoint output to disk. Same contract as today; the protocol makes it load-bearing.

### Drift probe signature contract

```python
def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live source.

    MUST read only headers / manifest / version endpoints — no full download.
    MUST be deterministic on unchanged source (modulo `fetched_at`).
    MUST raise DriftProbeError on transient network failure (runner classifies
    that distinctly from drift).
    """
```

Fingerprint shape:

| Key | Purpose | Comparator behavior |
|---|---|---|
| `probe_version` | Bumped when probe semantics change | Exact match |
| `source_version` | Source-reported version, best effort | Exact match |
| `headers` | `{file: [columns...]}` ordered column names | Ordered list (reordering = drift) |
| `checksums` | `{file: sha256-hex}` of header bytes only | Exact match |
| `fetched_at` | ISO timestamp | **Excluded from comparison** |
| `extras` | Optional `dict` for source-specific signals | Exact match if present |

`probe_version` is independent of plugin `version`. Bumping `probe_version` forces fingerprint regeneration.

## 8. Plugin loader & registry

New module `hvantk/core/plugin_loader.py`:

```python
class PluginRegistry:
    def __init__(self) -> None:
        self._providers: dict[str, Provider] = {}
        self._datasets: dict[str, DatasetSpec] = {}   # keyed by compound key
        self._load_errors: list[tuple[str, Exception]] = []

    def get_provider(self, name: str) -> Provider: ...
    def get_dataset(self, name: str) -> DatasetSpec: ...
    def list_providers(self) -> list[Provider]: ...
    def list_datasets(self, *, domain: str | None = None, backend: str | None = None) -> list[DatasetSpec]: ...
    def load_errors(self) -> list[tuple[str, Exception]]: ...

REGISTRY: PluginRegistry = _build_registry()    # module-level singleton, eager on first access
```

### Discovery flow

1. **Filesystem scan** of `hvantk/skills/*/plugin.yaml`, skipping directories starting with `_` (excludes `_conventions/`, `_hooks/`).
2. **Entry points** at group `hvantk.providers`, resolved via `importlib.metadata.entry_points`.

For each manifest:

1. Validate against `plugin_manifest.schema.json`. Failure → log to `load_errors`, skip plugin.
2. Check `runtime_requirements.python`. Missing dep → `load_errors`, skip plugin (preserves today's optional-`cptac` behavior).
3. Import the plugin package (`importlib.import_module("hvantk.skills.<provider>")`). `ImportError` → `load_errors`, skip.
4. Resolve each `datasets[].builder` and `datasets[].drift_probe` via `importlib.import_module` + `getattr`. Missing function → `load_errors`, skip just that dataset.
5. Construct the `Provider` dataclass from `(plugin.name, plugin.version, tuple(DatasetSpec...))`. Dataclass construction is the validation — type errors here surface as `PluginLoadError`.
6. Insert each `DatasetSpec` into `_datasets` and into the existing registry dicts:
   - `backend == "hail"` → `TABLE_BUILDERS[spec.name] = spec.builder`
   - `backend == "anndata"` → `MATRIX_BUILDERS[spec.name] = spec.builder`
   - `backend == "pandas"` → `PANDAS_BUILDERS[spec.name] = spec.builder` (new dict)

### Collision rule

- Two plugins claiming the same `name` → **hard error**, `PluginLoadError` raised, hvantk import fails.
- An out-of-tree entry-point plugin shadowing an in-tree provider name → **hard error**.
- Two out-of-tree plugins colliding → **hard error**.

Compound keys make dataset-name collisions structurally impossible across providers (the prefix differs); within a provider, YAML list keys are trivially unique.

### Timing

Registry is built **eagerly on first access**, not at `import hvantk`. The existing `hvantk/tables/registry.py` module becomes a thin lazy wrapper: `TABLE_BUILDERS["clinvar:variants"]` triggers a one-time registry build. Cost is sub-second for ~10 plugins because loading is metadata-only (no Hail init).

### Error handling

| Failure | Behavior |
|---|---|
| Manifest schema invalid | Skip plugin, append `load_errors` |
| Optional dep missing | Skip plugin, append `load_errors` |
| Builder/probe function not resolvable | Skip just that dataset, append `load_errors` |
| Module import fails | Skip plugin, append `load_errors` |
| Provider name collision | Hard error, `PluginLoadError` |
| Dataclass construction fails (type mismatch on DatasetSpec/Provider) | Skip plugin, append `load_errors` |

## 9. Drift probe runner

New module `hvantk/core/drift_runner.py`:

```python
def run_drift_check(dataset_name: str, *, timeout: int = 60) -> DriftResult:
    """Resolve dataset, invoke probe, diff against expected fingerprint."""
```

`DriftResult.status` is one of `"clean"`, `"drifted"`, `"probe_failed"`. The runner:

1. Resolves the `DatasetSpec` from `REGISTRY`.
2. Loads `expected = json.load(drift_fingerprint.json)`.
3. Invokes `observed = spec.drift_probe()` under a timeout (`signal.alarm` POSIX, `threading.Timer` fallback). Default 60s.
4. Strips `fetched_at` from both sides; deep-equal comparison.
5. On diff: structured diff captured (added/removed/changed keys).
6. On probe exception: captured as `probe_error`.

### CLI

```bash
hvantk drift cptac:phospho
hvantk drift --all
hvantk drift --domain proteomics
hvantk drift --json cptac:phospho
hvantk drift --regenerate cptac:phospho     # overwrite drift_fingerprint.json
```

Exit codes: `0` clean, `1` drifted, `2` probe_failed, `3` registry/loader error.

### Drift-driven update workflow

1. Scheduled CI job: `hvantk drift --all --json > drift_report.json`.
2. CI opens (or updates) a draft PR per `drifted` plugin, attaching the structured diff and a link to that plugin's `SKILL.md` §8 Update playbook.
3. Agent (Claude, etc.) reads `plugin.yaml` and `SKILL.md`, classifies the change as compatible / breaking / spurious, edits `builder.py` if needed, regenerates snapshots + `drift_fingerprint.json`, bumps `plugin.yaml.version`, pushes to the PR.
4. Human reviewer approves the diff narrative.

Probes are credential-free (HTTP only, no auth secrets) so CI can run them in any context.

## 10. CLI integration

Three CLI surfaces:

### Existing builder commands (registry-driven dispatch)

`hvantk/commands/make_table_cli.py` and `make_matrix_cli.py` are reworked to iterate the registry and auto-generate Click subcommands rather than hard-coded `if name == ...` branches:

```python
@click.group(name="mktable")
def mktable_group(): ...

def _attach_table_commands():
    for spec in REGISTRY.list_datasets(backend="hail"):
        @mktable_group.command(name=spec.name)
        @_raw_input_opt
        @_output_ht_opt
        def _cmd(input_path, output_path, _spec=spec, **kwargs):
            _spec.builder(input_path, output_path, **kwargs)

_attach_table_commands()
```

Plugin-specific kwargs come from the builder's signature via the introspection pattern already used in `tables/registry.py:create_table_adapter`.

`hvantk mktable --help` triggers full registry loading. Sub-second cost; by design.

### New plugin-system commands

```bash
hvantk plugins list
hvantk plugins describe <provider>
hvantk plugins errors
hvantk plugins validate <plugin.yaml-path>

hvantk drift <compound-key>
hvantk drift --all
hvantk drift --domain <domain>
hvantk drift --json <compound-key>
hvantk drift --regenerate <compound-key>
```

These live in `hvantk/commands/plugins_cli.py` and `hvantk/commands/drift_cli.py`.

### Per-plugin custom CLI commands

`plugin.yaml.cli:` block resolves to Click commands attached to the top-level `hvantk` group. Names MUST be prefixed with the provider (`cptac-phospho-download`), enforced at load time — same hard-fail behavior as provider name collisions.

### Files that move

| Current location | After refactor |
|---|---|
| `hvantk/commands/cptac_phospho_downloader.py` | `hvantk/skills/cptac/phospho/cli.py` |
| `hvantk/commands/peptideatlas_phospho_downloader.py` | `hvantk/skills/peptideatlas/phospho/cli.py` |
| `hvantk/commands/ucsc_downloader.py` | `hvantk/skills/ucsc-cellbrowser/<dataset>/cli.py` |

### Files that stay

| Module | Reason |
|---|---|
| `hvantk/commands/make_table_cli.py` | Becomes the registry-driven dispatcher |
| `hvantk/commands/make_matrix_cli.py` | Same, for matrix builders |
| `hvantk/commands/download_cli.py` | Thin umbrella listing all `*-download` plugin commands |
| `hvantk/commands/ptm_cli.py` | Not provider-specific |
| `hvantk/commands/enrichex_cli/*` | Not provider-specific |

### Backwards compatibility

- Bare-name lookups (`TABLE_BUILDERS["clinvar"]`, `hvantk mktable clinvar`) **break**. Every caller migrates to compound keys in this PR.
- Imports from `hvantk.commands.<provider>_downloader` **break**. No re-export shims; the refactor branch is the right moment to break that import path. Conventions doc is updated to point at plugin paths.

## 11. Packaging

`pyproject.toml` changes (Poetry syntax; project uses `poetry.core.masonry.api`):

```toml
[tool.poetry.plugins."hvantk.providers"]
# in-tree plugins register themselves the same way as external ones
clinvar = "hvantk.skills.clinvar"
hgnc = "hvantk.skills.hgnc"
# ... one entry per migrated provider

[tool.pytest.ini_options]
testpaths = ["hvantk/tests", "hvantk/skills"]

[tool.poetry]
exclude = ["hvantk/skills/**/tests/**"]   # keeps plugin.yaml + builder.py IN the wheel
```

Package data and wheel exclusion:

- Snapshots/fixtures (`hvantk/skills/**/tests/testdata/**`, `hvantk/skills/**/tests/snapshots/**`) are tracked in git and present in the sdist.
- `hvantk/skills/**/tests/**` is **excluded from the built wheel** so `pip install hvantk` users don't ship test data. Configured via Poetry's `[tool.poetry] exclude = ["hvantk/skills/**/tests/**"]` (the project uses `poetry.core.masonry.api` as its build backend). The same exclude block keeps `plugin.yaml` IN the wheel — it's a runtime artifact, not a test artifact.

A new test `hvantk/tests/test_plugin_packaging.py` asserts that a fresh wheel install can import every plugin without fixtures present.

## 12. Migration plan

### Phase 0 — Scaffolding (no behavior change)

1. Add `hvantk/core/plugin_api.py`, `plugin_loader.py`, `plugin_manifest.schema.json`, `drift_runner.py`.
2. Add `hvantk/commands/plugins_cli.py`, `drift_cli.py` (return empty results until plugins migrate).
3. Update `hvantk/tables/registry.py` so `TABLE_BUILDERS` / `MATRIX_BUILDERS` are lazy-built on first access.
4. Update `pyproject.toml`: empty `hvantk.providers` entry-points group (under `[tool.poetry.plugins."hvantk.providers"]` — Poetry's syntax for entry points), pytest testpaths, Poetry `exclude` for `hvantk/skills/**/tests/**`.

**Checkpoint:** `pytest hvantk/tests/` green. No provider behavior changed.

### Phase 1 — Per-provider migration (10 sub-commits, in order)

Order chosen so review effort ramps from simplest to most complex:

1. `hgnc` — single-dataset, smallest mover, canonical example.
2. `msigdb` — single-dataset, GMT format.
3. `gtex-eqtl` — single-dataset, parquet.
4. `gwas-catalog` — single-dataset, TSV with multi-row keys.
5. `insider` — single-dataset; recipe examples updated (currently keyed as `interactome`).
6. `clinvar` — single-dataset, VCF, special TSV-flatten branch.
7. `expression-atlas` — multi-dataset (per-accession); first many-datasets-per-provider case.
8. `peptideatlas` — single-dataset (phospho); first SKILL.md authored as part of move.
9. `cptac` — single-dataset (phospho), structured for future siblings; SKILL.md authored.
10. `ucsc-cellbrowser` — multi-dataset (3 collections); largest fixture/snapshot move.

Each sub-commit, atomically:

- `git mv` for builder, downloader, dataset class, test, fixtures, snapshots, SKILL.md (preserves blame).
- Author `plugin.yaml`.
- Author `drift_probe.py` + initial `drift_fingerprint.json` (run probe once, commit output).
- Update intra-plugin imports (`hvantk.tables.cptac` → `hvantk.skills.cptac.shared.parser`, etc.).
- Update `pyproject.toml` entry-points list with the new provider.
- Update recipe example JSONs to use compound keys where applicable.
- Run `pytest hvantk/skills/<provider>/` and `pytest hvantk/tests/` — both green.

### Phase 2 — Cleanup

1. Remove now-empty `hvantk/datasets/`, `hvantk/tables/<provider>.py` files, old paths under `hvantk/tests/snapshots/`, `hvantk/tests/testdata/raw/`.
2. Rewrite `hvantk/skills/_conventions/SKILL.md`:
   - §1 repository map → plugin folder layout
   - §6 registry → compound keys + loader-driven registration
   - §7 CLI pattern → `plugin.yaml.cli:` block
   - §8 test pattern → `hvantk/skills/<provider>/<dataset>/tests/`
   - §9 validation contract → plugin-relative paths
   - Add §12 Drift probe contract (fingerprint schema)
3. Update `docs_site/architecture.md`.
4. Update `CONTRIBUTING.md` with "add a new provider" walkthrough (link `hgnc` as canonical example).

### Phase 3 — Verification

1. Full `pytest hvantk/` green.
2. `hvantk drift --all` exits 0 for every provider against its committed fingerprint.
3. `hvantk plugins list` shows 10 providers, 0 load errors.
4. `hvantk mktable --help` enumerates all builders by compound key.
5. CHANGELOG entry with old-key → new-key migration table.

## 13. Risks

1. **Builder lambdas capturing module-level paths.** Several builders use lambdas that close over `input_path`. After the move, test fixture paths change; pytest discovery catches this immediately. Most common breakage during migration.
2. **`_conventions/SKILL.md` referencing pre-migration paths.** The conventions doc must be updated in lockstep with the last migration commit — a stale post-migration conventions doc is worse than no conventions doc.
3. **Wheel exclusion configuration.** Misconfiguring `[tool.setuptools.packages.find]` exclude rules only surfaces when someone `pip install`s, not during local pytest. Mitigated by `test_plugin_packaging.py` asserting wheel-installed plugin import.
4. **Compound-key migration in recipe JSONs.** Three example recipe files reference bare names (`clinvar`, `interactome`, `ucsc`). All three updated in the relevant per-provider commits. Any internal recipe consumed by tests must be discovered via grep in Phase 0 and tracked in the migration commits.
5. **`hl.import_vcf` / `hl.import_table` path arguments inside moved builders.** When builders move, relative paths in `hl.import_*` calls may break. Tests catch this; reviewer should still spot-check at least one Hail-backed plugin per commit.

## 14. Out of scope for this branch

- Out-of-tree plugin example package (follow-up doc PR).
- Scheduled CI drift workflow (follow-up infra PR — depends on cadence + reviewer routing).
- Folding `hvantk/resources/registry/<domain>/datasets.json` content into per-plugin manifests (future migration; today plugins reference via `source.registry_refs`).
- Provider-level umbrella SKILL.md per source (e.g., a `cptac/SKILL.md` describing CPTAC overall, separate from per-dataset SKILLs). Today's convention is one SKILL.md per dataset; revisit if it stops scaling.
