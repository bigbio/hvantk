# Refactoring Plan Review and Recommendations

## Executive Summary

After reviewing the proposed refactoring plan and the current codebase, this document provides a **streamlined, minimal-change approach** that focuses on high-impact improvements for library expansion and maintainability while avoiding unnecessary disruption.

## Current State Analysis

### What Already Exists ✅

1. **`core/` module** - Already contains:
   - `config.py` - Configuration management
   - `constants.py` - Shared constants
   - `hail_context.py` - Hail session management

2. **`data/` module** - Already contains:
   - `dataset.py` - Dataset handling
   - `file_utils.py` - File utilities
   - `data_streamer.py` - Data streaming

3. **`tables/` module** - Contains all table builders:
   - `table_builders.py` - All variant/gene annotation builders
   - `matrix_builders.py` - Expression matrix builders
   - `ucsc.py`, `expression_atlas.py`, etc.

4. **`commands/` module** - Well-organized CLI:
   - Individual command files for each feature
   - `hgc_cli.py` - HGC commands
   - `make_table_cli.py`, `make_matrix_cli.py` - Builder CLIs
   - Batch processing support

5. **`hgc/` module** - Mature, feature-complete joint genotyping

### What Needs Improvement 🔧

1. **Module naming** - `tables/` is too generic, doesn't reflect domain
2. **Builder organization** - All builders in single files (200+ lines each)
3. **Protocol definitions** - No formal contracts for builders/streamers
4. **Documentation** - Architecture not clearly documented
5. **Test organization** - All tests in single directory, not module-aligned

## Proposed Changes (Minimal, High-Impact)

### Phase 1: Add Protocol Definitions (Week 1)
**Impact**: High - Establishes contracts for future builders
**Risk**: Low - No existing code needs to change

**Tasks**:
1. Create `hvantk/core/protocols.py` with:
   - `Builder` protocol
   - `Streamer` protocol  
   - `Downloader` protocol
2. Add docstrings and type hints
3. Update documentation

**Files to Create**:
- `hvantk/core/protocols.py` (~100 lines)

### Phase 2: Reorganize Builders by Domain (Week 2)
**Impact**: High - Makes codebase more navigable
**Risk**: Medium - Requires updating imports

**Tasks**:
1. Create domain-specific builder directories:
   ```
   hvantk/
   ├── builders/
   │   ├── __init__.py (re-export all for backward compatibility)
   │   ├── variants/      # ClinVar, dbNSFP, gnomAD variants, CCR
   │   ├── genes/         # Ensembl, GeVIR, gnomAD gene metrics
   │   ├── proteins/      # INSIDER (interactome)
   │   └── expression/    # UCSC, GTEx, bulk RNA-seq
   ```

2. Move existing builders from `tables/table_builders.py` into domain folders
3. Keep `hvantk/tables/__init__.py` as compatibility shim with deprecation notice
4. Update imports in command files

**Backward Compatibility**:
```python
# hvantk/tables/table_builders.py
import warnings
from hvantk.builders.variants.clinvar import create_clinvar_tb
# ... (re-export all with deprecation warning)

warnings.warn(
    "Importing from hvantk.tables.table_builders is deprecated. "
    "Use hvantk.builders.variants/genes/proteins/expression instead.",
    DeprecationWarning,
    stacklevel=2
)
```

### Phase 3: Improve Documentation (Week 3)
**Impact**: High - Makes library more accessible
**Risk**: None - Documentation only

**Tasks**:
1. Create `docs/ARCHITECTURE.md` documenting current structure
2. Update `docs/USAGE.md` with domain-organized examples
3. Add API reference section to README
4. Create `docs/CONTRIBUTING.md` for developers

### Phase 4: Test Organization (Week 4)
**Impact**: Medium - Easier to find and run relevant tests
**Risk**: Low - Just moving files

**Tasks**:
1. Organize tests by module:
   ```
   hvantk/tests/
   ├── unit/
   │   ├── builders/
   │   │   ├── test_variants.py
   │   │   ├── test_genes.py
   │   │   ├── test_proteins.py
   │   │   └── test_expression.py
   │   ├── test_core.py
   │   └── test_commands.py
   ├── integration/
   │   └── test_workflows.py
   └── testdata/ (keep as-is)
   ```

## What We Should NOT Do ❌

### 1. Create New "Kit" Modules (VAK/GAK/PAK/OMEX)
**Reason**: Adds complexity without clear benefit
- Current users don't think in "kits"
- Adds extra layer of indirection
- `builders/variants/` is clearer than `vak/builders/`
- No functional advantage

### 2. Major CLI Restructuring
**Reason**: Breaking changes without user demand
- Current CLI is working well
- `hvantk mktable clinvar` is clear and intuitive
- Changing to `hvantk vak build clinvar` is less clear
- No complaints about current UX

### 3. Pipeline/Recipe Engine (from scratch)
**Reason**: Premature complexity
- Batch processing already exists (`mktable-batch`, `mkmatrix-batch`)
- JSON/YAML recipe support already works
- Building a full orchestration engine is overkill
- Can add features incrementally as needed

### 4. Data Catalog System
**Reason**: Already exists in `resources/`
- `resources/catalog.yaml` already exists
- `resources/registry/` has organized datasets
- No need to reinvent

### 5. Force Migration Schedule
**Reason**: Artificial deadline creates risk
- 7-week timeline is arbitrary
- Better to migrate incrementally as needed
- No urgency driving this

## Recommended Implementation Plan

### Week 1: Protocols and Documentation
- [ ] Create `hvantk/core/protocols.py`
- [ ] Add type hints and docstrings
- [ ] Create `docs/ARCHITECTURE.md`
- [ ] Update README with current architecture

### Week 2: Builder Reorganization
- [ ] Create `hvantk/builders/` directory structure
- [ ] Move variant builders to `builders/variants/`
- [ ] Move gene builders to `builders/genes/`
- [ ] Move protein builders to `builders/proteins/`
- [ ] Move expression builders to `builders/expression/`
- [ ] Add backward compatibility shims in `hvantk/tables/`
- [ ] Update imports in command files

### Week 3: Testing and Polish
- [ ] Reorganize tests to match module structure
- [ ] Run full test suite
- [ ] Update documentation
- [ ] Create migration guide for external users

### Week 4: Review and Release
- [ ] Code review
- [ ] Performance validation (no regressions)
- [ ] Update CHANGELOG
- [ ] Tag release

## Success Metrics

✅ **Must Have**:
- All existing tests pass
- All existing CLI commands work unchanged
- No performance degradation
- Clear architecture documentation

✅ **Nice to Have**:
- Protocol definitions for new builders
- Domain-organized builder code
- Better organized tests
- Migration guide for advanced users

## Conclusion

The original refactoring plan was comprehensive but overly ambitious. This streamlined approach:

1. **Preserves what works**: CLI, HGC, batch processing
2. **Improves what matters**: Code organization, documentation, protocols
3. **Avoids unnecessary change**: No new "kit" modules, no pipeline engine, no CLI restructuring
4. **Delivers real value**: Easier to add new builders, clearer codebase, better docs

**Estimated effort**: 3-4 weeks vs 7 weeks in original plan
**Risk level**: Low vs Medium-High in original plan
**User impact**: Minimal (backward compatible) vs High (breaking changes) in original plan
