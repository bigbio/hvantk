"""Assert intra-core directional rules:

    core/io      may import   core/models
    core/models  may NOT import core/io
    core/utils   may NOT import core/models or core/io  (utils stays generic)

Same shape as test_dependency_directions.py: AST walk every import.
"""
from __future__ import annotations

import ast
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _imports_in(subpath: str) -> list[tuple[Path, str]]:
    out: list[tuple[Path, str]] = []
    root = PACKAGE_ROOT / "core" / subpath
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


def _has(forbidden: str, imports: list[tuple[Path, str]]) -> list[tuple[Path, str]]:
    return [(p, m) for p, m in imports if m == forbidden or m.startswith(forbidden + ".")]


def test_core_models_does_not_import_core_io():
    bad = _has("hvantk.core.io", _imports_in("models"))
    assert not bad, (
        "core/models must not import from core/io. Offenders:\n"
        + "\n".join(f"  {p.relative_to(PACKAGE_ROOT)} -> {m}" for p, m in bad)
    )


def test_core_utils_does_not_import_core_models_or_io():
    bad = _has("hvantk.core.models", _imports_in("utils")) + \
          _has("hvantk.core.io", _imports_in("utils"))
    assert not bad, (
        "core/utils must stay generic (no core/models or core/io imports). "
        f"Offenders: {bad}"
    )


def test_provenance_unknown_is_quarantined():
    """Provenance.unknown is for the legacy shim and tests only.

    If a real plugin or algorithm calls it, drift detection silently fails.
    Allowlist: hvantk/core/io/_legacy.py and anything under hvantk/tests/.
    """
    pkg = PACKAGE_ROOT
    allowed = {
        pkg / "core" / "io" / "_legacy.py",
        pkg / "core" / "models" / "provenance.py",  # the definition itself
    }
    offenders: list[Path] = []
    for py in pkg.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        if py in allowed:
            continue
        if "tests" in py.parts:
            continue  # tests can use it freely
        source = py.read_text()
        if "Provenance.unknown" in source:
            offenders.append(py)
    assert not offenders, (
        "Provenance.unknown is restricted to the legacy shim and tests. "
        "Offenders:\n" + "\n".join(f"  {p.relative_to(pkg)}" for p in offenders)
    )


# Reverse-rule guard for issue #120: if hvantk/core/constants.py is
# ever recreated, it must not declare plugin-specific constants.
def test_core_constants_has_no_plugin_specific_blocks():
    """The core/constants.py file was deleted by issue #120's closing PR.

    It may be recreated in the future for legitimate cross-plugin
    constants, but must never reintroduce plugin-prefixed names that
    were moved out to skills/<plugin>/shared/constants.py.

    This test:
      - passes when the file is absent (current state);
      - passes when the file is present and has no plugin-prefixed names;
      - fails when the file is present and any forbidden prefix appears.

    Forbidden prefixes correspond to the plugins whose constants were
    relocated by #120: clingen, gencc, cosmic_cgc, clinvar, hgnc, ucsc,
    expression_atlas, ensembl, alphagenome.
    """
    forbidden_prefixes = (
        "CLINGEN_",
        "GENCC_",
        "COSMIC_",
        "CLINVAR_",
        "HGNC_",
        "UCSC_",
        "EXPRESSION_ATLAS_",
        "ENSEMBL_BIOMART_",
        "ALPHAGENOME_",
    )
    core_constants = PACKAGE_ROOT / "core" / "constants.py"
    if not core_constants.is_file():
        # File deleted by #120's closing PR; nothing to guard against.
        return
    source = core_constants.read_text()
    offenders = [
        prefix for prefix in forbidden_prefixes
        if any(line.lstrip().startswith(prefix) for line in source.splitlines())
    ]
    assert not offenders, (
        "hvantk/core/constants.py must not declare plugin-specific constants. "
        f"Found prefixes: {offenders}. Move each block to "
        "hvantk/skills/<plugin>/shared/constants.py per issue #120."
    )


# Reverse-rule guard for issue #121: core/streamers/ must not import from skills/.
def test_core_streamers_imports_no_skills():
    """core/streamers/*.py must not contain ``from hvantk.skills``."""
    streamers_dir = PACKAGE_ROOT / "core" / "streamers"
    if not streamers_dir.is_dir():
        return
    offenders: list[str] = []
    for py in streamers_dir.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        for i, line in enumerate(py.read_text().splitlines(), 1):
            s = line.lstrip()
            if s.startswith("from hvantk.skills") or s.startswith("import hvantk.skills"):
                offenders.append(f"{py.relative_to(PACKAGE_ROOT)}:{i}: {s}")
    assert not offenders, (
        "core/streamers/ must not import from hvantk.skills. "
        f"Offenders: {offenders}. See issue #121."
    )


# Reverse-rule guard for issue #121: no plugin-specific classes in core/streamers/.
def test_core_streamers_has_no_plugin_specific_classes():
    """No class declared in core/streamers/ may carry a plugin-name prefix."""
    forbidden_prefixes = (
        "ClinGen", "ClinVar", "GenCC", "COSMIC", "CosmicCGC", "HGNC",
        "Ensembl", "MSigDB", "DbNSFP", "UCSC", "AlphaGenome",
        "ExpressionAtlas", "GTEx", "PQTL", "GWAS",
        "GnomAD", "GeVIR", "UniProt", "PeptideAtlas", "CPTAC", "Insider",
    )
    streamers_dir = PACKAGE_ROOT / "core" / "streamers"
    if not streamers_dir.is_dir():
        return
    offenders: list[str] = []
    for py in streamers_dir.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        for i, line in enumerate(py.read_text().splitlines(), 1):
            s = line.lstrip()
            if s.startswith("class "):
                name = s[len("class "):]
                if any(name.startswith(p) for p in forbidden_prefixes):
                    offenders.append(f"{py.relative_to(PACKAGE_ROOT)}:{i}: {s}")
    assert not offenders, (
        "core/streamers/ must not declare plugin-specific classes. "
        f"Concrete plugin streamers belong in skills/<plugin>/streamers.py. "
        f"Offenders: {offenders}. See issue #121."
    )


# Reverse-rule guard for issue #121: core/utils/ must hold only platform-generic files.
def test_core_utils_has_no_provider_specific_files():
    """core/utils/ must not hold provider-specific files."""
    forbidden_prefixes = (
        "clinvar_", "mondo_", "hgnc_", "cosmic_", "gencc_", "clingen_",
        "gtex_", "expression_atlas_", "ucsc_", "alphagenome_",
        "peptideatlas_", "uniprot_", "cptac_", "gene_disease_",
    )
    utils_dir = PACKAGE_ROOT / "core" / "utils"
    if not utils_dir.is_dir():
        return
    offenders: list[str] = []
    for py in utils_dir.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        if any(py.name.startswith(p) for p in forbidden_prefixes):
            offenders.append(str(py.relative_to(PACKAGE_ROOT)))
    assert not offenders, (
        "core/utils/ must not hold provider-specific files. "
        f"Offenders: {offenders}. Move each to skills/<plugin>/ per issue #121."
    )
