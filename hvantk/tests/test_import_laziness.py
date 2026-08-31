"""Regression guards for the lazy-import work that took ``hvantk --help`` 10.4 s -> 0.1 s.

That speed-up is invisible: nothing fails when it regresses, the CLI just gets
slow again. A single eager ``import pandas`` added to ``hvantk/core/models/__init__.py``
-- or one ``from hvantk.algorithms... import`` at module scope in ``hvantk/hvantk.py``
-- silently restores the old cost, because Python imports a package before any of
its submodules. These tests are the thing that notices.

Every check runs in a subprocess. The pytest session has already imported pandas,
anndata and often Hail by the time any test runs, so an in-process
``sys.modules`` assertion would pass unconditionally and be worse than no test.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

# Modules whose import cost is the whole point: pandas ~0.4 s, anndata ~0.8 s,
# hail ~5.4 s. Kept as a tuple so the failure message can name the culprit.
HEAVY = ("pandas", "anndata", "hail")


def _is_loaded(module: str, loaded: set[str]) -> bool:
    return module in loaded


def _modules_after(code: str) -> set[str]:
    """Run ``code`` in a clean interpreter, return the full module names it loaded."""
    probe = (
        f"{code}\n"
        "import sys\n"
        "print(' '.join(sorted(sys.modules)))\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, f"probe failed:\n{proc.stderr}"
    return set(proc.stdout.split())


def test_importing_provenance_does_not_pull_pandas_or_anndata():
    """The plugin loader imports exactly this, on every ``plugins``/``catalog`` call.

    ``Provenance`` is a plain dataclass, but it lives in a package whose siblings
    import pandas and anndata -- so an eager ``__init__`` makes this dataclass
    cost ~1.2 s.
    """
    loaded = _modules_after("from hvantk.core.models import Provenance")
    offenders = sorted(set(HEAVY) & loaded)
    assert not offenders, (
        f"hvantk.core.models eagerly imported {offenders}. Every name in that "
        f"package must stay behind the PEP 562 __getattr__ in its __init__.py."
    )


def test_importing_the_artifact_union_is_opt_in():
    """``Artifact`` needs all four classes, so it *does* cost pandas + anndata.

    That is fine and expected -- it is why the union lives in its own module and
    is registered in ``_EXPORTS`` rather than built at package import time. This
    test pins the arrangement: asking for ``Artifact`` pays, not asking does not.
    If this ever passes without pandas, the union has been hollowed out.
    """
    loaded = _modules_after("from hvantk.core.models import Artifact")
    assert "pandas" in loaded, (
        "Artifact no longer pulls pandas -- the union is probably no longer a "
        "real union of the concrete classes."
    )


def _subcommand_modules() -> set[str]:
    from hvantk.hvantk import _LAZY_COMMANDS

    return {module for module, _attr, _help in _LAZY_COMMANDS.values()}


def test_building_the_cli_imports_no_subcommand_module():
    """Importing the CLI must resolve nothing.

    Asserting on a proxy like "hail is not loaded" is too weak: most subcommand
    modules are themselves cheap (``hvantk.tools.hgc`` imports in 0.02 s and pulls
    nothing), so an eager import of one of those would sail past a hail-only
    check while still breaking the pattern. The invariant is that *no* entry in
    ``_LAZY_COMMANDS`` is imported, which is what this asserts.
    """
    loaded = _modules_after("from hvantk.hvantk import cli")
    eager = sorted(m for m in _subcommand_modules() if _is_loaded(m, loaded))
    assert not eager, f"CLI import eagerly loaded subcommand modules: {eager}"


def test_listing_commands_imports_no_subcommand_module():
    """``hvantk --help`` renders every short help; it may only read the registry.

    Resolving each command to ask for its help would import all 18 modules --
    ~11 s, including Hail via ``download_cli`` and ``enrichex`` -- which is
    exactly the cost ``_LAZY_COMMANDS`` exists to avoid.
    """
    loaded = _modules_after(
        "from click.testing import CliRunner\n"
        "from hvantk.hvantk import cli\n"
        "r = CliRunner().invoke(cli, ['--help'])\n"
        "assert r.exit_code == 0, r.output\n"
    )
    eager = sorted(m for m in _subcommand_modules() if _is_loaded(m, loaded))
    assert not eager, f"`hvantk --help` imported subcommand modules: {eager}"


def test_artifact_union_covers_every_artifact_class():
    """The union must list every concrete artifact, or ``isinstance`` silently lies.

    Adding a fifth artifact type and forgetting the union makes
    ``isinstance(new_artifact, Artifact)`` return False, which reads as "not an
    artifact" everywhere the union is used to narrow.
    """
    import typing

    from hvantk.core.models import Artifact
    from hvantk.core.models.annotation_table import AnnotationTable
    from hvantk.core.models.expression_matrix import ExpressionMatrix
    from hvantk.core.models.gene_set import GeneSet
    from hvantk.core.models.variant_matrix import VariantMatrix

    expected = {AnnotationTable, ExpressionMatrix, VariantMatrix, GeneSet}
    assert set(typing.get_args(Artifact)) == expected


def test_loader_rejects_the_artifact_union_as_a_manifest_artifact_type(tmp_path):
    """``artifact_type:`` in a plugin.yaml must name a concrete class.

    ``Artifact`` resolves through the same ``getattr(models, name)`` lookup as
    the concrete types, so without a guard a manifest could declare it. That
    would make ``run_builder_for_spec``'s isinstance check accept *any* artifact,
    silently voiding the manifest's declaration -- and because ``types.UnionType``
    has no ``__name__``, the contract-violation message would itself raise
    AttributeError rather than report the mismatch.
    """
    import shutil

    from hvantk.core.plugin.api import PluginLoadError
    from hvantk.core.plugin.loader import PluginRegistry

    src = Path(__file__).parent / "testdata" / "raw" / "plugins" / "fake_plugin"
    dst = tmp_path / "fake_plugin"
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__"))

    manifest = dst / "plugin.yaml"
    manifest.write_text(
        manifest.read_text().replace(
            "    backend: hail\n", "    backend: hail\n    artifact_type: Artifact\n"
        )
    )

    reg = PluginRegistry()
    reg.load_from_directory(dst)
    with pytest.raises(PluginLoadError, match="concrete artifact class"):
        reg.get_dataset("fake:default")


def test_loader_still_accepts_a_concrete_artifact_type(tmp_path):
    """The guard must not reject the legitimate case it sits next to."""
    import shutil

    from hvantk.core.models import AnnotationTable
    from hvantk.core.plugin.loader import PluginRegistry

    src = Path(__file__).parent / "testdata" / "raw" / "plugins" / "fake_plugin"
    dst = tmp_path / "fake_plugin"
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__"))

    manifest = dst / "plugin.yaml"
    manifest.write_text(
        manifest.read_text().replace(
            "    backend: hail\n",
            "    backend: hail\n    artifact_type: AnnotationTable\n",
        )
    )

    reg = PluginRegistry()
    reg.load_from_directory(dst)
    assert reg.get_dataset("fake:default").artifact_type is AnnotationTable
