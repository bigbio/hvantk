"""Smoke test: Phases E/F/G/H entry-points have @algorithm metadata."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from hvantk.core.models.backends import Backend, get_algorithm_meta


# ---------- Phase E: enrichex ----------

def test_enrichex_burden_analysis_has_metadata():
    from hvantk.algorithms.enrichex.burden import run_burden_analysis

    meta = get_algorithm_meta(run_burden_analysis)
    assert meta.name == "burden_analysis"


def test_enrichex_stratified_burden_has_metadata():
    from hvantk.algorithms.enrichex.burden import run_stratified_burden_analysis

    meta = get_algorithm_meta(run_stratified_burden_analysis)
    assert meta.name == "stratified_burden_analysis"


# ---------- Phase F: expression ----------

def test_expression_tissue_specificity_has_metadata():
    from hvantk.algorithms.expression.tissue_specificity import compute_specificity

    meta = get_algorithm_meta(compute_specificity)
    assert meta.name == "tissue_specificity"


# ---------- Phase G: qtlcascade (already had decorators) ----------

def test_qtlcascade_has_at_least_one_decorated_function():
    """qtlcascade had @algorithm decorators before this PR; verify they still register."""
    import importlib
    import inspect

    mod = importlib.import_module("hvantk.algorithms.qtlcascade.cascade")
    decorated_count = sum(
        1 for _, obj in inspect.getmembers(mod)
        if callable(obj) and hasattr(obj, "_algorithm_meta")
    )
    assert decorated_count > 0


# ---------- Phase H: hgc ----------

def test_hgc_combine_gvcfs_has_metadata():
    from hvantk.algorithms.hgc.combiners import combine_gvcfs

    meta = get_algorithm_meta(combine_gvcfs)
    assert meta.name == "combine_gvcfs"


def test_hgc_combine_vdses_has_metadata():
    from hvantk.algorithms.hgc.combiners import combine_vdses

    meta = get_algorithm_meta(combine_vdses)
    assert meta.name == "combine_vdses"


def test_hgc_convert_vds_to_mt_has_metadata():
    from hvantk.algorithms.hgc.converters import convert_vds_to_mt

    meta = get_algorithm_meta(convert_vds_to_mt)
    assert meta.name == "convert_vds_to_mt"


# ---------- Cross-phase: every domain stays skill-free ----------

def test_algorithm_meta_captures_inputs_outputs():
    from hvantk.core.models import AnnotationTable
    from hvantk.core.models.backends import Backend, algorithm, get_algorithm_meta

    @algorithm(
        name="typed",
        backends=[Backend.PANDAS],
        inputs={"ann": AnnotationTable},
        outputs={"result": AnnotationTable},
    )
    def typed(ann):
        return ann

    meta = get_algorithm_meta(typed)
    assert meta.inputs == {"ann": AnnotationTable}
    assert meta.outputs == {"result": AnnotationTable}
    assert meta.required_backend is None


def test_algorithm_meta_captures_required_backend():
    from hvantk.core.models.backends import Backend, algorithm, get_algorithm_meta

    @algorithm(
        name="hail_only",
        backends=[Backend.HAIL],
        required_backend="hail",
    )
    def hail_only():
        pass

    meta = get_algorithm_meta(hail_only)
    assert meta.required_backend == "hail"


def test_existing_algorithm_meta_still_works():
    """Backward compat: algorithms without inputs/outputs kwargs continue to work."""
    from hvantk.core.models.backends import Backend, algorithm, get_algorithm_meta

    @algorithm(name="legacy_shape", backends=[Backend.HAIL])
    def legacy():
        pass

    meta = get_algorithm_meta(legacy)
    assert meta.name == "legacy_shape"
    assert meta.inputs == {}
    assert meta.outputs == {}


@pytest.mark.parametrize("domain", ["enrichex", "expression", "qtlcascade", "hgc"])
def test_domain_has_no_skill_imports(domain):
    root = Path(__file__).resolve().parents[1] / "algorithms" / domain
    bad = []
    for py in root.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("hvantk.skills"):
                        bad.append((py.name, alias.name))
            elif isinstance(node, ast.ImportFrom) and node.module:
                if node.module.startswith("hvantk.skills"):
                    bad.append((py.name, node.module))
    assert not bad, f"algorithms/{domain}/ has skill imports: {bad}"
