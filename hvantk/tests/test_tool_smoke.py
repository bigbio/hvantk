"""Smoke test: every hvantk/tools/**/*.tool.yaml validates and loads cleanly."""

from __future__ import annotations

from hvantk.core.tool import loader as tool_loader


def test_all_tool_manifests_load_without_errors():
    tool_loader.reset_registry_for_tests()
    reg = tool_loader.get_registry()
    errors = reg.load_errors()
    assert errors == [], "Unexpected tool manifest load errors:\n" + "\n".join(
        f"  {p}: {e}" for p, e in errors
    )


def test_all_expected_tool_domains_have_at_least_one_tool():
    """Smoke: every domain in the Literal enum should have at least one manifest."""
    tool_loader.reset_registry_for_tests()
    reg = tool_loader.get_registry()
    domains = {t.domain for t in reg.list_tools()}
    # NOTE: "annotation" intentionally absent — its only tool (annotate_features,
    # a legacy unwired argparse script) was removed; the domain has no CLI tool.
    expected = {
        "plugins",
        "expression",
        "genesets",
        "ptm",
        "ancestry",
        "qtl",
        "enrichex",
        "hgc",
        "infra",
    }
    missing = expected - domains
    assert not missing, f"Domains with no manifested tool: {missing}"


def test_lazy_commands_match_tool_manifests():
    """``_LAZY_COMMANDS`` (hvantk/hvantk.py) must agree with the tool manifests.

    ``hvantk/hvantk.py`` keeps its own ``(module, function, short help)`` copy
    of every top-level command so that ``hvantk --help`` never has to import
    yaml/jsonschema (see the module docstring there). That copy silently
    diverged from the manifests for 14/18 commands (#304): ``hvantk --help``
    and ``hvantk tools list`` described the same command differently. The fix
    keeps the copy in ``hvantk.py`` for the runtime cost reason, but makes
    divergence a test failure instead of something nobody notices -- this test
    pays the yaml/jsonschema import cost so ``hvantk --help`` does not have to.
    """
    from hvantk.hvantk import _LAZY_COMMANDS

    tool_loader.reset_registry_for_tests()
    reg = tool_loader.get_registry()
    by_name = {t.name: t for t in reg.list_tools()}

    mismatches = []
    for name, (module, attr, description) in sorted(_LAZY_COMMANDS.items()):
        spec = by_name.get(name)
        if spec is None:
            mismatches.append(
                f"{name}: registered in _LAZY_COMMANDS but no hvantk/tools/**/*.tool.yaml "
                "manifest declares it"
            )
            continue
        expected = (module, attr, description)
        actual = (spec.cli_module, spec.cli_callable, spec.description)
        if actual != expected:
            mismatches.append(
                f"{name}: _LAZY_COMMANDS {expected!r} != manifest {actual!r} "
                f"({spec.manifest_path})"
            )
    assert not mismatches, (
        "_LAZY_COMMANDS disagrees with the tool manifests:\n  "
        + "\n  ".join(mismatches)
    )
