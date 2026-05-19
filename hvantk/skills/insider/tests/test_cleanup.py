"""Lightweight tests for INSIDER temp-file cleanup."""

from __future__ import annotations

import importlib.util
import sys
import types
import uuid
from pathlib import Path

import pytest

# Resolve the plugin builder source file path without importing it (the import
# pulls in Hail). The test fixture below loads the file with lightweight stubs.
INSIDER_BUILDER_PATH = (
    Path(__file__).resolve().parents[1] / "builder.py"
)


@pytest.fixture
def insider_builder_module(monkeypatch):
    """Load hvantk/skills/insider/builder.py with lightweight dependency stubs.

    The real module imports ``hl``, ``_create_table_base``,
    ``_parse_insider_bed_to_temp_tsv``, and ``_cleanup_temp_file`` from
    ``hvantk.core.builders.table``. We stub those at import time so the
    builder can be loaded without touching Hail or the heavy table-builders
    module, then monkeypatch them per-test to assert cleanup ordering.
    """

    def register(name: str, module: types.ModuleType) -> types.ModuleType:
        monkeypatch.setitem(sys.modules, name, module)
        return module

    hail = register("hail", types.ModuleType("hail"))
    hail.Table = object
    hail.tint32 = object()
    hail.tfloat64 = object()
    hail.tstr = object()
    hail.import_table = lambda *args, **kwargs: None
    hail.locus_interval = lambda *args, **kwargs: None

    # Stub hvantk.core.builders.table so the plugin builder's `from ... import`
    # succeeds without loading the real (Hail-heavy) module. The per-test
    # monkeypatching below rebinds these names on the loaded plugin module.
    register("hvantk", types.ModuleType("hvantk"))
    register("hvantk.core", types.ModuleType("hvantk.core"))
    register("hvantk.core.builders", types.ModuleType("hvantk.core.builders"))
    table_builders = register(
        "hvantk.core.builders.table",
        types.ModuleType("hvantk.core.builders.table"),
    )
    table_builders._create_table_base = lambda **kwargs: None
    table_builders._parse_insider_bed_to_temp_tsv = lambda input_path: "stub.tsv"
    table_builders._cleanup_temp_file = lambda path: None

    module_name = f"insider_builder_under_test_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, INSIDER_BUILDER_PATH)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, module_name, module)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class _FakeImportedTable:
    def __init__(self):
        self.contig = "chr1"
        self.start = 1
        self.end = 2

    def annotate(self, **kwargs):
        return self

    def select(self, *fields):
        return self


def test_create_interactome_tb_cleans_up_temp_tsv_on_success(
    insider_builder_module, monkeypatch
):
    cleanup_calls = []
    parse_calls = []

    fake_hail = types.SimpleNamespace(
        tint32=object(),
        import_table=lambda *args, **kwargs: _FakeImportedTable(),
        locus_interval=lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(insider_builder_module, "hl", fake_hail)
    monkeypatch.setattr(
        insider_builder_module,
        "_parse_insider_bed_to_temp_tsv",
        lambda input_path: parse_calls.append(input_path) or "insider-temp.tsv",
    )
    monkeypatch.setattr(
        insider_builder_module,
        "_cleanup_temp_file",
        lambda path: cleanup_calls.append(path),
    )

    def succeed_after_import(**kwargs):
        kwargs["import_func"]()
        return "ok"

    monkeypatch.setattr(insider_builder_module, "_create_table_base", succeed_after_import)

    result = insider_builder_module.create_interactome_tb(
        input_path="insider.bed",
        output_path="insider.ht",
    )

    assert result == "ok"
    assert parse_calls == ["insider.bed"]
    assert cleanup_calls == ["insider-temp.tsv"]


def test_create_interactome_tb_cleans_up_temp_tsv_on_failure(
    insider_builder_module, monkeypatch
):
    cleanup_calls = []
    parse_calls = []

    fake_hail = types.SimpleNamespace(
        tint32=object(),
        import_table=lambda *args, **kwargs: _FakeImportedTable(),
        locus_interval=lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(insider_builder_module, "hl", fake_hail)
    monkeypatch.setattr(
        insider_builder_module,
        "_parse_insider_bed_to_temp_tsv",
        lambda input_path: parse_calls.append(input_path) or "insider-temp.tsv",
    )
    monkeypatch.setattr(
        insider_builder_module,
        "_cleanup_temp_file",
        lambda path: cleanup_calls.append(path),
    )

    def fail_after_import(**kwargs):
        kwargs["import_func"]()
        raise RuntimeError("checkpoint failed")

    monkeypatch.setattr(insider_builder_module, "_create_table_base", fail_after_import)

    with pytest.raises(RuntimeError, match="checkpoint failed"):
        insider_builder_module.create_interactome_tb(
            input_path="insider.bed",
            output_path="insider.ht",
        )

    assert parse_calls == ["insider.bed"]
    assert cleanup_calls == ["insider-temp.tsv"]
