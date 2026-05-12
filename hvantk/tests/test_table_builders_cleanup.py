"""Lightweight tests for INSIDER temp-file cleanup."""

from __future__ import annotations

import importlib.util
import sys
import types
import uuid
from pathlib import Path

import pytest

TABLE_BUILDERS_PATH = Path(__file__).resolve().parents[1] / "tables" / "table_builders.py"


@pytest.fixture
def table_builders_module(monkeypatch):
    """Load table_builders.py with lightweight dependency stubs."""

    def register(name: str, module: types.ModuleType) -> types.ModuleType:
        monkeypatch.setitem(sys.modules, name, module)
        return module

    hail = register("hail", types.ModuleType("hail"))
    hail.Table = object
    hail.tint32 = object()
    hail.tfloat64 = object()
    hail.tstr = object()
    hail.utils = types.SimpleNamespace(
        new_temp_file=lambda *args, **kwargs: f"{uuid.uuid4().hex}.tsv"
    )
    hail.import_table = lambda *args, **kwargs: None
    hail.locus_interval = lambda *args, **kwargs: None

    register("hvantk", types.ModuleType("hvantk"))
    register("hvantk.core", types.ModuleType("hvantk.core"))
    register("hvantk.data", types.ModuleType("hvantk.data"))
    register("hvantk.utils", types.ModuleType("hvantk.utils"))

    table_utils = register(
        "hvantk.utils.table_utils", types.ModuleType("hvantk.utils.table_utils")
    )
    table_utils.get_row_fields = lambda ht: []
    table_utils.build_rename_map = lambda *args, **kwargs: {}
    table_utils.str_to_bool = lambda value: value

    metadata = register(
        "hvantk.core.metadata", types.ModuleType("hvantk.core.metadata")
    )
    metadata.build_table_metadata = lambda source_name, input_path, ht: {
        "source_name": source_name,
        "input_path": input_path,
    }

    constants = register(
        "hvantk.core.constants", types.ModuleType("hvantk.core.constants")
    )
    for name in (
        "ENSEMBL_BIOMART_FIELDS",
        "CLINGEN_GENE_DISEASE_FIELDS",
        "CLINGEN_CLASSIFICATION_LEVELS",
        "GENCC_SUBMISSION_FIELDS",
        "GENCC_CLASSIFICATION_LEVELS",
        "HGNC_GENE_FIELDS",
        "HGNC_PIPE_SEPARATED_FIELDS",
        "COSMIC_CGC_FIELDS",
        "COSMIC_CGC_CLASSIFICATION_LEVELS",
        "COSMIC_MUTATION_CONTEXTS",
    ):
        setattr(constants, name, {})

    file_utils = register(
        "hvantk.data.file_utils", types.ModuleType("hvantk.data.file_utils")
    )
    file_utils.resolve_compression = lambda path: path

    genome = register("hvantk.utils.genome", types.ModuleType("hvantk.utils.genome"))
    genome.contig_recoding = lambda *args, **kwargs: {}

    module_name = f"table_builders_under_test_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, TABLE_BUILDERS_PATH)
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
    table_builders_module, monkeypatch
):
    cleanup_calls = []
    parse_calls = []

    fake_hail = types.SimpleNamespace(
        tint32=object(),
        import_table=lambda *args, **kwargs: _FakeImportedTable(),
        locus_interval=lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(table_builders_module, "hl", fake_hail)
    monkeypatch.setattr(
        table_builders_module,
        "_parse_insider_bed_to_temp_tsv",
        lambda input_path: parse_calls.append(input_path) or "insider-temp.tsv",
    )
    monkeypatch.setattr(
        table_builders_module,
        "_cleanup_temp_file",
        lambda path: cleanup_calls.append(path),
    )

    def succeed_after_import(**kwargs):
        kwargs["import_func"]()
        return "ok"

    monkeypatch.setattr(table_builders_module, "_create_table_base", succeed_after_import)

    result = table_builders_module.create_interactome_tb(
        input_path="insider.bed",
        output_path="insider.ht",
    )

    assert result == "ok"
    assert parse_calls == ["insider.bed"]
    assert cleanup_calls == ["insider-temp.tsv"]


def test_create_interactome_tb_cleans_up_temp_tsv_on_failure(
    table_builders_module, monkeypatch
):
    cleanup_calls = []
    parse_calls = []

    fake_hail = types.SimpleNamespace(
        tint32=object(),
        import_table=lambda *args, **kwargs: _FakeImportedTable(),
        locus_interval=lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(table_builders_module, "hl", fake_hail)
    monkeypatch.setattr(
        table_builders_module,
        "_parse_insider_bed_to_temp_tsv",
        lambda input_path: parse_calls.append(input_path) or "insider-temp.tsv",
    )
    monkeypatch.setattr(
        table_builders_module,
        "_cleanup_temp_file",
        lambda path: cleanup_calls.append(path),
    )

    def fail_after_import(**kwargs):
        kwargs["import_func"]()
        raise RuntimeError("checkpoint failed")

    monkeypatch.setattr(table_builders_module, "_create_table_base", fail_after_import)

    with pytest.raises(RuntimeError, match="checkpoint failed"):
        table_builders_module.create_interactome_tb(
            input_path="insider.bed",
            output_path="insider.ht",
        )

    assert parse_calls == ["insider.bed"]
    assert cleanup_calls == ["insider-temp.tsv"]
