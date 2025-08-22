"""
Registry and adapters for table builders.

Allows running named builders with a unified interface from recipes.

Contract per entry:
- name: str (e.g., "clinvar", "interactome", "gevir", "gnomad-metrics", "ensembl-gene")
- input_path: str
- output_path: str
- params: dict (optional) – builder-specific parameters

This module adapts hvantk.tables.table_builders functions to this contract.
"""
from __future__ import annotations

import logging
from typing import Callable, Dict, Any

logger = logging.getLogger(__name__)

# Lazy import helpers to avoid importing Hail at CLI import time

# Table builders adapters

def _clinvar_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_clinvar_tb

    params = params or {}
    create_clinvar_tb(
        input_path=input_path,
        output_path=output_path,
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
        reference_genome=str(params.get("reference_genome", "GRCh38")),
    )


def _interactome_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_interactome_tb

    params = params or {}
    create_interactome_tb(
        input_path=input_path,
        output_path=output_path,
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
        reference_genome=str(params.get("reference_genome", "GRCh38")),
    )


def _gevir_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_gevir_tb

    params = params or {}
    # fields may be list or comma-separated string
    fields = params.get("fields")
    if isinstance(fields, str):
        fields = [f.strip() for f in fields.split(",") if f.strip()]
    create_gevir_tb(
        input_path=input_path,
        output_path=output_path,
        fields=fields,
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
    )


def _gnomad_metrics_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_gnomad_constraint_gene_metrics_tb

    params = params or {}
    fields = params.get("fields")
    if isinstance(fields, str):
        fields = [f.strip() for f in fields.split(",") if f.strip()]
    create_gnomad_constraint_gene_metrics_tb(
        input_path=input_path,
        output_path=output_path,
        fields=fields,
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
    )


def _ensembl_gene_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_ensembl_gene_tb

    params = params or {}
    fields = params.get("fields")
    if isinstance(fields, str):
        fields = [f.strip() for f in fields.split(",") if f.strip()]
    create_ensembl_gene_tb(
        input_path=input_path,
        output_path=output_path,
        fields=fields,
        canonical=bool(params.get("canonical", True)),
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
    )


def _dbnsfp_adapter(input_path: str, output_path: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.table_builders import create_dbnsfp_tb

    params = params or {}
    prefixes = params.get("group_prefixes")
    if isinstance(prefixes, str):
        prefixes = [p.strip() for p in prefixes.split(",") if p.strip()]

    create_dbnsfp_tb(
        input_path=input_path,
        output_path=output_path,
        reference_genome=str(params.get("reference_genome", "GRCh38")),
        overwrite=bool(params.get("overwrite", False)),
        export_tsv=bool(params.get("export_tsv", False)),
        min_partitions=int(params.get("min_partitions", 200)),
        force_bgz=bool(params.get("force_bgz", True)),
        parse_transcript_scores=bool(params.get("parse_transcript_scores", True)),
        group_prefixes=prefixes,
    )


TABLE_BUILDERS: Dict[str, Callable[[str, str, Dict[str, Any] | None], None]] = {
    "clinvar": _clinvar_adapter,
    "interactome": _interactome_adapter,
    "gevir": _gevir_adapter,
    "gnomad-metrics": _gnomad_metrics_adapter,
    "ensembl-gene": _ensembl_gene_adapter,
    "dbnsfp": _dbnsfp_adapter,
}


def run_table_builder(name: str, input_path: str, output_path: str, params: Dict[str, Any] | None = None) -> None:
    """Run a registered table builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in TABLE_BUILDERS:
        raise KeyError(f"Unknown table builder: {name}")
    logger.info(f"Running builder '{name}' with input={input_path} output={output_path} params={params}")
    TABLE_BUILDERS[name](input_path, output_path, params or {})


# Matrix builders adapters

def _ucsc_mt_adapter(inputs: Dict[str, str], output_mt: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.matrix_builders import build_ucsc_mt

    params = params or {}
    required = ["expression_matrix", "metadata"]
    missing = [k for k in required if k not in inputs or not inputs[k]]
    if missing:
        raise ValueError(f"Missing required inputs for 'ucsc': {missing}")

    build_ucsc_mt(
        expression_matrix_path=inputs["expression_matrix"],
        metadata_path=inputs["metadata"],
        output_mt=output_mt,
        gene_column=str(params.get("gene_column", "gene")),
        metadata_index_col=int(params.get("metadata_index_col", 0)),
        delimiter=str(params.get("delimiter", "\t")),
        min_partitions=int(params.get("min_partitions", 50)),
        force_bgz=bool(params.get("force_bgz", True)),
        split_gene_field=bool(params.get("split_gene_field", True)),
        overwrite=bool(params.get("overwrite", False)),
    )


def _expression_atlas_mt_adapter(inputs: Dict[str, str], output_mt: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.matrix_builders import build_expression_atlas_mt

    params = params or {}
    required = ["expression_matrix", "sdrf"]
    missing = [k for k in required if k not in inputs or not inputs[k]]
    if missing:
        raise ValueError(f"Missing required inputs for 'expression-atlas': {missing}")

    build_expression_atlas_mt(
        expression_matrix_path=inputs["expression_matrix"],
        sdrf_file=inputs["sdrf"],
        output_mt=output_mt,
        gene_column=str(params.get("gene_column", "Gene ID")),
        sample_id_column=str(params.get("sample_id_column", "sample_id")),
        delimiter=str(params.get("delimiter", "\t")),
        min_partitions=int(params.get("min_partitions", 50)),
        force_bgz=bool(params.get("force_bgz", False)),
        overwrite=bool(params.get("overwrite", False)),
    )


def _cptac_mt_adapter(inputs: Dict[str, str], output_mt: str, params: Dict[str, Any] | None = None):
    from hvantk.tables.matrix_builders import build_cptac_mt

    params = params or {}
    required = ["expression", "metadata"]
    missing = [k for k in required if k not in inputs or not inputs[k]]
    if missing:
        raise ValueError(f"Missing required inputs for 'cptac': {missing}")

    # Optional lists: categorical_cols, numeric_cols can be comma-separated strings
    def _split_list(val):
        if val is None:
            return None
        if isinstance(val, list):
            return val
        if isinstance(val, str):
            parts = [x.strip() for x in val.split(",") if x.strip()]
            return parts or None
        return None

    build_cptac_mt(
        expression_path=inputs["expression"],
        metadata_path=inputs["metadata"],
        output_mt=output_mt,
        gene_id_col=str(params.get("gene_id_col", "GeneID")),
        gene_name_col=params.get("gene_name_col", "Gene Name"),
        sample_id_col=str(params.get("sample_id_col", "SampleID")),
        expression_col=str(params.get("expression_col", "Expression")),
        categorical_cols=_split_list(params.get("categorical_cols")),
        numeric_cols=_split_list(params.get("numeric_cols")),
        overwrite=bool(params.get("overwrite", False)),
    )


MATRIX_BUILDERS: Dict[str, Callable[[Dict[str, str], str, Dict[str, Any] | None], None]] = {
    "ucsc": _ucsc_mt_adapter,
    "expression-atlas": _expression_atlas_mt_adapter,
    "cptac": _cptac_mt_adapter,
}


def run_matrix_builder(name: str, inputs: Dict[str, str], output_mt: str, params: Dict[str, Any] | None = None) -> None:
    """Run a registered matrix builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in MATRIX_BUILDERS:
        raise KeyError(f"Unknown matrix builder: {name}")
    logger.info(f"Running matrix builder '{name}' with inputs={inputs} output={output_mt} params={params}")
    MATRIX_BUILDERS[name](inputs, output_mt, params or {})
