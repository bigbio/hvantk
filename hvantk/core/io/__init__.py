"""Artifact serialization.

Dispatch:
  *.parquet            -> AnnotationTable (pandas backend)
  *.ht/                -> AnnotationTable (hail backend)
  *.h5ad               -> ExpressionMatrix (anndata backend) [Task 14]
  *.geneset.json       -> GeneSet [Task 15]

Every saved artifact gets a sidecar <path>.provenance.json. Load returns
the artifact with its manifest re-attached as Provenance, or
Provenance.unknown(...) if no manifest is found [Task 16].
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from hvantk.core.io._errors import ArtifactTypeError, SchemaIdMismatchError
from hvantk.core.io._formats import (
    load_annotation_table_ht,
    load_annotation_table_parquet,
    load_expression_matrix_h5ad,
    load_gene_set_json,
    save_annotation_table_ht,
    save_annotation_table_parquet,
    save_expression_matrix_h5ad,
    save_gene_set_json,
)
from hvantk.core.io._manifest import read_manifest, write_manifest
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.gene_set import GeneSet


def save(artifact: Any, path: str | Path) -> None:
    path = Path(path)
    if isinstance(artifact, AnnotationTable):
        if path.suffix == ".parquet":
            save_annotation_table_parquet(artifact, path)
        elif path.suffix == ".ht" or path.name.endswith(".ht/"):
            save_annotation_table_ht(artifact, path)
        else:
            raise ArtifactTypeError(
                f"AnnotationTable save: unrecognized extension for {path}"
            )
        write_manifest(artifact.provenance, path)
        return
    if isinstance(artifact, ExpressionMatrix):
        if path.suffix == ".h5ad":
            save_expression_matrix_h5ad(artifact, path)
        else:
            raise ArtifactTypeError(
                f"ExpressionMatrix save: unrecognized extension for {path}"
            )
        write_manifest(artifact.provenance, path)
        return
    if isinstance(artifact, GeneSet):
        if path.name.endswith(".geneset.json"):
            save_gene_set_json(artifact, path)
        else:
            raise ArtifactTypeError(
                f"GeneSet save: must end in .geneset.json, got {path}"
            )
        write_manifest(artifact.provenance, path)
        return
    raise ArtifactTypeError(
        f"save: no handler for artifact type {type(artifact).__name__}"
    )


def load(path: str | Path) -> Any:
    path = Path(path)
    provenance = read_manifest(path)
    if provenance is None:
        raise ArtifactTypeError(
            f"load: no provenance manifest at {path}; legacy shim lands in Task 16"
        )
    if path.name.endswith(".geneset.json"):
        return load_gene_set_json(path, provenance)
    if path.suffix == ".parquet":
        return load_annotation_table_parquet(path, provenance)
    if path.suffix == ".ht" or path.name.endswith(".ht/"):
        return load_annotation_table_ht(path, provenance)
    if path.suffix == ".h5ad":
        return load_expression_matrix_h5ad(path, provenance)
    raise ArtifactTypeError(f"load: unrecognized extension for {path}")
