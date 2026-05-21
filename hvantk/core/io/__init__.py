"""Artifact serialization.

Dispatch:
  *.parquet            -> AnnotationTable (pandas backend)
  *.ht/                -> AnnotationTable (hail backend)
  *.h5ad               -> ExpressionMatrix (anndata backend) [Task 14]
  *.geneset.json       -> GeneSet [Task 15]

Every saved artifact gets a sidecar <path>.provenance.json. Load returns
the artifact with its manifest re-attached as Provenance, or a legacy
unknown provenance if no manifest is found [Task 16].
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from hvantk.core.io._errors import ArtifactTypeError, SchemaIdMismatchError
from hvantk.core.io._formats import (
    load_annotation_table_ht,
    load_annotation_table_parquet,
    load_expression_matrix_h5ad,
    load_expression_matrix_mt,
    load_gene_set_json,
    save_annotation_table_ht,
    save_annotation_table_parquet,
    save_expression_matrix_h5ad,
    save_expression_matrix_mt,
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
        elif path.suffix == ".mt" or path.name.endswith(".mt/"):
            save_expression_matrix_mt(artifact, path)
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


def load(path: str | Path, *, expected_schema_id: str | None = None) -> Any:
    path = Path(path)
    provenance = read_manifest(path)
    if provenance is None:
        from hvantk.core.io._legacy import unknown_provenance_for
        provenance = unknown_provenance_for(path)
    elif expected_schema_id is not None and provenance.schema_id != expected_schema_id:
        raise SchemaIdMismatchError(
            f"schema_id mismatch at {path}: manifest has {provenance.schema_id!r}, "
            f"caller expected {expected_schema_id!r}"
        )
    if path.name.endswith(".geneset.json"):
        return load_gene_set_json(path, provenance)
    if path.suffix == ".parquet":
        return load_annotation_table_parquet(path, provenance)
    if path.suffix == ".ht" or path.name.endswith(".ht/"):
        return load_annotation_table_ht(path, provenance)
    if path.suffix == ".h5ad":
        return load_expression_matrix_h5ad(path, provenance)
    if path.suffix == ".mt" or path.name.endswith(".mt/"):
        return load_expression_matrix_mt(path, provenance)
    raise ArtifactTypeError(f"load: unrecognized extension for {path}")


def load_native(
    path: str | Path,
    *,
    expected_schema_id: str | None = None,
) -> "tuple[Any, Any]":
    """Load an artifact and return its underlying native object plus provenance.

    Useful for algorithms that are legitimately backend-native (e.g. Hail
    genotype workflows, distributed cascade pipelines) and don't need the
    portable Artifact query API. The native object is returned directly —
    no conversion when the artifact's backend matches the file format.

    Returns
    -------
    native_obj
        The native object underlying the artifact:
          - .parquet  -> pandas.DataFrame
          - .ht/      -> hail.Table
          - .h5ad     -> anndata.AnnData
          - .mt/      -> hail.MatrixTable
          - .geneset.json -> list[str] of gene IDs
    provenance
        The Provenance carried by the artifact (or the legacy-shim
        placeholder for files with no sidecar manifest — see
        ``hvantk/core/io/_legacy.py``).

    Raises
    ------
    SchemaIdMismatchError
        If expected_schema_id is given and the manifest's schema_id differs.
    ArtifactTypeError
        If the file extension isn't a known artifact format.

    Examples
    --------
    >>> ht, prov = core_io.load_native("variants.ht")
    >>> # ht is hl.Table; use the full Hail surface
    >>> result = ht.filter(ht.AC > 0).select("locus", "alleles")
    >>> # Save the result with chained provenance
    >>> core_io.save_native(result, "filtered.ht", provenance=Provenance(
    ...     ..., parents=(prov,)
    ... ))
    """
    artifact = load(path, expected_schema_id=expected_schema_id)

    if isinstance(artifact, AnnotationTable):
        if artifact.backend == "hail":
            native = artifact.to_hail()       # zero-cost: returns self._table
        else:
            native = artifact.to_pandas()     # zero-cost for pandas-backed
    elif isinstance(artifact, ExpressionMatrix):
        if artifact.backend == "hail-mt":
            native = artifact.to_hail_mt()    # zero-cost
        else:
            native = artifact.to_anndata()    # zero-cost for anndata-backed
    elif isinstance(artifact, GeneSet):
        native = artifact.to_list()
    else:
        raise ArtifactTypeError(
            f"load_native: unhandled artifact type {type(artifact).__name__}"
        )
    return native, artifact.provenance


def save_native(
    native_obj: Any,
    path: str | Path,
    *,
    provenance: Any,
) -> None:
    """Save a raw native object with provenance, bypassing the artifact wrapper.

    Symmetric to load_native: the caller has a native object (hl.Table,
    hl.MatrixTable, pd.DataFrame, anndata.AnnData, or list[str] for gene sets)
    and wants to persist it with a Provenance sidecar manifest.

    The format is inferred from the path extension:
      - .parquet  -> expects pandas.DataFrame
      - .ht/      -> expects hail.Table (calls ht.write)
      - .h5ad     -> expects anndata.AnnData
      - .mt/      -> expects hail.MatrixTable (calls mt.write)
      - .geneset.json -> expects iterable[str]

    Provenance is required — there's no unknown-placeholder shortcut here;
    the explicit intent of save_native is to carry meaningful provenance
    from algorithm derivations.

    Examples
    --------
    >>> filtered_ht = ht.filter(ht.AC > 0)
    >>> core_io.save_native(filtered_ht, "filtered.ht", provenance=Provenance(
    ...     plugin="qtlcascade", dataset="qtlcascade:filtered",
    ...     plugin_version="0.1.0", source_fingerprint="...",
    ...     schema_id="qtlcascade-filtered-v1",
    ...     build_timestamp=datetime.now(timezone.utc),
    ...     builder_commit=None,
    ...     parents=(input_provenance,),
    ... ))
    """
    path = Path(path)

    # Wrap the native object in the matching artifact and call save().
    # The artifact constructors (from_hail / from_pandas / from_anndata /
    # from_hail_mt) are zero-cost — they just bind self._table / self._matrix
    # to the passed object without copying.
    if path.suffix == ".parquet":
        # pandas.DataFrame expected
        artifact: Any = AnnotationTable.from_pandas(native_obj, provenance=provenance)
    elif path.suffix == ".ht" or path.name.endswith(".ht/"):
        # hail.Table expected
        artifact = AnnotationTable.from_hail(native_obj, provenance=provenance)
    elif path.suffix == ".h5ad":
        # anndata.AnnData expected
        artifact = ExpressionMatrix.from_anndata(native_obj, provenance=provenance)
    elif path.suffix == ".mt" or path.name.endswith(".mt/"):
        # hail.MatrixTable expected
        artifact = ExpressionMatrix.from_hail_mt(native_obj, provenance=provenance)
    elif path.name.endswith(".geneset.json"):
        members = frozenset(native_obj)
        name = path.stem.replace(".geneset", "")
        artifact = GeneSet(name=name, provenance=provenance, _members=members)
    else:
        raise ArtifactTypeError(
            f"save_native: unrecognized extension for {path}"
        )
    save(artifact, path)
