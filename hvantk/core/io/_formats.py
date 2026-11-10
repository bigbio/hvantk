"""Per-format handlers. Dispatch table lives in __init__.py."""

from __future__ import annotations

import gzip
import json as _json
import math
import numbers
from pathlib import Path
from typing import Any

import pandas as pd

from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.gene_set import GeneSet
from hvantk.core.models.provenance import Provenance
from hvantk.core.models.variant_matrix import VariantMatrix


def save_annotation_table_parquet(ann: AnnotationTable, path: Path) -> None:
    df = ann.to_pandas()
    df.to_parquet(path, index=False)


def load_annotation_table_parquet(
    path: Path, provenance: Provenance
) -> AnnotationTable:
    df = pd.read_parquet(path)
    return AnnotationTable.from_pandas(df, provenance=provenance)


_TSV_FORBIDDEN = ("\t", "\n", "\r")
_TSV_ARRAY_TYPES = (list, tuple, set, frozenset)


def _tsv_is_missing(v: Any) -> bool:
    if v is None or v is pd.NA or v is pd.NaT:
        return True
    if hasattr(v, "item") and not isinstance(v, (str, bytes)):
        try:
            v = v.item()  # np.float32('nan') is not a python float until here
        except (ValueError, TypeError):
            return False
    return isinstance(v, float) and math.isnan(v)


def _tsv_cell(v: Any, dtype: str) -> str:
    """Render one value as a TSV cell (see ``save_annotation_table_tsv``).

    Raises ValueError for an array element that is itself an array or contains
    ``|`` (the joined cell would be ambiguous); the caller adds column and key.
    """
    if hasattr(v, "tolist") and getattr(v, "ndim", 0) > 0:
        v = v.tolist()  # numpy array -> list
    if isinstance(v, _TSV_ARRAY_TYPES):
        items = v
        if isinstance(items, (set, frozenset)):
            items = sorted(items, key=lambda x: (_tsv_is_missing(x), str(x)))
        parts = []
        for x in items:
            if hasattr(x, "tolist") and getattr(x, "ndim", 0) > 0:
                x = x.tolist()
            if isinstance(x, _TSV_ARRAY_TYPES):
                raise ValueError("nested array element")
            part = _tsv_cell(x, "str")
            if "|" in part:
                raise ValueError(f"array element {part!r} contains '|'")
            parts.append(part)
        return "|".join(parts)
    if _tsv_is_missing(v):
        return ""
    if hasattr(v, "item") and not isinstance(v, (str, bytes)):
        v = v.item()  # numpy scalar -> python scalar
    if isinstance(v, bool):  # before int: bool is an int subclass
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # pandas holds an int column with a missing value as float64
        if dtype == "int" and v.is_integer():
            return str(int(v))
        return repr(v)
    return str(v)


def _tsv_sort_rows(
    columns: dict[str, list[Any]], key_fields: list[str], schema: dict[str, str]
) -> list[int]:
    """Row order for the export: stable, missing last, numbers numerically.

    Sorts by ``key_fields`` (those present in ``columns``), then by every
    remaining field left to right so ties are broken deterministically; with no
    key field, by every field. Works on per-column python values, so large ints
    are never upcast.
    """
    fields = list(columns)
    order = [k for k in key_fields if k in columns]
    order += [f for f in fields if f not in order]
    plans = []
    for c in order:
        vals = [v for v in columns[c] if not _tsv_is_missing(v)]
        numeric = all(
            isinstance(v, numbers.Real) and not isinstance(v, bool) for v in vals
        )
        plans.append((c, numeric, schema.get(c, "str")))

    def key(i: int) -> tuple:
        out = []
        for c, numeric, dtype in plans:
            v = columns[c][i]
            if _tsv_is_missing(v):
                out.append((1, 0, ""))
            elif numeric:
                out.append((0, v.item() if hasattr(v, "item") else v, ""))
            else:
                try:
                    text = _tsv_cell(v, dtype)
                except ValueError:  # refused later, with column and key
                    text = repr(v)
                out.append((0, 0, text))
        return tuple(out)

    return sorted(range(len(columns[fields[0]]) if fields else 0), key=key)


def save_annotation_table_tsv(ann: AnnotationTable, path: Path) -> None:
    """Export an AnnotationTable as deterministic gzipped TSV.

    Header is the table's fields in order. Rows are sorted in pandas for both
    backends (``to_pandas()`` row order is not a documented contract): by the
    hail key fields present in the frame (``to_pandas`` flattens struct keys, so
    a key such as ``locus`` may be absent), then by every remaining field left
    to right so ties are deterministic; a pandas table has no key, so by every
    field. The sort is stable, numbers sort numerically, missing values last.
    Array cells join their elements with ``|`` (sets sorted first); a missing
    value is an empty cell; booleans are ``true``/``false``; integers print as
    plain decimals; floats as ``repr``. An empty array, a missing value and
    ``[""]`` all write the same empty cell. A cell containing a tab, newline or
    carriage return, an array element containing ``|``, or a nested array is
    refused (ValueError naming the column and row key) before anything is
    written. Two saves of one table are byte-identical (UTF-8, LF line ends,
    gzip ``mtime=0`` and no filename in the header).
    """
    df = ann.to_pandas()
    fields = [str(c) for c in df.columns]
    df.columns = fields
    columns = {c: df[c].tolist() for c in fields}  # python values, no upcasting
    hail_key = [str(k) for k in ann._table.key] if ann.backend == "hail" else []
    key_fields = [k for k in hail_key if k in columns]
    id_fields = key_fields or fields[:1]
    lines = ["\t".join(fields)]
    for i in _tsv_sort_rows(columns, key_fields, ann.schema):
        cells = []
        for c in fields:
            try:
                cell = _tsv_cell(columns[c][i], ann.schema.get(c, "str"))
                if any(ch in cell for ch in _TSV_FORBIDDEN):
                    raise ValueError(f"cell {cell!r} holds a tab, newline or CR")
            except ValueError as e:
                key = ", ".join(
                    f"{k}={_tsv_cell(columns[k][i], 'str')}"
                    if not isinstance(columns[k][i], _TSV_ARRAY_TYPES)
                    else f"{k}=<array>"
                    for k in id_fields
                )
                raise ValueError(
                    f"TSV export: column {c!r} (row key {key}): {e}"
                ) from None
            cells.append(cell)
        lines.append("\t".join(cells))
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mtime=0, mode="wb", fileobj=raw) as gz:
            gz.write(payload)


def save_annotation_table_ht(ann: AnnotationTable, path: Path) -> None:
    ht = ann.to_hail()
    ht.write(str(path), overwrite=True)


def load_annotation_table_ht(path: Path, provenance: Provenance) -> AnnotationTable:
    import hail as hl

    ht = hl.read_table(str(path))
    return AnnotationTable.from_hail(ht, provenance=provenance)


def save_expression_matrix_h5ad(em: ExpressionMatrix, path: Path) -> None:
    em.to_anndata().write_h5ad(str(path))


def load_expression_matrix_h5ad(path: Path, provenance: Provenance) -> ExpressionMatrix:
    import anndata as ad

    adata = ad.read_h5ad(str(path))
    return ExpressionMatrix.from_anndata(adata, provenance=provenance)


def save_variant_matrix_mt(vm: "VariantMatrix", path: Path) -> None:
    """Write a VariantMatrix as a Hail MatrixTable directory (``.mt/``)."""
    mt = vm.to_hail_mt()
    mt.write(str(path), overwrite=True)


def load_variant_matrix_mt(path: Path, provenance: Provenance) -> "VariantMatrix":
    """Load a Hail MatrixTable directory and wrap it as a VariantMatrix."""
    import hail as hl

    mt = hl.read_matrix_table(str(path))
    return VariantMatrix.from_hail_mt(mt, provenance=provenance)


def save_gene_set_json(gs: GeneSet, path: Path) -> None:
    path.write_text(
        _json.dumps(
            {
                "name": gs.name,
                "members": sorted(gs.to_list()),
            },
            indent=2,
        )
    )


def load_gene_set_json(path: Path, provenance: Provenance) -> GeneSet:
    payload = _json.loads(path.read_text())
    return GeneSet(
        name=payload["name"],
        provenance=provenance,
        _members=frozenset(payload["members"]),
    )
