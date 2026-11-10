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


def _tsv_is_missing(v: Any) -> bool:
    if v is None or v is pd.NA or v is pd.NaT:
        return True
    return isinstance(v, float) and math.isnan(v)


def _tsv_cell(v: Any, dtype: str) -> str:
    """Render one value as a TSV cell (see ``save_annotation_table_tsv``)."""
    if hasattr(v, "tolist") and getattr(v, "ndim", 0) > 0:
        v = v.tolist()  # numpy array -> list
    if isinstance(v, (list, tuple, set, frozenset)):
        items = v
        if isinstance(items, (set, frozenset)):
            items = sorted(items, key=lambda x: (_tsv_is_missing(x), str(x)))
        return "|".join(_tsv_cell(x, "str") for x in items)
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


def _tsv_sort_key(df: pd.DataFrame, cols: list[str], schema: dict[str, str]) -> Any:
    """Per-row sort key over ``cols``: missing last, numbers numerically."""
    plans = []
    for c in cols:
        vals = [v for v in df[c] if not _tsv_is_missing(v)]
        numeric = all(
            isinstance(v, numbers.Real) and not isinstance(v, bool) for v in vals
        )
        plans.append((c, numeric, schema.get(c, "str")))

    def key(row: Any) -> tuple:
        out = []
        for c, numeric, dtype in plans:
            v = row[c]
            if _tsv_is_missing(v):
                out.append((1, 0.0, ""))
            elif numeric:
                out.append((0, float(v), ""))
            else:
                out.append((0, 0.0, _tsv_cell(v, dtype)))
        return tuple(out)

    return key


def save_annotation_table_tsv(ann: AnnotationTable, path: Path) -> None:
    """Export an AnnotationTable as deterministic gzipped TSV.

    Header is the table's fields in order. Rows are sorted by the key fields
    (hail: ``ht.key``, which a keyed Table already iterates in; pandas has no key,
    so every field left to right, stably, missing values last). Array cells join
    their elements with ``|`` (sets sorted first); a missing value is an empty
    cell; booleans are ``true``/``false``; integers print as plain decimals;
    floats as ``repr``. A cell containing a tab, newline or carriage return is
    refused (ValueError naming the column and row key) before anything is
    written. Two saves of one table are byte-identical (UTF-8, LF line ends,
    gzip ``mtime=0`` and no filename in the header).
    """
    df = ann.to_pandas()
    fields = [str(c) for c in df.columns]
    df.columns = fields
    hail_key = [str(k) for k in ann._table.key] if ann.backend == "hail" else []
    rows = [r for _, r in df.iterrows()]
    if not hail_key:  # keyed hail tables already iterate in key order
        rows.sort(key=_tsv_sort_key(df, fields, ann.schema))  # stable
    id_fields = hail_key or fields[:1]
    lines = ["\t".join(fields)]
    for r in rows:
        cells = []
        for c in fields:
            cell = _tsv_cell(r[c], ann.schema.get(c, "str"))
            if any(ch in cell for ch in _TSV_FORBIDDEN):
                key = ", ".join(f"{k}={_tsv_cell(r[k], 'str')}" for k in id_fields)
                raise ValueError(
                    f"TSV export: column {c!r} (row key {key}) holds a tab, "
                    f"newline or carriage return: {cell!r}"
                )
            cells.append(cell)
        lines.append("\t".join(cells))
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
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
