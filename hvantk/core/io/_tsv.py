"""Deterministic gzipped-TSV export of an AnnotationTable (save only).

Format
------
- UTF-8, LF line ends, gzip with ``mtime=0`` and no file name in the header, so two
  saves of the same table are byte-identical.
- Header: the field names in order; struct fields are flattened to ``a.b``. In a
  hail table an array (or set) of structs becomes one array field per struct field,
  also named ``a.b``, whose elements line up.
- Rows of a hail table are sorted by Hail: by the key fields, then by every other
  field left to right (loci in genomic order, missing values last). A pandas table
  has no key, so its rows are sorted in Python by every field left to right: numbers
  numerically, missing values last, ties broken by the written text.
- Cells: a missing value (or NaN) is empty; booleans are ``true``/``false``; integers
  are plain decimals; floats use Python's ``repr`` (float32 at float32 precision, so
  ``0.1`` rather than ``0.10000000149011612``); loci and intervals use Hail's text
  form (``chr1:100``).
- Array and set cells join their elements with ``|`` (a pandas set is sorted first).
  Inside an element ``%`` is written ``%25`` and ``|`` is written ``%7C``, so a reader
  splits on ``|`` and percent-decodes each element. An empty array, a missing value
  and ``[""]`` all write an empty cell.
- There is no quoting: a ``"`` is written as is, so read the file with quote handling
  off (polars ``quote_char=None``, pandas ``quoting=csv.QUOTE_NONE``).
- Refused with a ValueError naming the column and the row key, before anything is
  written at the output path: a tab, newline or carriage return in a field name or a
  value, a dict or struct value, and an array element that is itself an array,
  struct or dict.

Neither backend materialises the table a second time: a hail table is sorted and
JSON-encoded by Hail into a temporary UTF-8 text file that is read back row by row from
local disk (copied there first when Hail's temporary directory is remote), and the
output is written to ``<path>.partial`` and renamed into place only once complete.
"""

from __future__ import annotations

import gzip
import json
import math
import numbers
import os
import tempfile
import urllib.parse
from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np
import pandas as pd

from hvantk.core.models.annotation_table import AnnotationTable

_FORBIDDEN = ("\t", "\n", "\r")
_LINES_PER_WRITE = 4096


def _is_missing(v: Any) -> bool:
    """True for None, pd.NA, NaT and NaN (python or numpy)."""
    if v is None or v is pd.NA or v is pd.NaT:
        return True
    return isinstance(v, (float, np.floating)) and math.isnan(v)


def _check_text(text: str) -> str:
    """Return ``text``, refusing a tab, newline or carriage return in it."""
    if any(ch in text for ch in _FORBIDDEN):
        raise ValueError(f"value {text!r} holds a tab, newline or carriage return")
    return text


def _scalar(v: Any) -> str:
    """Render one non-container value."""
    if _is_missing(v):
        return ""
    if isinstance(v, np.float32):
        return str(v)  # shortest float32 text, same style as repr
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, bool):  # before int: bool is an int subclass
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, dict):
        raise ValueError("a struct or dict value")
    return _check_text(str(v))


def _element(x: Any) -> str:
    """Render one array element, percent-encoding ``%`` and ``|``."""
    if isinstance(x, (list, tuple, set, frozenset, dict, np.ndarray)):
        raise ValueError("an array element that is an array, struct or dict")
    return _scalar(x).replace("%", "%25").replace("|", "%7C")


def _cell(v: Any) -> str:
    """Render one value as a TSV cell (format in the module docstring)."""
    if isinstance(v, np.ndarray):
        v = v[()] if v.ndim == 0 else list(v)  # list() keeps numpy scalar types
    if isinstance(v, (set, frozenset)):
        return "|".join(t for _, t in sorted((_is_missing(x), _element(x)) for x in v))
    if isinstance(v, (list, tuple)):
        return "|".join(_element(x) for x in v)
    return _scalar(v)


def _safe_cell(v: Any) -> str:
    """Render a key value for an error message without raising."""
    try:
        return _cell(v)
    except ValueError:
        return repr(v)


def _refusal(column: str, key: str, err: ValueError) -> ValueError:
    return ValueError(f"TSV export: column {column!r} (row key {key}): {err}")


def _check_header(fields: list[str]) -> None:
    for f in fields:
        try:
            _check_text(f)
        except ValueError as e:
            raise ValueError(f"TSV export: field name {f!r}: {e}") from None


def _write(path: Path, header: list[str], rows: Iterable[list[str]]) -> None:
    """Stream ``rows`` to a gzipped TSV; nothing appears at ``path`` on failure."""
    partial = path.with_name(path.name + ".partial")
    try:
        with (
            open(partial, "wb") as raw,
            gzip.GzipFile(filename="", mtime=0, mode="wb", fileobj=raw) as gz,
        ):
            buf = ["\t".join(header)]
            for cells in rows:
                buf.append("\t".join(cells))
                if len(buf) >= _LINES_PER_WRITE:
                    gz.write(("\n".join(buf) + "\n").encode("utf-8"))
                    buf = []
            if buf:
                gz.write(("\n".join(buf) + "\n").encode("utf-8"))
        os.replace(partial, path)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise


# --- pandas backend: sort in Python, render every cell once ---


def _column_values(series: pd.Series) -> list[Any]:
    """Python values of one column; float32 values stay float32."""
    values = series.tolist()  # python ints keep every digit
    if str(series.dtype).lower() == "float32":
        values = [None if _is_missing(x) else np.float32(x) for x in values]
    return values


def _sort_key(v: Any, text: str, numeric: bool) -> tuple:
    if _is_missing(v):
        return (1, 0, "")
    if numeric:
        return (0, float(v) if isinstance(v, np.floating) else v, text)
    return (0, 0, text)


def _pandas_rows(df: pd.DataFrame) -> tuple[list[str], Iterator[list[str]]]:
    fields = [str(c) for c in df.columns]
    _check_header(fields)
    columns = [_column_values(df.iloc[:, j]) for j in range(len(fields))]
    rendered, keys = [], []
    for name, values in zip(fields, columns, strict=True):
        cells = []
        for i, v in enumerate(values):
            try:
                cells.append(_cell(v))
            except ValueError as e:
                raise _refusal(
                    name, f"{fields[0]}={_safe_cell(columns[0][i])}", e
                ) from None
        present = [v for v in values if not _is_missing(v)]
        numeric = all(
            isinstance(v, numbers.Real) and not isinstance(v, (bool, np.bool_))
            for v in present
        )
        rendered.append(cells)
        keys.append(
            [_sort_key(v, t, numeric) for v, t in zip(values, cells, strict=True)]
        )
    order = sorted(range(len(df)), key=lambda i: tuple(k[i] for k in keys))
    return fields, ([col[i] for col in rendered] for i in order)


# --- hail backend: Hail sorts and JSON-encodes, Python streams the rows back ---


def _float_kind(dtype: Any) -> tuple[bool, bool] | None:
    """(is a container, is float32) for float fields and float arrays/sets."""
    import hail as hl

    floats = (hl.tfloat32, hl.tfloat64)
    if dtype in floats:
        return (False, dtype == hl.tfloat32)
    if isinstance(dtype, (hl.tarray, hl.tset)) and dtype.element_type in floats:
        return (True, dtype.element_type == hl.tfloat32)
    return None


def _decode_float(v: Any, f32: bool) -> Any:
    if v is None:
        return None
    v = float(v)  # NaN and +/-Infinity arrive as the JSON strings "NaN", "Infinity"
    return np.float32(v) if f32 else v


def _struct_field(name: str):
    return lambda s: s[name]


def _split_struct_arrays(ht: Any) -> Any:
    """Turn each array (or set) of structs into one array field per struct field.

    ``lsdb: array<struct{name, url}>`` becomes ``lsdb.name`` and ``lsdb.url``, two
    arrays whose elements line up; deeper nesting is left for the renderer to refuse.
    """
    import hail as hl

    exprs = {}
    for f in ht.row:
        e, t = ht[f], ht[f].dtype
        if isinstance(t, (hl.tarray, hl.tset)) and isinstance(
            t.element_type, hl.tstruct
        ):
            arr = hl.array(e) if isinstance(t, hl.tset) else e
            for sub in t.element_type.fields:
                exprs[f"{f}.{sub}"] = arr.map(_struct_field(sub))
        else:
            exprs[f] = e
    return ht.select(**exprs)


def _hail_export(ann: AnnotationTable) -> tuple[str, list[str], list, list[int]]:
    """Sort and JSON-encode a hail table into a temporary text file.

    Returns the temp path, the flattened field names, each field's float kind and
    the positions of the fields that identify a row in an error message.
    """
    import hail as hl

    ht = ann.to_hail()
    key = list(ht.key)
    order = key + [f for f in ht.row if f not in key]
    ht = ht.order_by(*[ht[f] for f in order]).flatten()  # unkeyed, order kept
    ht = _split_struct_arrays(ht)
    fields = list(ht.row)
    _check_header(fields)
    key_pos = [
        j
        for j, f in enumerate(fields)
        if any(f == k or f.startswith(k + ".") for k in key)
    ]
    kinds = [_float_kind(ht[f].dtype) for f in fields]
    exprs = {}
    for f in fields:
        e = ht[f]
        if isinstance(e.dtype, hl.tdict):
            raise ValueError(f"TSV export: column {f!r}: a dict field is refused")
        if isinstance(e.dtype, (hl.tlocus, hl.tinterval)):
            e = hl.str(e)
        exprs[f] = hl.json(e)  # never missing: a missing value is "null"
    tmp = hl.utils.new_temp_file(prefix="tsv_export", extension="tsv")
    ht.select(**exprs).export(tmp, header=False, parallel=None)
    return tmp, fields, kinds, key_pos or [0]


def _local_export(tmp: str, workdir: str) -> str:
    """A local path to Hail's export: the file itself when it is local, else a copy.

    The export is read on the local filesystem, never line by line through
    ``hl.hadoop_open``: each read there crosses into the JVM, which ran at about
    0.2 MB/s (an hour for a 210,000 x 312 table) and could stall for good.
    """
    parsed = urllib.parse.urlparse(tmp)
    if parsed.scheme in ("", "file"):
        return parsed.path if parsed.scheme else tmp
    import hail as hl

    local = os.path.join(workdir, "export.tsv")
    hl.hadoop_copy(tmp, "file://" + local)
    return local


def _hail_rows(
    local: str, fields: list[str], kinds: list, key_pos: list[int]
) -> Iterator[list[str]]:
    """Stream the JSON-encoded export back as rendered rows, in Hail's order."""
    # Hail writes the export as UTF-8 ("RORα" must not come back as "RORÎ±").
    with open(local, encoding="utf-8") as fh:
        for line in fh:
            values = [json.loads(t) for t in line.rstrip("\n").split("\t")]
            for j, kind in enumerate(kinds):
                if kind is None:
                    continue
                container, f32 = kind
                if container and values[j] is not None:
                    values[j] = [_decode_float(x, f32) for x in values[j]]
                elif not container:
                    values[j] = _decode_float(values[j], f32)
            cells = []
            for name, v in zip(fields, values, strict=True):
                try:
                    cells.append(_cell(v))
                except ValueError as e:
                    key = ", ".join(
                        f"{fields[k]}={_safe_cell(values[k])}" for k in key_pos
                    )
                    raise _refusal(name, key, e) from None
            yield cells


def save_annotation_table_tsv(ann: AnnotationTable, path: Path) -> None:
    """Export an AnnotationTable as deterministic gzipped TSV (see module docstring)."""
    if ann.backend != "hail":
        fields, rows = _pandas_rows(ann._table)
        _write(Path(path), fields, rows)
        return
    from hvantk.core.utils.hail_helpers import cleanup_temp_file

    tmp, fields, kinds, key_pos = _hail_export(ann)
    try:
        with tempfile.TemporaryDirectory(prefix="hvantk-tsv-") as workdir:
            local = _local_export(tmp, workdir)
            _write(Path(path), fields, _hail_rows(local, fields, kinds, key_pos))
    finally:
        cleanup_temp_file(tmp)
