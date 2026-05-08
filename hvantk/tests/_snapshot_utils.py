"""Snapshot utilities for skill round-trip tests.

Provides JSON-stable schema serialization, key-matched row extraction,
and a regeneration helper used when builder output legitimately changes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Iterable

__all__ = [
    "hail_schema_to_dict",
    "collect_sample_rows",
    "load_snapshot",
    "regenerate_snapshots",
]


def _hail_type_to_str(t: Any) -> str | dict:
    """Render a Hail type as a JSON-stable string or nested dict."""
    import hail as hl

    if isinstance(t, hl.tstruct):
        return {name: _hail_type_to_str(t[name]) for name in t}
    return str(t)


def hail_schema_to_dict(table_or_mt: Any) -> dict:
    """Serialize a Hail Table or MatrixTable schema to a JSON-stable dict."""
    import hail as hl

    if isinstance(table_or_mt, hl.MatrixTable):
        return {
            "row_key": list(table_or_mt.row_key),
            "col_key": list(table_or_mt.col_key),
            "row": _hail_type_to_str(table_or_mt.row.dtype),
            "col": _hail_type_to_str(table_or_mt.col.dtype),
            "entry": _hail_type_to_str(table_or_mt.entry.dtype),
        }
    return {
        "key": list(table_or_mt.key),
        "row": _hail_type_to_str(table_or_mt.row.dtype),
    }


def collect_sample_rows(table: Any, keys: list[dict]) -> list[dict]:
    """Collect rows whose key fields match one of the provided dicts.

    Returns rows in the same order as `keys`. Comparison is on key fields only.
    """
    import hail as hl

    key_field_names = list(table.key)
    collected = table.collect()
    by_key: dict[tuple, dict] = {}
    for row in collected:
        row_dict = dict(row)
        key_tuple = tuple(_to_jsonable(row_dict[k]) for k in key_field_names)
        by_key[key_tuple] = {
            "key": {k: _to_jsonable(row_dict[k]) for k in key_field_names},
            "row": {k: _to_jsonable(v) for k, v in row_dict.items() if k not in key_field_names},
        }

    out: list[dict] = []
    for k in keys:
        key_tuple = tuple(_to_jsonable(k[name]) for name in key_field_names)
        if key_tuple in by_key:
            out.append(by_key[key_tuple])
    return out


def _to_jsonable(value: Any) -> Any:
    """Convert Hail-collected values into JSON-stable Python primitives."""
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, set):
        return sorted(_to_jsonable(v) for v in value)
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    # Hail Locus and similar
    if hasattr(value, "contig") and hasattr(value, "position"):
        return f"{value.contig}:{value.position}"
    return str(value)


def load_snapshot(path: str | Path) -> Any:
    """Load a snapshot JSON file."""
    return json.loads(Path(path).read_text())


def regenerate_snapshots(
    builder_fn: Callable[..., Any],
    fixture_path: str,
    snapshot_dir: str | Path,
    keys: Iterable[dict],
    builder_kwargs: dict | None = None,
) -> None:
    """Run the builder against the fixture and write canonical snapshots.

    Writes:
      - <snapshot_dir>/schema.json
      - <snapshot_dir>/sample_rows.json
    """
    import hail as hl
    from tempfile import TemporaryDirectory

    snapshot_dir = Path(snapshot_dir)
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory() as tmp:
        output_path = str(Path(tmp) / "out.ht")
        builder_fn(input_path=fixture_path, output_path=output_path, **(builder_kwargs or {}))
        ht = hl.read_table(output_path)

        schema = hail_schema_to_dict(ht)
        (snapshot_dir / "schema.json").write_text(json.dumps(schema, indent=2, sort_keys=True))

        rows = collect_sample_rows(ht, keys=list(keys))
        (snapshot_dir / "sample_rows.json").write_text(json.dumps(rows, indent=2, sort_keys=True))
