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
    "anndata_schema_to_dict",
    "anndata_sample_rows",
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


def _to_hashable(value: Any) -> Any:
    """Convert a JSON-converted value into a hashable form for dict-key use.

    Variant tables key on `(locus, alleles)` where alleles is `array<str>`,
    which `_to_jsonable` returns as a Python list. Lists are unhashable, so
    they cannot appear inside a tuple used as a dict key without conversion.
    """
    if isinstance(value, list):
        return tuple(_to_hashable(v) for v in value)
    if isinstance(value, dict):
        return tuple(sorted((k, _to_hashable(v)) for k, v in value.items()))
    return value


def _jsonable_to_hail_python(value: Any, dtype: Any) -> Any:
    """Convert JSON-stable snapshot keys back into Hail-compatible values."""
    import hail as hl

    if value is None:
        return None
    if isinstance(dtype, hl.tlocus):
        contig, position = str(value).rsplit(":", 1)
        rg = getattr(dtype.reference_genome, "name", dtype.reference_genome)
        return hl.Locus(contig, int(position), reference_genome=rg)
    if isinstance(dtype, hl.tinterval):
        # Stored as {"start": "<contig>:<pos>", "end": "<contig>:<pos>"} per
        # _to_jsonable. Reconstruct each point via the point-type's branch
        # (recurses into tlocus above for locus<rg>-typed intervals).
        start = _jsonable_to_hail_python(value["start"], dtype.point_type)
        end = _jsonable_to_hail_python(value["end"], dtype.point_type)
        return hl.Interval(start=start, end=end, includes_start=True, includes_end=False)
    if isinstance(dtype, hl.tarray):
        return [_jsonable_to_hail_python(v, dtype.element_type) for v in value]
    if isinstance(dtype, hl.tset):
        return {_jsonable_to_hail_python(v, dtype.element_type) for v in value}
    if isinstance(dtype, hl.ttuple):
        return tuple(
            _jsonable_to_hail_python(v, t)
            for v, t in zip(value, dtype.types)
        )
    if isinstance(dtype, hl.tstruct):
        return {
            name: _jsonable_to_hail_python(value[name], dtype[name])
            for name in dtype
        }
    return value


def collect_sample_rows(table: Any, keys: list[dict]) -> list[dict]:
    """Collect rows whose key fields match one of the provided dicts.

    Returns rows in the same order as `keys`. Comparison is on key fields only.
    Raises KeyError if any requested key is missing from the table — this is
    intentional, since silently dropping missing keys would cause confusing
    snapshot diffs.
    """
    import hail as hl

    if not keys:
        return []

    key_field_names = list(table.key)
    key_dtype = table.key.dtype
    key_rows = [
        {
            name: _jsonable_to_hail_python(k[name], key_dtype[name])
            for name in key_field_names
        }
        for k in keys
    ]
    requested_keys = hl.Table.parallelize(
        key_rows,
        schema=key_dtype,
        key=key_field_names,
    )
    collected = table.semi_join(requested_keys).collect()
    by_key: dict[tuple, dict] = {}
    for row in collected:
        row_dict = dict(row)
        key_tuple = tuple(
            _to_hashable(_to_jsonable(row_dict[k]))
            for k in key_field_names
        )
        by_key[key_tuple] = {
            "key": {k: _to_jsonable(row_dict[k]) for k in key_field_names},
            "row": {
                k: _to_jsonable(v)
                for k, v in row_dict.items()
                if k not in key_field_names
            },
        }

    out: list[dict] = []
    for k in keys:
        key_tuple = tuple(_to_hashable(_to_jsonable(k[name])) for name in key_field_names)
        if key_tuple not in by_key:
            raise KeyError(
                f"collect_sample_rows: requested key {k!r} not found in table; "
                f"available keys: {list(by_key.keys())[:5]} (showing up to 5)"
            )
        out.append(by_key[key_tuple])
    return out


def _to_jsonable(value: Any) -> Any:
    """Convert Hail-collected values into JSON-stable Python primitives.

    Hail Locus is rendered as "<contig>:<position>". Hail Struct is rendered
    as a nested dict by recursively converting each field. Other Hail-specific
    types (Call, Interval) fall through to repr via str(). Sets are sorted by
    their JSON-converted values; this assumes elements are orderable scalars
    (Hail set element types are typed scalars).
    """
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
    try:
        import hail as hl
        if isinstance(value, hl.Locus):
            return f"{value.contig}:{value.position}"
        if isinstance(value, hl.Struct):
            return {k: _to_jsonable(v) for k, v in dict(value).items()}
        if isinstance(value, hl.Interval):
            # Render as {"start": "<contig>:<pos>", "end": "<contig>:<pos>"} —
            # symmetric with locus serialization, unambiguous, round-trippable
            # via _jsonable_to_hail_python's hl.tinterval branch. Assumes the
            # standard half-open BED convention (includes_start=True,
            # includes_end=False).
            return {"start": _to_jsonable(value.start), "end": _to_jsonable(value.end)}
    except ImportError:
        pass
    return str(value)


def anndata_schema_to_dict(adata: Any) -> dict:
    """Serialize an AnnData object's schema to a JSON-stable dict.

    Captures `n_obs`, `n_vars`, sorted `obs_columns` / `var_columns`, the
    dtype and Python-class name of `X`, and the list of layer names. Plain
    `int` casts are used (numpy ints are not JSON-stable).
    """
    import anndata as ad  # noqa: F401  (lazy-import for parity with hl)

    if adata.X is None:
        x_dtype: Any = None
        x_format = "None"
    else:
        x_dtype = str(adata.X.dtype)
        x_format = type(adata.X).__name__

    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_columns": sorted(adata.obs.columns.tolist()),
        "var_columns": sorted(adata.var.columns.tolist()),
        "X_dtype": x_dtype,
        "X_format": x_format,
        # anndata >=0.13 reports a None default layer key (list(adata.layers.keys())
        # == [None]) where <=0.11 reported []; drop it so snapshots are version-stable.
        "layers": sorted(k for k in adata.layers.keys() if k is not None),
    }


def anndata_sample_rows(adata: Any, n: int = 5) -> dict:
    """Collect a JSON-stable head sample from an AnnData.

    Returns a dict with three keys:
      - `obs_head`: first n rows of `adata.obs` as records (index column preserved).
      - `var_head`: first n rows of `adata.var` as records (index column preserved).
      - `X_corner`: first n×n slice of `adata.X` as a nested Python list, or
        `None` if `adata.X is None`. Sparse matrices are densified via `.toarray()`.
    """
    import anndata as ad  # noqa: F401

    obs_head = adata.obs.head(n).reset_index().to_dict(orient="records")
    var_head = adata.var.head(n).reset_index().to_dict(orient="records")

    if adata.X is None:
        x_corner: Any = None
    else:
        sub = adata.X[:n, :n]
        if hasattr(sub, "toarray"):
            sub = sub.toarray()
        x_corner = sub.tolist()

    return {"obs_head": obs_head, "var_head": var_head, "X_corner": x_corner}


def load_snapshot(path: str | Path) -> Any:
    """Load a snapshot JSON file."""
    return json.loads(Path(path).read_text())


def phase_b_snapshot_adapter(builder_fn, dataset_name: str):
    """Adapt a Phase B builder for the snapshot-test calling convention.

    The snapshot helpers below were designed around the legacy Phase A
    ``(input_path, output_path, **kw) -> hl.Table`` shape. Phase B builders
    take ``(parsed_input, ctx, **params) -> Artifact``. This factory bridges
    the two: the returned callable constructs a deterministic ``BuildContext``,
    invokes the Phase B builder, persists via ``artifact.save()``, and returns
    the inner Hail Table so ``regenerate_snapshots`` can introspect it.

    ``dataset_name`` should be the plugin's compound dataset key
    (e.g. ``"clinvar:variants"``).
    """
    from hvantk.core.models.build_context import BuildContext

    plugin = dataset_name.split(":", 1)[0]

    def _adapter(input_path, output_path, **kw):
        # Strip kwargs the Phase A signature accepted but Phase B does not.
        kw.pop("overwrite", None)
        kw.pop("export_tsv", None)
        ctx = BuildContext(
            plugin=plugin,
            dataset=dataset_name,
            plugin_version="test",
            source_fingerprint="sha256:test",
            builder_commit=None,
        )
        artifact = builder_fn(parsed_input=input_path, ctx=ctx, **kw)
        artifact.save(output_path)
        return artifact.to_hail()

    return _adapter


def regenerate_snapshots(
    builder_fn: Callable[..., Any],
    fixture_path: str,
    snapshot_dir: str | Path,
    keys: Iterable[dict] | None = None,
    builder_kwargs: dict | None = None,
    input_path_kwarg: str = "input_path",
) -> None:
    """Run the builder against the fixture and write canonical snapshots.

    Dispatches on the builder return type:
      - If the builder returns an `anndata.AnnData`, writes anndata-shape
        snapshots via `anndata_schema_to_dict` + `anndata_sample_rows`.
      - Otherwise, falls back to reading a Hail Table from `output_path` and
        writing Hail-shape snapshots via `hail_schema_to_dict` +
        `collect_sample_rows`.

    Parameters
    ----------
    builder_fn:
        Builder function under test.
    fixture_path:
        Path to the input fixture file. Passed to the builder under the kwarg
        named by `input_path_kwarg`.
    snapshot_dir:
        Directory where `schema.json` and `sample_rows.json` are written.
    keys:
        For Hail-Table builders, the list of key dicts to extract via
        `collect_sample_rows`. Ignored on the AnnData path. Defaults to an
        empty list when omitted.
    builder_kwargs:
        Extra kwargs to pass to the builder. If `output_path` is present in
        this mapping, it is used as-is (typical for AnnData builders that
        write `.h5ad`). Otherwise a temporary `.ht` path is created and used.
    input_path_kwarg:
        Name of the builder's input-path argument. Defaults to `"input_path"`,
        matching the Hail Table builder convention. The UCSC snapshot wrapper
        uses `"expression_matrix_path"`, for example.

    Writes:
      - <snapshot_dir>/schema.json
      - <snapshot_dir>/sample_rows.json
    """
    from tempfile import TemporaryDirectory

    snapshot_dir = Path(snapshot_dir)
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory() as tmp:
        # Build call kwargs. The Hail path defaults output_path to <tmp>/out.ht;
        # the AnnData path's output_path is supplied by the caller via builder_kwargs.
        call_kwargs: dict[str, Any] = {input_path_kwarg: fixture_path}
        if builder_kwargs and "output_path" in builder_kwargs:
            call_kwargs.update(builder_kwargs)
        else:
            call_kwargs["output_path"] = str(Path(tmp) / "out.ht")
            if builder_kwargs:
                call_kwargs.update(builder_kwargs)

        result = builder_fn(**call_kwargs)

        # Dispatch on return type.
        try:
            import anndata as ad

            if isinstance(result, ad.AnnData):
                schema = anndata_schema_to_dict(result)
                (snapshot_dir / "schema.json").write_text(
                    json.dumps(schema, indent=2, sort_keys=True)
                )
                rows = anndata_sample_rows(result)
                (snapshot_dir / "sample_rows.json").write_text(
                    json.dumps(rows, indent=2, sort_keys=True)
                )
                return
        except ImportError:
            pass

        # Hail Table fallback (existing logic).
        import hail as hl

        ht = hl.read_table(call_kwargs["output_path"])
        schema = hail_schema_to_dict(ht)
        (snapshot_dir / "schema.json").write_text(json.dumps(schema, indent=2, sort_keys=True))

        rows = collect_sample_rows(ht, keys=list(keys or []))
        (snapshot_dir / "sample_rows.json").write_text(json.dumps(rows, indent=2, sort_keys=True))
