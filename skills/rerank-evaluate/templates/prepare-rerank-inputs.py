"""Put a cohort table, a feature axis, and a label set into one identifier space for
``hvantk rerank``, which does no identifier mapping of its own
(`skills/rerank-evaluate/references/inputs-and-identifiers.md`).

Three subcommands:

    prepare-rerank-inputs.py labels --geneset LABELS.json [--group NAME ...] \\
        [--restrict-to TABLE --restrict-column COL] --out labels.txt
        Pure Python, no Hail: union every group of a GeneSetCollection JSON (or only
        the named --group values), optionally restrict to a universe, and write one
        gene per line, sorted.

    prepare-rerank-inputs.py cohort --table TABLE --key-column COL \\
        --gene-index GENE_INDEX.ht --hgnc-tsv HGNC.txt --out OUT.tsv \\
        [--manifest-in COHORT.yaml --manifest-out OUT.yaml]
        Needs Hail: map TABLE's --key-column (Ensembl gene_id) to an approved HGNC
        symbol via GENE_INDEX's own hgnc_id field, add it as a new 'symbol' column,
        and write every original column alongside it. Rows with no approved symbol,
        or whose symbol is shared by more than one gene id, are dropped and counted.
        With --manifest-in/--manifest-out, also write a symbol-keyed cohort manifest
        pointing at OUT, carrying the original's prior and cohort_axes over unchanged.

    prepare-rerank-inputs.py features --ht TABLE.ht --fields F1,F2,... \\
        --hgnc-tsv HGNC.txt [--gene-index GENE_INDEX.ht] \\
        [--restrict-to TABLE --restrict-column COL] --out OUT.tsv
        Needs Hail: read TABLE.ht (any single-key gene-keyed hvantk table), select
        --fields, map its key to an approved HGNC symbol the same way as `cohort`
        (via --gene-index if TABLE.ht itself carries no hgnc_id field; defaults to
        TABLE.ht), and write a TSV with a literal 'gene' column (the symbol) plus the
        requested fields -- the shape hvantk.algorithms.rerank.catalog.builders.table_axis
        requires. --restrict-to/--restrict-column narrow the output to a universe
        (e.g. the cohort's own mapped gene list).

The approved-symbol map (`build_symbol_map`) is shared by `cohort` and `features` so a
gene id that drops out of one drops out of the other identically: it reads the HGNC
complete-set TSV's raw `hgnc_id`/`symbol`/`status` columns (Approved rows only -- the
same columns `hvantk.algorithms.rerank.blocks.load_gene_groups` reads for --blocks,
NOT the built `hgnc:lookup` table, whose own builder renames `symbol` to
`gene_symbol`), drops gene ids with no approved symbol, then drops every gene id
whose symbol is shared by more than one gene id (both sides, not just one).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


# --- shared: gene_id -> approved HGNC symbol, 1:1 only ------------------------------


def build_symbol_map(gene_index_path: str, hgnc_tsv_path: str):
    """``(mapping, report)``: ``mapping`` is ``{gene_id: symbol}``, unambiguous only.

    Needs Hail (reads *gene_index_path* as a Hail Table); the caller is expected to
    have called ``init_hail`` already.
    """
    import hail as hl
    import pandas as pd

    ht = hl.read_table(gene_index_path)
    key_fields = list(ht.key)
    if len(key_fields) != 1:
        raise ValueError(
            f"{gene_index_path}: expected a single-field key (a gene-keyed table); "
            f"got {key_fields}"
        )
    key_field = key_fields[0]
    if "hgnc_id" not in ht.row:
        raise ValueError(
            f"{gene_index_path}: table has no 'hgnc_id' field; point --gene-index at "
            "a table that carries one (e.g. the gene index built by "
            "`hvantk annotate spine`)"
        )
    gi = ht.select("hgnc_id").to_pandas().rename(columns={key_field: "gene_id"})
    n_total = len(gi)

    hgnc = pd.read_csv(hgnc_tsv_path, sep="\t", dtype=str, low_memory=False)
    missing_cols = [c for c in ("hgnc_id", "symbol") if c not in hgnc.columns]
    if missing_cols:
        raise ValueError(f"{hgnc_tsv_path}: missing column(s) {missing_cols}")
    if "status" in hgnc.columns:
        hgnc = hgnc[hgnc["status"] == "Approved"]
    hgnc = hgnc.drop_duplicates("hgnc_id", keep="first")
    symbol_by_hgnc = dict(zip(hgnc["hgnc_id"], hgnc["symbol"]))

    gi["symbol"] = gi["hgnc_id"].map(symbol_by_hgnc)
    no_symbol = gi["symbol"].isna()
    n_no_symbol = int(no_symbol.sum())
    mapped = gi.loc[~no_symbol, ["gene_id", "symbol"]]

    dup_symbols = set(mapped["symbol"][mapped["symbol"].duplicated(keep=False)])
    ambiguous_mask = mapped["symbol"].isin(dup_symbols)
    n_ambiguous_dropped = int(ambiguous_mask.sum())
    clean = mapped.loc[~ambiguous_mask]

    mapping = dict(zip(clean["gene_id"], clean["symbol"]))
    report = {
        "n_total": n_total,
        "n_no_symbol": n_no_symbol,
        "n_distinct_ambiguous_symbols": len(dup_symbols),
        "n_ambiguous_dropped": n_ambiguous_dropped,
        "n_clean": len(mapping),
    }
    return mapping, report


def _report_line(report: dict) -> str:
    return (
        f"gene index symbol map: {report['n_clean']}/{report['n_total']} gene(s) "
        f"mapped to an unambiguous approved symbol ({report['n_no_symbol']} with no "
        f"approved symbol, {report['n_ambiguous_dropped']} dropped across "
        f"{report['n_distinct_ambiguous_symbols']} ambiguous symbol(s))"
    )


def _sep_for(path) -> str:
    return "\t" if str(path).endswith((".tsv", ".bgz", ".gz")) else ","


def _read_column_csv(path: str, column: str) -> set:
    """One column of a plain (uncompressed) TSV/CSV, as a set of stripped strings."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=_sep_for(path))
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: no column {column!r} (have {reader.fieldnames})")
        return {
            row[column].strip()
            for row in reader
            if row.get(column) and row[column].strip()
        }


# --- labels ---------------------------------------------------------------------------


def cmd_labels(args) -> int:
    from hvantk.core.utils.gene_sets import GeneSetCollection

    collection = GeneSetCollection.load(args.geneset)
    names = args.group if args.group else collection.names()
    unknown = [n for n in names if collection.get(n) is None]
    if unknown:
        raise ValueError(
            f"{args.geneset}: no group(s) {unknown}; available: {collection.names()}"
        )
    genes: set = set()
    for name in names:
        genes |= collection.get(name).genes
    n_before = len(genes)

    n_after = n_before
    if args.restrict_to:
        if not args.restrict_column:
            raise ValueError("--restrict-to requires --restrict-column")
        restrict = _read_column_csv(args.restrict_to, args.restrict_column)
        genes &= restrict
        n_after = len(genes)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for g in sorted(genes):
            fh.write(g + "\n")
    print(
        f"labels: {len(names)} group(s) ({', '.join(sorted(names))}), {n_before} "
        f"gene(s) before restriction, {n_after} after -> {out_path}"
    )
    return 0


# --- cohort -----------------------------------------------------------------------


def cmd_cohort(args) -> int:
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import pandas as pd

    mapping, report = build_symbol_map(args.gene_index, args.hgnc_tsv)
    print(_report_line(report))

    df = pd.read_csv(args.table, sep=_sep_for(args.table), dtype={args.key_column: str})
    if args.key_column not in df.columns:
        raise ValueError(f"{args.table}: no column {args.key_column!r}")
    df["symbol"] = df[args.key_column].map(mapping)
    n_before = len(df)
    kept = df[df["symbol"].notna()].reset_index(drop=True)
    n_dropped = n_before - len(kept)
    print(
        f"cohort: kept {len(kept)}/{n_before} row(s) ({n_dropped} dropped: no "
        "approved symbol, or a symbol shared by more than one gene id)"
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept.to_csv(out_path, sep=_sep_for(out_path), index=False)
    print(f"wrote {out_path}")

    if args.manifest_out:
        if not args.manifest_in:
            raise ValueError("--manifest-out requires --manifest-in")
        import yaml

        from hvantk.algorithms.cohort.spec import load_cohort

        manifest = load_cohort(args.manifest_in)
        doc = {
            "name": manifest.name,
            "key": "symbol",
            "key_column": "symbol",
            # Kept verbatim, exactly as the caller spelled --out: a cohort manifest's
            # `table` is resolved by whoever reads it later, relative to ITS working
            # directory, not to the manifest file's own location
            # (hvantk/algorithms/cohort/frame.py), so the caller should pass the same
            # path form (absolute, ideally) it will use to read the config later.
            "table": str(args.out),
            "prior": {
                "column": manifest.prior.column,
                "direction": manifest.prior.direction,
            },
        }
        if manifest.cohort_axes:
            doc["cohort_axes"] = [
                {"axis": a.axis, "columns": list(a.columns)}
                for a in manifest.cohort_axes
            ]
        manifest_out = Path(args.manifest_out)
        manifest_out.parent.mkdir(parents=True, exist_ok=True)
        manifest_out.write_text(yaml.safe_dump(doc, sort_keys=False))
        print(f"wrote {manifest_out} (key=symbol, table={doc['table']})")
    return 0


# --- features -----------------------------------------------------------------------


def cmd_features(args) -> int:
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    gene_index_path = args.gene_index or args.ht
    mapping, report = build_symbol_map(gene_index_path, args.hgnc_tsv)
    print(_report_line(report))

    ht = hl.read_table(args.ht)
    key_fields = list(ht.key)
    if len(key_fields) != 1:
        raise ValueError(
            f"{args.ht}: expected a single-field key (a gene-keyed table); got "
            f"{key_fields}"
        )
    key_field = key_fields[0]
    fields = [f.strip() for f in args.fields.split(",") if f.strip()]
    if not fields:
        raise ValueError("--fields must name at least one column")
    missing = [f for f in fields if f not in ht.row]
    if missing:
        raise ValueError(f"{args.ht}: no field(s) {missing}")

    df = ht.select(*fields).to_pandas().rename(columns={key_field: "gene_id"})
    df["gene"] = df["gene_id"].map(mapping)
    n_before = len(df)
    df = df[df["gene"].notna()].reset_index(drop=True)
    n_no_symbol = n_before - len(df)

    n_restricted = 0
    if args.restrict_to:
        if not args.restrict_column:
            raise ValueError("--restrict-to requires --restrict-column")
        restrict = _read_column_csv(args.restrict_to, args.restrict_column)
        n_before_restrict = len(df)
        df = df[df["gene"].isin(restrict)].reset_index(drop=True)
        n_restricted = n_before_restrict - len(df)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df[["gene", *fields]].to_csv(out_path, sep=_sep_for(out_path), index=False)
    print(
        f"features: {len(df)} gene(s) written ({n_no_symbol} dropped: no approved/"
        f"unambiguous symbol"
        + (f"; {n_restricted} dropped by --restrict-to" if args.restrict_to else "")
        + f") -> {out_path}"
    )
    return 0


# --- CLI -----------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="command", required=True)

    lb = sub.add_parser("labels", help="flatten a GeneSetCollection to a symbol list")
    lb.add_argument("--geneset", required=True)
    lb.add_argument("--group", action="append", default=None, metavar="NAME")
    lb.add_argument("--restrict-to", default=None)
    lb.add_argument("--restrict-column", default=None)
    lb.add_argument("--out", required=True)
    lb.set_defaults(fn=cmd_labels)

    co = sub.add_parser("cohort", help="map a cohort table's gene_id to a symbol")
    co.add_argument("--table", required=True)
    co.add_argument("--key-column", required=True)
    co.add_argument("--gene-index", required=True)
    co.add_argument("--hgnc-tsv", required=True)
    co.add_argument("--out", required=True)
    co.add_argument("--manifest-in", default=None)
    co.add_argument("--manifest-out", default=None)
    co.set_defaults(fn=cmd_cohort)

    fe = sub.add_parser("features", help="build a symbol-keyed feature axis TSV")
    fe.add_argument("--ht", required=True)
    fe.add_argument("--fields", required=True, help="comma-separated column names")
    fe.add_argument("--hgnc-tsv", required=True)
    fe.add_argument(
        "--gene-index",
        default=None,
        help="defaults to --ht; set separately when --ht itself carries no hgnc_id",
    )
    fe.add_argument("--restrict-to", default=None)
    fe.add_argument("--restrict-column", default=None)
    fe.add_argument("--out", required=True)
    fe.set_defaults(fn=cmd_features)

    args = ap.parse_args()
    try:
        return args.fn(args)
    except Exception as exc:  # noqa: BLE001 - a CLI entry point, not a library call
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
