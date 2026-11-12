#!/usr/bin/env python3
"""Describe a training set built by ``hvantk training-set build`` and write its receipt.

Reads the four outputs of one build -- ``<prefix>.ht`` (with its ``.provenance.json``
sidecar), ``<prefix>.tsv.gz``, ``<prefix>.dictionary.tsv`` and ``<prefix>.manifest.json``
-- and reports the row funnel, the label balance, the tier counts, the excluded columns,
the gene-list mapping and each input's identifying provenance. The receipt binds those
numbers to the SHA-256 of the ``.ht``'s provenance sidecar and of the ``.tsv.gz``.

It also checks that the four outputs agree with each other: the same rows in the
``.ht``, the TSV and the manifest; the same columns, in the same order, in the ``.ht``,
the TSV header and the dictionary; a tier on every training feature; the manifest's
tier counts equal to the dictionary's. Feature columns typed ``str`` are listed under
``string_features``: a score read from a dbNSFP built before rankscores were parsed as
numbers arrives as text without any error.

    python describe-training-set.py --prefix OUT/clinvar_heart \
        --receipt OUT/clinvar_heart.receipt.json

Runs inside the image (needs Hail: it counts the ``.ht``'s rows). Exits 1 when the
outputs disagree; ``string_features`` alone does not fail it.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
from pathlib import Path

TIERS = ("L1", "L2", "L3", "L4")
#: Identifying fields copied from each input's provenance in the manifest.
INPUT_FIELDS = ("dataset", "schema_id", "builder_commit", "source_fingerprint")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _tsv_header_and_rows(path: Path) -> tuple[list[str], int]:
    with gzip.open(path, "rt", encoding="utf-8", newline="\n") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = sum(1 for _ in fh)
    return header, rows


def describe(prefix: str) -> dict:
    """Every number the receipt carries, read from the four outputs of ``prefix``."""
    import hail as hl

    base = Path(prefix)
    ht_path = base.with_name(base.name + ".ht")
    tsv_path = base.with_name(base.name + ".tsv.gz")
    dict_path = base.with_name(base.name + ".dictionary.tsv")
    manifest_path = base.with_name(base.name + ".manifest.json")
    sidecar = ht_path.with_name(ht_path.name + ".provenance.json")

    ht = hl.read_table(str(ht_path))
    columns = list(ht.row)
    types = {name: str(ht[name].dtype) for name in columns}
    ht_rows = ht.count()

    tsv_header, tsv_rows = _tsv_header_and_rows(tsv_path)
    with open(dict_path, newline="", encoding="utf-8") as fh:
        dictionary = list(csv.DictReader(fh, delimiter="\t"))
    manifest = json.loads(manifest_path.read_text())

    features = [r for r in dictionary if r["training_feature"] == "yes"]
    tier_of = {r["column"]: r["tier"] for r in features}
    counts = manifest.get("counts", {})
    inputs = {}
    for name, prov in (manifest.get("input_provenance") or {}).items():
        inputs[name] = None if prov is None else {k: prov.get(k) for k in INPUT_FIELDS}
    return {
        "prefix": str(base.resolve()),
        "key": list(ht.key),
        "ht_rows": ht_rows,
        "ht_columns": columns,
        "tsv_header": tsv_header,
        "tsv_rows": tsv_rows,
        "tsv_sha256": sha256_file(tsv_path),
        "sidecar_sha256": sha256_file(sidecar) if sidecar.exists() else None,
        "schema_id": json.loads(sidecar.read_text()).get("schema_id")
        if sidecar.exists()
        else None,
        "dictionary_columns": [r["column"] for r in dictionary],
        "features": len(features),
        "features_without_tier": [c for c, t in tier_of.items() if not t],
        "dictionary_tier_counts": {
            t: sum(v == t for v in tier_of.values()) for t in TIERS
        },
        "string_features": sorted(c for c in tier_of if types.get(c) == "str"),
        "spec": manifest.get("spec"),
        "funnel": {
            k: counts.get(k)
            for k in (
                "clinvar_rows",
                "after_variant_filter",
                "labelled",
                "labelled_min_stars",
                "rows",
            )
        },
        "rows": counts.get("rows"),
        "label_counts": counts.get("label", {}),
        "label_by_gene_list": counts.get("label_by_gene_list", {}),
        "review_stars": counts.get("review_stars", {}),
        "no_gene_id": counts.get("no_gene_id"),
        "gene_rule": counts.get("gene_rule", {}),
        "tier_counts": manifest.get("tiers", {}),
        "excluded": manifest.get("excluded", []),
        "gene_lists": {
            name: {k: v for k, v in report.items() if k != "unmapped"}
            for name, report in (manifest.get("gene_lists") or {}).items()
        },
        "inputs": inputs,
    }


def consistency_problems(desc: dict) -> list[str]:
    """Every disagreement between the four outputs; empty means consistent."""
    rows = desc["ht_rows"]
    manifest_tiers = {t: desc["tier_counts"].get(t, 0) for t in TIERS}
    checks = [
        (desc["sidecar_sha256"] is None, "the .ht has no .provenance.json sidecar"),
        (
            desc["key"] != ["locus", "alleles"],
            f".ht key is {desc['key']}, not ['locus', 'alleles']",
        ),
        (rows == 0, "the .ht has no rows"),
        (desc["tsv_rows"] != rows, f"TSV rows {desc['tsv_rows']} != .ht rows {rows}"),
        (desc["rows"] != rows, f"manifest rows {desc['rows']} != .ht rows {rows}"),
        (
            sum(desc["label_counts"].values()) != rows,
            f"manifest label counts {desc['label_counts']} do not sum to the .ht rows",
        ),
        (
            desc["tsv_header"] != desc["ht_columns"],
            "TSV header differs from the .ht columns (names or order)",
        ),
        (
            desc["dictionary_columns"] != desc["ht_columns"],
            "dictionary columns differ from the .ht columns (names or order)",
        ),
        (
            bool(desc["features_without_tier"]),
            "training features without a tier: "
            + ", ".join(desc["features_without_tier"]),
        ),
        (
            manifest_tiers != desc["dictionary_tier_counts"],
            f"manifest tier counts {manifest_tiers} != dictionary tier counts "
            f"{desc['dictionary_tier_counts']}",
        ),
    ]
    return [message for failed, message in checks if failed]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prefix", required=True, help="the build's --output prefix")
    ap.add_argument("--receipt", default=None, help="write the receipt JSON here")
    args = ap.parse_args()

    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    desc = describe(args.prefix)
    problems = consistency_problems(desc)

    print(f"prefix:      {desc['prefix']}")
    print(f"rows:        {desc['ht_rows']} (funnel {desc['funnel']})")
    print(f"labels:      {desc['label_counts']}")
    print(f"tiers:       {desc['tier_counts']} of {desc['features']} features")
    print(f"excluded:    {desc['excluded'] or 'none'}")
    print(f"no_gene_id:  {desc['no_gene_id']}")
    for name, report in desc["gene_lists"].items():
        print(f"gene list:   {name} {report}")
    print(f"str features: {desc['string_features'] or 'none'}")
    print(f"tsv sha256:  {desc['tsv_sha256']}")
    print(f"sidecar:     {desc['sidecar_sha256'] or 'MISSING'}")
    for p in problems:
        print(f"PROBLEM: {p}")

    if args.receipt:
        receipt = {
            **{k: v for k, v in desc.items() if k not in ("ht_columns", "tsv_header")},
            "columns": len(desc["ht_columns"]),
            "consistent": not problems,
            "problems": problems,
            "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
        Path(args.receipt).write_text(json.dumps(receipt, indent=2) + "\n")
        print(f"receipt:     {args.receipt}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
