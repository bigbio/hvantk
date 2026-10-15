#!/usr/bin/env python3
"""Describe a built gene index and write its receipt.

Reports the gene count, the HGNC mapping rate, per-chromosome and per-biotype counts,
the genes with no HGNC record, and a content digest over every row and field -- the
checksum the receipt binds to, because the gene index itself carries no
``.provenance.json`` sidecar of its own (``references/inputs-and-release.md``).
Optionally records the two input datasets' own sidecar digests (their
``.provenance.json`` IS parsed, for a few identifying fields -- plugin, dataset,
schema_id, plugin_version, build_timestamp, source_fingerprint -- but never re-built
into a full Provenance object or re-validated against anything), so the receipt also
answers "what two inputs made this".

    python describe-gene-index.py --gene-index GI.ht [--gene-table G.ht] [--hgnc H.ht] \
        --receipt RECEIPT.json

Runs inside the image (needs Hail). Exits 1 on an empty table, a key other than
``gene_id``, or a missing ``SPINE_FIELDS`` name.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

#: Up to this many example genes with no HGNC record are recorded in the receipt.
N_HGNC_EXAMPLES = 10


def norm(value):
    """Stable-string floats, pass everything else through -- see
    ``skills/_conventions``' shared digest convention. The gene index's own fields are
    all strings, ints or arrays of strings (no floats), so this only matters if a future
    field adds one.
    """
    if isinstance(value, float):
        return format(value, ".10g")
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [norm(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(norm(v) for v in value)
    if isinstance(value, dict):
        return {k: norm(v) for k, v in value.items()}
    return value


def content_sha256(ht) -> str:
    """SHA-256 over every row and every field, order-insensitive.

    Deliberately not restricted to ``SPINE_FIELDS``: this is the receipt's own full
    fidelity digest, not the eval's contract-only digest (``eval/check.py`` has that
    one), so it also catches a field the contract does not yet know about.
    """
    from hvantk.tests._snapshot_utils import _to_jsonable

    rows = [norm(_to_jsonable(dict(row))) for row in ht.collect()]
    strings = sorted(json.dumps(r, sort_keys=True, separators=(",", ":")) for r in rows)
    return hashlib.sha256("\n".join(strings).encode()).hexdigest()


def sidecar_info(path_str) -> dict | None:
    """The sha256 of ``<path_str>.provenance.json`` plus a few identifying fields --
    just enough to bind the receipt to exactly that input build, never re-built into a
    full Provenance object. ``builder_commit`` is deliberately not one of these fields:
    it is always null when the input was built through ``hvantk reprocess`` (the normal
    path), so ``plugin_version``, ``build_timestamp`` and ``source_fingerprint`` are
    what actually distinguish one build from another here
    (``references/inputs-and-release.md``)."""
    if not path_str:
        return None
    sidecar = Path(str(path_str) + ".provenance.json")
    if not sidecar.is_file():
        return {"path": str(sidecar), "present": False}
    payload = json.loads(sidecar.read_text())
    return {
        "path": str(sidecar),
        "present": True,
        "sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
        "plugin": payload.get("plugin"),
        "dataset": payload.get("dataset"),
        "schema_id": payload.get("schema_id"),
        "plugin_version": payload.get("plugin_version"),
        "build_timestamp": payload.get("build_timestamp"),
        "source_fingerprint": payload.get("source_fingerprint"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-index", required=True)
    ap.add_argument("--gene-table", default=None, help="ensembl-gene:structure .ht")
    ap.add_argument("--hgnc", default=None, help="hgnc:lookup .ht")
    ap.add_argument("--receipt", required=True)
    args = ap.parse_args()

    from hvantk.algorithms.annotation.spine import SPINE_FIELDS, spine_mapping_rate
    from hvantk.core.utils.hail_context import init_hail

    init_hail(quiet=True)
    import hail as hl

    ht = hl.read_table(args.gene_index)

    if list(ht.key) != ["gene_id"]:
        print(f"ERROR: expected key ['gene_id'], got {list(ht.key)}", file=sys.stderr)
        return 1
    missing = [f for f in SPINE_FIELDS if f not in ht.row]
    if missing:
        print(
            f"ERROR: gene index at {args.gene_index} is missing SPINE_FIELDS: {missing}",
            file=sys.stderr,
        )
        return 1

    n_genes = ht.count()
    if n_genes == 0:
        print(f"ERROR: gene index at {args.gene_index} has zero rows", file=sys.stderr)
        return 1

    mapping_rate = spine_mapping_rate(ht)
    chromosome_counts = ht.aggregate(hl.agg.counter(ht.chromosome))
    biotype_counts = ht.aggregate(hl.agg.counter(ht.gene_biotype))
    n_without_hgnc = ht.aggregate(hl.agg.count_where(~hl.is_defined(ht.hgnc_id)))
    examples_rows = (
        ht.filter(~hl.is_defined(ht.hgnc_id))
        .select("gene_name")
        .head(N_HGNC_EXAMPLES)
        .collect()
    )
    hgnc_examples = [
        {"gene_id": row.gene_id, "gene_name": row.gene_name} for row in examples_rows
    ]
    digest = content_sha256(ht)
    chromosome_counts_sorted = dict(sorted(chromosome_counts.items()))
    biotype_counts_sorted = dict(sorted(biotype_counts.items()))
    gene_table_sidecar = sidecar_info(args.gene_table)
    hgnc_sidecar = sidecar_info(args.hgnc)

    receipt = {
        "gene_index": str(args.gene_index),
        "n_genes": n_genes,
        "mapping_rate": mapping_rate,
        "chromosome_counts": chromosome_counts_sorted,
        "biotype_counts": biotype_counts_sorted,
        "n_without_hgnc": n_without_hgnc,
        "missing_hgnc_examples": hgnc_examples,
        "content_sha256": digest,
        "gene_table_sidecar": gene_table_sidecar,
        "hgnc_sidecar": hgnc_sidecar,
    }
    Path(args.receipt).write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")

    print(f"gene index: {n_genes} genes -> {args.gene_index}")
    print(f"HGNC mapping rate: {mapping_rate:.4%} (full precision: {mapping_rate!r})")
    print(f"chromosome counts: {chromosome_counts_sorted}")
    print(f"biotype counts: {biotype_counts_sorted}")
    print(f"without HGNC: {n_without_hgnc}; examples: {hgnc_examples}")
    print(f"content_sha256: {digest}")
    print(f"gene-table sidecar: {gene_table_sidecar}")
    print(f"hgnc sidecar: {hgnc_sidecar}")
    print(f"receipt written: {args.receipt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
