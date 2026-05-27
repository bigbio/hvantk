"""Hail helpers shared between the eQTL (GTEx) and pQTL (Fang) Phase B builders.

GTEx and Fang use the same variant-ID format (``chr1_1000050_C_T_b38``) and
the same tissue-file scanning convention, so the parsing helpers live here
rather than being duplicated. Plugin-specific imports
(``_import_eqtl_gtex_*``, ``_import_pqtl_gtex_fang``) live in each plugin's
own ``builder.py``.

Previously colocated in ``hvantk/core/builders/table.py``; relocated when
that file was retired (issue #114).
"""

from __future__ import annotations

from pathlib import Path

import hail as hl


def parse_gtex_variant_id(ht, variant_id_field="variant_id", reference_genome="GRCh38"):
    """Parse GTEx variant IDs into ``locus`` and ``alleles``.

    Format: ``chr1_1000050_C_T_b38`` — the build suffix is discarded.
    Used by both eQTL and pQTL builders (GTEx/Fang share the same ID format).

    Contig names are normalised to match the reference genome: GRCh38 contigs
    use the ``chr`` prefix, GRCh37 contigs omit it.
    """
    parts = ht[variant_id_field].split("_")
    raw_contig = parts[0]
    bare = hl.if_else(raw_contig.startswith("chr"), raw_contig[3:], raw_contig)
    contig = hl.if_else(
        hl.literal(reference_genome).startswith("GRCh38"),
        "chr" + bare,
        bare,
    )
    return ht.annotate(
        locus=hl.locus(contig, hl.int32(parts[1]), reference_genome=reference_genome),
        alleles=hl.array([parts[2], parts[3]]),
    )


def strip_ensembl_version(gene_id_expr):
    """Strip Ensembl version suffix (``ENSG00000000003.15`` → ``ENSG00000000003``)."""
    return gene_id_expr.split("\\.")[0]


def scan_tissue_files(input_path, extensions):
    """Return ``[(file_path, tissue_name), ...]`` from *input_path*.

    Tissue name is inferred from the filename prefix before the first dot
    (e.g. ``Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz`` →
    ``Brain_Cortex``). Accepts a single file or a directory.
    """
    p = Path(input_path)
    if p.is_file():
        return [(str(p), p.stem.split(".")[0])]
    if not p.is_dir():
        raise FileNotFoundError(f"Not a file or directory: {input_path}")
    matches = []
    for ext in extensions:
        matches.extend(sorted(p.glob(f"*{ext}")))
    if not matches:
        raise FileNotFoundError(f"No files matching {extensions} in {input_path}")
    return [(str(f), f.stem.split(".")[0]) for f in matches]
