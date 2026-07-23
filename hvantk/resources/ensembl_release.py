"""The single Ensembl release pin for the whole toolkit.

Everything that reads Ensembl gene models -- the ``ensembl-gene`` plugin, the PTM
coordinate mapper, and the gene spine -- imports the release from here. Two definitions
that drift apart give different CDS lengths, different MANE Select assignments and
different coordinates for the same gene, which is a silent correctness bug rather than a
loud one.

It lives in ``resources`` rather than in the plugin because the dependency rule is
``tools -> skills -> algorithms -> core``, with ``resources`` as substrate alongside
``core``. Both a skill (``ensembl_gene``) and an algorithm (``ptm``) consume this pin;
putting it in the skill would make ``algorithms`` depend on ``skills`` and invert the
rule.
"""
from __future__ import annotations

ENSEMBL_RELEASE = "113"

ENSEMBL_GTF_FILENAME = f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"

ENSEMBL_GTF_URL = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}"
    f"/gtf/homo_sapiens/{ENSEMBL_GTF_FILENAME}"
)
