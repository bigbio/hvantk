"""Plugin-local constants for the gnomad-metrics skill.

gnomAD constraint download coordinates. The tables live in the public gnomAD
release bucket (``gcp-public-data--gnomad``) over plain HTTPS, no auth.

Two releases are supported by the downloader:

- ``v2.1.1`` (GRCh37) — flat column names, one row per gene, includes a
  ``gene_id`` column. ``by_gene`` is the table hvantk standardises on
  (pLI / oe_lof / LOEUF / mis_z, keyed by ``gene_id``). ``by_transcript`` is
  the same metrics at transcript resolution.
- ``v4.0`` (GRCh38) — one row per transcript (with a ``mane_select`` flag),
  **dotted** column names (``lof.pLI``, ``lof.oe_ci.upper`` = LOEUF,
  ``mis.z_score`` …) and **no ``gene_id`` column**. Building it needs a
  non-default key (see ``build_gnomad_metrics_metrics`` ``key`` param).

Note: the v2.1.1 path segment is ``2.1.1`` (no ``v`` prefix); v4.0 is ``v4.0``.
gnomAD did not re-release constraint for v4.1, so v4.0 is the newest.
"""

# Public gnomAD release bucket (HTTPS, no auth).
GNOMAD_RELEASE_BASE_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release"
)

# version -> table -> object path (relative to GNOMAD_RELEASE_BASE_URL).
GNOMAD_CONSTRAINT_TABLES = {
    "v2.1.1": {
        "by_gene": "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        "by_transcript": (
            "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
        ),
    },
    "v4.0": {
        "constraint_metrics": ("v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv"),
    },
}

# Default release + per-release default table.
DEFAULT_VERSION = "v2.1.1"
DEFAULT_TABLE = {
    "v2.1.1": "by_gene",
    "v4.0": "constraint_metrics",
}


def resolve_table(version: str, table: str | None = None) -> str:
    """Validate ``version``/``table`` and return the resolved table name.

    Raises
    ------
    ValueError
        If ``version`` is unknown, or ``table`` is not valid for ``version``.
    """
    if version not in GNOMAD_CONSTRAINT_TABLES:
        raise ValueError(
            f"Unknown gnomAD constraint version {version!r}; "
            f"choose from {sorted(GNOMAD_CONSTRAINT_TABLES)}"
        )
    tables = GNOMAD_CONSTRAINT_TABLES[version]
    if table is None:
        table = DEFAULT_TABLE[version]
    if table not in tables:
        raise ValueError(
            f"Unknown table {table!r} for gnomAD {version}; "
            f"choose from {sorted(tables)}"
        )
    return table


def constraint_url(version: str, table: str | None = None) -> str:
    """Return the download URL for a gnomAD constraint table."""
    table = resolve_table(version, table)
    return f"{GNOMAD_RELEASE_BASE_URL}/{GNOMAD_CONSTRAINT_TABLES[version][table]}"


def constraint_filename(version: str, table: str | None = None) -> str:
    """Return the canonical basename for a gnomAD constraint table."""
    return constraint_url(version, table).rsplit("/", 1)[-1]
