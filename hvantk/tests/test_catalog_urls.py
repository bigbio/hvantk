"""Every catalog-declared download URL must be pinned, not a moving "latest" alias.

A ``/latest/`` path segment resolves to whatever the upstream currently calls latest,
which can change bytes under a committed accession with no record of what changed.
#380 fixed exactly this for ``gwas_catalog`` (pinned to its dated release archive) and
marked ``msigdb``'s URL provenance-only; nothing asserted the rule, so a future catalog
edit could reintroduce it with no test noticing.

Reads the catalogs directly off disk rather than through ``HvantkRegistry``: the
registry's per-bucket views are keyed and reshaped for search/listing, and a ``url``
can be nested (``files[].url``) as well as top-level on an entry, so a plain
recursive walk of the raw JSON is the simplest way to see every occurrence regardless
of where the schema puts it. Stdlib-only, no Hail, so this runs in the default pytest
selection.
"""

from __future__ import annotations

import json
from pathlib import Path

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"

CATALOG_PATHS = sorted(SKILLS_DIR.glob("*/catalog/datasets.json"))


def _iter_urls(node):
    """Yield every string value bound to a ``url`` key anywhere under ``node``.

    Recurses through dicts and lists alike, since a ``url`` can sit directly on a
    dataset entry or nested one level down inside its ``files`` list.
    """
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "url" and isinstance(value, str):
                yield value
            else:
                yield from _iter_urls(value)
    elif isinstance(node, list):
        for item in node:
            yield from _iter_urls(item)


def test_discovery_finds_the_known_catalogs():
    """Guard the guard: if the glob silently found nothing, the check below would
    pass for the wrong reason."""
    names = {p.parent.parent.name for p in CATALOG_PATHS}
    assert "gwas_catalog" in names, (
        f"expected gwas_catalog's catalog among the discovered ones, found: {names}"
    )
    assert len(CATALOG_PATHS) >= 10, (
        f"expected at least 10 catalog/datasets.json files, found {len(CATALOG_PATHS)}"
    )


def test_the_scanner_sees_a_known_url():
    """Guard the guard: if `_iter_urls` silently found nothing, the check below would
    pass for the wrong reason."""
    gwas = next(p for p in CATALOG_PATHS if p.parent.parent.name == "gwas_catalog")
    urls = list(_iter_urls(json.loads(gwas.read_text())))
    assert any("ftp.ebi.ac.uk" in u for u in urls), urls


def test_no_catalog_url_points_at_a_moving_latest_alias():
    """Checked as a path SEGMENT (``/latest/``), not a bare substring: the clinvar
    catalog's ``ClinVar_latest`` is an accession, not a URL, and gwas_catalog's own
    description prose names the ``releases/latest/`` endpoint it deliberately did NOT
    pin to (#380) -- neither should trip this.
    """
    offenders: list[str] = []
    for catalog_path in CATALOG_PATHS:
        data = json.loads(catalog_path.read_text())
        for url in _iter_urls(data):
            if "/latest/" in url:
                offenders.append(f"{catalog_path.relative_to(SKILLS_DIR)}: {url}")
    assert not offenders, (
        "catalog url pinned to a moving /latest/ alias instead of a specific release "
        "(see #380):\n  " + "\n  ".join(offenders)
    )
