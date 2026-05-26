"""
MONDO Disease Ontology parser and utilities.

Provides tools for parsing the MONDO ontology and categorizing diseases
based on their ontological relationships. Extends the generic OBO parser
with MONDO-specific categories and download logic.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from hvantk.core.utils.obo_parser import BaseOboOntology

logger = logging.getLogger(__name__)


# Top-level disease categories in MONDO for clinical relevance
# These are high-level groupings based on organ systems, etiology, etc.
MONDO_DISEASE_CATEGORIES = {
    # Organ/System-based categories
    "MONDO:0004995": "cardiovascular disease",
    "MONDO:0005071": "nervous system disease",
    "MONDO:0005066": "metabolic disease",
    "MONDO:0005046": "immune system disease",
    "MONDO:0005328": "eye disease",
    "MONDO:0005087": "respiratory system disease",
    "MONDO:0005151": "endocrine system disease",
    "MONDO:0002025": "psychiatric disorder",
    "MONDO:0005042": "gastrointestinal system disease",
    "MONDO:0100545": "hereditary neurological disease",
    "MONDO:0021166": "auditory system disease",
    "MONDO:0002081": "musculoskeletal system disease",
    "MONDO:0024623": "skin disease",
    "MONDO:0005240": "kidney disease",
    "MONDO:0005154": "liver disease",
    "MONDO:0003150": "bone disease",
    "MONDO:0002051": "hematopoietic system disease",
    "MONDO:0005039": "reproductive system disease",
    # Cancer/neoplasm
    "MONDO:0005070": "neoplasm",
    "MONDO:0004992": "cancer",
    # Genetic/hereditary
    "MONDO:0003847": "hereditary disease",
    "MONDO:0000426": "autosomal dominant disease",
    "MONDO:0006025": "autosomal recessive disease",
    "MONDO:0010383": "X-linked disease",
    # Developmental
    "MONDO:0021147": "developmental disorder",
    "MONDO:0002320": "congenital nervous system disorder",
    "MONDO:0019042": "malformation syndrome",
    "MONDO:0100500": "neurodevelopmental disorder",
    # Infectious
    "MONDO:0005550": "infectious disease",
    # Syndromic
    "MONDO:0002254": "syndromic disease",
    # Connective tissue
    "MONDO:0023603": "hereditary connective tissue disorder",
    # Cardiogenetic
    "MONDO:0100547": "cardiogenetic disease",
}


class MondoOntology(BaseOboOntology):
    """
    Parser and utility class for the MONDO disease ontology.

    Extends :class:`BaseOboOntology` with MONDO-specific defaults:
    only ``MONDO:`` terms are loaded, and ``categorize`` / ``categorize_multiple``
    default to :data:`MONDO_DISEASE_CATEGORIES` when no custom categories
    are provided.

    Parameters
    ----------
    obo_path : str, optional
        Path to the MONDO OBO file. If not provided, must call ``load()`` later.
    """

    def __init__(self, obo_path: Optional[str] = None):
        super().__init__(obo_path=obo_path, id_prefix="MONDO:")

    def categorize(
        self,
        mondo_id: str,
        categories: Optional[Dict[str, str]] = None,
    ) -> List[Tuple[str, str]]:
        """
        Categorize a disease by finding which top-level categories it belongs to.

        Parameters
        ----------
        mondo_id : str
            The MONDO ID to categorize.
        categories : dict, optional
            Custom category mapping {MONDO_ID: category_name}.
            If None, uses MONDO_DISEASE_CATEGORIES.

        Returns
        -------
        list of tuples
            List of (category_id, category_name) tuples for matching categories.
        """
        if categories is None:
            categories = MONDO_DISEASE_CATEGORIES

        mondo_id = self._normalize_mondo_id(mondo_id)
        return super().categorize(mondo_id, categories)

    def categorize_multiple(
        self,
        mondo_ids: List[str],
        categories: Optional[Dict[str, str]] = None,
    ) -> Dict[str, List[str]]:
        """
        Categorize multiple diseases and return them grouped by category.

        Parameters
        ----------
        mondo_ids : list
            List of MONDO IDs to categorize.
        categories : dict, optional
            Custom category mapping {MONDO_ID: category_name}.

        Returns
        -------
        dict
            Mapping of category_name -> list of MONDO IDs in that category.
        """
        if categories is None:
            categories = MONDO_DISEASE_CATEGORIES

        return super().categorize_multiple(mondo_ids, categories)

    def _normalize_mondo_id(self, mondo_id: str) -> str:
        """Normalize a MONDO ID to the standard format."""
        mondo_id = mondo_id.strip()
        if mondo_id.startswith("MONDO:"):
            return mondo_id
        # Handle numeric IDs without prefix
        if mondo_id.isdigit():
            return f"MONDO:{mondo_id.zfill(7)}"
        return f"MONDO:{mondo_id}"


def download_mondo_obo(output_path: str, overwrite: bool = False) -> str:
    """
    Download the latest MONDO OBO file.

    Parameters
    ----------
    output_path : str
        Path to save the downloaded file.
    overwrite : bool
        Whether to overwrite existing file.

    Returns
    -------
    str
        Path to the downloaded file.
    """
    import urllib.request
    import ssl

    output_path = Path(output_path)

    if output_path.exists() and not overwrite:
        logger.info(f"MONDO OBO file already exists at {output_path}")
        return str(output_path)

    url = (
        "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo"
    )
    logger.info(f"Downloading MONDO ontology from {url}")

    # Create SSL context that doesn't verify (for environments with cert issues)
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE

    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        urllib.request.urlretrieve(url, output_path)
    except Exception:
        # Try with unverified SSL
        with urllib.request.urlopen(url, context=ssl_context) as response:
            with open(output_path, "wb") as f:
                f.write(response.read())

    logger.info(f"Downloaded MONDO ontology to {output_path}")
    return str(output_path)
