"""
MONDO Disease Ontology parser and utilities.

Provides tools for parsing the MONDO ontology and categorizing diseases
based on their ontological relationships.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

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


class MondoOntology:
    """
    Parser and utility class for the MONDO disease ontology.

    Provides methods to:
    - Parse MONDO OBO format files
    - Navigate the disease hierarchy
    - Categorize diseases based on their ancestors
    """

    def __init__(self, obo_path: Optional[str] = None):
        """
        Initialize the MONDO ontology parser.

        Parameters
        ----------
        obo_path : str, optional
            Path to the MONDO OBO file. If not provided, must call load() later.
        """
        self.terms: Dict[str, dict] = {}  # id -> {name, is_a, ...}
        self.children: Dict[str, Set[str]] = defaultdict(set)  # parent -> children
        self.name_to_id: Dict[str, str] = {}  # name -> id (lowercase)
        self._loaded = False

        if obo_path:
            self.load(obo_path)

    def load(self, obo_path: str) -> None:
        """
        Load and parse a MONDO OBO file.

        Parameters
        ----------
        obo_path : str
            Path to the MONDO OBO file.
        """
        logger.info(f"Loading MONDO ontology from {obo_path}")

        self.terms = {}
        self.children = defaultdict(set)
        self.name_to_id = {}

        current_term = None
        in_term_section = False

        with open(obo_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    # Save previous term if any
                    if current_term and "id" in current_term:
                        self._add_term(current_term)
                    current_term = {"is_a": []}
                    in_term_section = True
                elif line.startswith("[") and line.endswith("]"):
                    # Other section (Typedef, etc.) - save current term and skip
                    if current_term and "id" in current_term:
                        self._add_term(current_term)
                    current_term = None
                    in_term_section = False
                elif in_term_section and current_term is not None:
                    if line.startswith("id: "):
                        current_term["id"] = line[4:].strip()
                    elif line.startswith("name: "):
                        current_term["name"] = line[6:].strip()
                    elif line.startswith("is_a: "):
                        # Parse "is_a: MONDO:0000001 ! disease"
                        parent_part = line[6:].split("!")[0].strip()
                        # Handle annotations like {source="..."}
                        parent_id = parent_part.split()[0].strip()
                        if parent_id.startswith("MONDO:"):
                            current_term["is_a"].append(parent_id)
                    elif line.startswith("is_obsolete: true"):
                        current_term["obsolete"] = True

        # Don't forget the last term
        if current_term and "id" in current_term:
            self._add_term(current_term)

        self._loaded = True
        logger.info(f"Loaded {len(self.terms)} MONDO terms")

    def _add_term(self, term: dict) -> None:
        """Add a parsed term to the internal data structures."""
        term_id = term.get("id")
        if not term_id or not term_id.startswith("MONDO:"):
            return

        # Skip obsolete terms
        if term.get("obsolete"):
            return

        self.terms[term_id] = term

        # Build parent-child relationships
        for parent_id in term.get("is_a", []):
            self.children[parent_id].add(term_id)

        # Build name lookup
        name = term.get("name", "").lower()
        if name:
            self.name_to_id[name] = term_id

    def get_name(self, mondo_id: str) -> Optional[str]:
        """Get the name for a MONDO ID."""
        term = self.terms.get(mondo_id)
        return term.get("name") if term else None

    def get_parents(self, mondo_id: str) -> List[str]:
        """Get direct parent IDs for a term."""
        term = self.terms.get(mondo_id)
        return term.get("is_a", []) if term else []

    def get_children(self, mondo_id: str) -> Set[str]:
        """Get direct children IDs for a term."""
        return self.children.get(mondo_id, set())

    def get_ancestors(self, mondo_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all ancestors (parents, grandparents, etc.) of a term.

        Parameters
        ----------
        mondo_id : str
            The MONDO ID to find ancestors for.
        include_self : bool
            Whether to include the term itself in the result.

        Returns
        -------
        set
            Set of ancestor MONDO IDs.
        """
        ancestors = set()
        if include_self:
            ancestors.add(mondo_id)

        to_visit = list(self.get_parents(mondo_id))
        while to_visit:
            parent = to_visit.pop()
            if parent not in ancestors:
                ancestors.add(parent)
                to_visit.extend(self.get_parents(parent))

        return ancestors

    def get_descendants(self, mondo_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all descendants (children, grandchildren, etc.) of a term.

        Parameters
        ----------
        mondo_id : str
            The MONDO ID to find descendants for.
        include_self : bool
            Whether to include the term itself in the result.

        Returns
        -------
        set
            Set of descendant MONDO IDs.
        """
        descendants = set()
        if include_self:
            descendants.add(mondo_id)

        to_visit = list(self.get_children(mondo_id))
        while to_visit:
            child = to_visit.pop()
            if child not in descendants:
                descendants.add(child)
                to_visit.extend(self.get_children(child))

        return descendants

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

        # Normalize the MONDO ID
        mondo_id = self._normalize_mondo_id(mondo_id)

        # Get all ancestors of this term
        ancestors = self.get_ancestors(mondo_id, include_self=True)

        # Find which categories this term falls under
        matches = []
        for cat_id, cat_name in categories.items():
            if cat_id in ancestors:
                matches.append((cat_id, cat_name))

        return matches

    def categorize_multiple(
        self,
        mondo_ids: List[str],
        categories: Optional[Dict[str, str]] = None,
    ) -> Dict[str, List[str]]:
        """
        Categorize multiple diseases and return genes grouped by category.

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

        result = defaultdict(list)

        for mondo_id in mondo_ids:
            cats = self.categorize(mondo_id, categories)
            for cat_id, cat_name in cats:
                result[cat_name].append(mondo_id)

        return dict(result)

    def _normalize_mondo_id(self, mondo_id: str) -> str:
        """Normalize a MONDO ID to the standard format."""
        mondo_id = mondo_id.strip()
        if mondo_id.startswith("MONDO:"):
            return mondo_id
        # Handle numeric IDs without prefix
        if mondo_id.isdigit():
            return f"MONDO:{mondo_id.zfill(7)}"
        return f"MONDO:{mondo_id}"

    def search_by_name(self, query: str, exact: bool = False) -> List[Tuple[str, str]]:
        """
        Search for terms by name.

        Parameters
        ----------
        query : str
            Search query.
        exact : bool
            If True, require exact match. Otherwise, substring match.

        Returns
        -------
        list of tuples
            List of (mondo_id, name) tuples matching the query.
        """
        query = query.lower()
        results = []

        for term_id, term in self.terms.items():
            name = term.get("name", "").lower()
            if exact:
                if name == query:
                    results.append((term_id, term.get("name")))
            else:
                if query in name:
                    results.append((term_id, term.get("name")))

        return results


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

    url = "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo"
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
            with open(output_path, 'wb') as f:
                f.write(response.read())

    logger.info(f"Downloaded MONDO ontology to {output_path}")
    return str(output_path)
