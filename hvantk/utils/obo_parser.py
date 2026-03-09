"""
Generic OBO format ontology parser.

Provides a base class for parsing OBO ontology files and navigating
term hierarchies. Ontology-specific subclasses (e.g., MondoOntology)
extend this with domain-specific categories and utilities.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class BaseOboOntology:
    """
    Generic parser for OBO format ontology files.

    Parses standard OBO format (``[Term]`` sections with ``id:``, ``name:``,
    ``is_a:`` relationships) and provides hierarchy traversal methods.

    Parameters
    ----------
    obo_path : str, optional
        Path to an OBO file. If provided, the file is loaded immediately.
    id_prefix : str, optional
        If set, only terms whose ID starts with this prefix are kept
        (e.g., ``"MONDO:"``, ``"HP:"``, ``"GO:"``). Terms with other
        prefixes are silently skipped during parsing. If None, all terms
        are loaded.
    """

    def __init__(
        self,
        obo_path: Optional[str] = None,
        id_prefix: Optional[str] = None,
    ):
        self.id_prefix = id_prefix
        self.terms: Dict[str, dict] = {}  # id -> {name, is_a, ...}
        self.children: Dict[str, Set[str]] = defaultdict(set)
        self.name_to_id: Dict[str, str] = {}  # lowercase name -> id
        self._loaded = False

        if obo_path:
            self.load(obo_path)

    @property
    def loaded(self) -> bool:
        """Whether an OBO file has been loaded."""
        return self._loaded

    def load(self, obo_path: str) -> None:
        """
        Load and parse an OBO file.

        Parameters
        ----------
        obo_path : str
            Path to the OBO file.
        """
        logger.info(f"Loading ontology from {obo_path}")

        self.terms = {}
        self.children = defaultdict(set)
        self.name_to_id = {}

        current_term = None
        in_term_section = False

        with open(obo_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    if current_term and "id" in current_term:
                        self._add_term(current_term)
                    current_term = {"is_a": []}
                    in_term_section = True
                elif line.startswith("[") and line.endswith("]"):
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
                        # Parse "is_a: PREFIX:0000001 ! some term"
                        parent_part = line[6:].split("!")[0].strip()
                        parent_id = parent_part.split()[0].strip()
                        current_term["is_a"].append(parent_id)
                    elif line.startswith("is_obsolete: true"):
                        current_term["obsolete"] = True

        # Don't forget the last term
        if current_term and "id" in current_term:
            self._add_term(current_term)

        self._loaded = True
        logger.info(f"Loaded {len(self.terms)} terms")

    def _add_term(self, term: dict) -> None:
        """Add a parsed term to the internal data structures."""
        term_id = term.get("id")
        if not term_id:
            return

        # Apply prefix filter if configured
        if self.id_prefix and not term_id.startswith(self.id_prefix):
            return

        # Skip obsolete terms
        if term.get("obsolete"):
            return

        # Filter is_a parents by prefix if configured
        if self.id_prefix:
            term["is_a"] = [
                p for p in term.get("is_a", [])
                if p.startswith(self.id_prefix)
            ]

        self.terms[term_id] = term

        for parent_id in term.get("is_a", []):
            self.children[parent_id].add(term_id)

        name = term.get("name", "").lower()
        if name:
            self.name_to_id[name] = term_id

    def get_name(self, term_id: str) -> Optional[str]:
        """Get the name for a term ID."""
        term = self.terms.get(term_id)
        return term.get("name") if term else None

    def get_parents(self, term_id: str) -> List[str]:
        """Get direct parent IDs for a term."""
        term = self.terms.get(term_id)
        return term.get("is_a", []) if term else []

    def get_children(self, term_id: str) -> Set[str]:
        """Get direct children IDs for a term."""
        return self.children.get(term_id, set())

    def get_ancestors(self, term_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all ancestors (parents, grandparents, etc.) of a term.

        Parameters
        ----------
        term_id : str
            The term ID to find ancestors for.
        include_self : bool
            Whether to include the term itself in the result.

        Returns
        -------
        set
            Set of ancestor term IDs.
        """
        ancestors = set()
        if include_self:
            ancestors.add(term_id)

        to_visit = list(self.get_parents(term_id))
        while to_visit:
            parent = to_visit.pop()
            if parent not in ancestors:
                ancestors.add(parent)
                to_visit.extend(self.get_parents(parent))

        return ancestors

    def get_descendants(self, term_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all descendants (children, grandchildren, etc.) of a term.

        Parameters
        ----------
        term_id : str
            The term ID to find descendants for.
        include_self : bool
            Whether to include the term itself in the result.

        Returns
        -------
        set
            Set of descendant term IDs.
        """
        descendants = set()
        if include_self:
            descendants.add(term_id)

        to_visit = list(self.get_children(term_id))
        while to_visit:
            child = to_visit.pop()
            if child not in descendants:
                descendants.add(child)
                to_visit.extend(self.get_children(child))

        return descendants

    def categorize(
        self,
        term_id: str,
        categories: Dict[str, str],
    ) -> List[Tuple[str, str]]:
        """
        Categorize a term by finding which category ancestors it belongs to.

        Parameters
        ----------
        term_id : str
            The term ID to categorize.
        categories : dict
            Category mapping ``{term_id: category_name}``.

        Returns
        -------
        list of tuples
            List of ``(category_id, category_name)`` for matching categories.
        """
        ancestors = self.get_ancestors(term_id, include_self=True)

        matches = []
        for cat_id, cat_name in categories.items():
            if cat_id in ancestors:
                matches.append((cat_id, cat_name))

        return matches

    def categorize_multiple(
        self,
        term_ids: List[str],
        categories: Dict[str, str],
    ) -> Dict[str, List[str]]:
        """
        Categorize multiple terms and group them by category.

        Parameters
        ----------
        term_ids : list
            List of term IDs to categorize.
        categories : dict
            Category mapping ``{term_id: category_name}``.

        Returns
        -------
        dict
            Mapping of ``category_name`` -> list of term IDs in that category.
        """
        result = defaultdict(list)

        for tid in term_ids:
            cats = self.categorize(tid, categories)
            for cat_id, cat_name in cats:
                result[cat_name].append(tid)

        return dict(result)

    def search_by_name(
        self, query: str, exact: bool = False
    ) -> List[Tuple[str, str]]:
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
            List of ``(term_id, name)`` tuples matching the query.
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
