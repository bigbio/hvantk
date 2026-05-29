"""Gene set input parsing and validation for prepare-geneset.

Parses a headerless two-column TSV (gene_set_name<TAB>gene_symbol) into a
dictionary of {set_name: List[str]}, with format and identifier validation.

Also provides HGNC-based symbol validation and alias resolution using the
existing ``_load_hgnc_symbol_maps()`` infrastructure from
``hvantk.core.utils.gene_aliases``.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import logging
import re

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ID type detection
# ---------------------------------------------------------------------------

_ENSEMBL_GENE_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$")
_ENSEMBL_TRANSCRIPT_RE = re.compile(r"^ENST\d{11}(\.\d+)?$")
_ENTREZ_RE = re.compile(r"^\d+$")
_MOUSE_SYMBOL_RE = re.compile(r"^[A-Z][a-z][a-z0-9]+$")


def detect_id_type(token: str) -> str:
    """Classify a gene identifier token.

    Returns one of: ``"ensembl_gene"``, ``"ensembl_transcript"``,
    ``"entrez"``, ``"mouse_symbol"``, ``"whitespace"``, ``"symbol"``.
    """
    if " " in token or "\t" in token:
        return "whitespace"
    if _ENSEMBL_GENE_RE.match(token):
        return "ensembl_gene"
    if _ENSEMBL_TRANSCRIPT_RE.match(token):
        return "ensembl_transcript"
    if _ENTREZ_RE.match(token):
        return "entrez"
    if _MOUSE_SYMBOL_RE.match(token) and not re.search(r"[A-Z]{2}", token):
        return "mouse_symbol"
    return "symbol"


def validate_gene_ids(
    genes: List[str],
) -> Tuple[List[str], Dict[str, List[str]]]:
    """Check a list of gene tokens for non-symbol identifiers.

    Returns
    -------
    valid : List[str]
        Tokens classified as ``"symbol"``.
    problems : Dict[str, List[str]]
        Mapping of problem type to list of offending tokens.
    """
    valid: List[str] = []
    problems: Dict[str, List[str]] = {}

    for gene in genes:
        id_type = detect_id_type(gene)
        if id_type == "symbol":
            valid.append(gene)
        else:
            problems.setdefault(id_type, []).append(gene)

    return valid, problems


_ERROR_TEMPLATES = {
    "ensembl_gene": (
        "Found Ensembl gene IDs (e.g., {example}). This command requires "
        "HGNC gene symbols. Convert first using biomart, hvantk GeneMapper, "
        "or https://www.genenames.org/tools/multi-symbol-checker/"
    ),
    "ensembl_transcript": (
        "Found Ensembl transcript IDs (e.g., {example}). This command "
        "requires HGNC gene symbols, not transcript IDs. Convert first "
        "using biomart, hvantk GeneMapper, or "
        "https://www.genenames.org/tools/multi-symbol-checker/"
    ),
    "entrez": (
        "Found numeric-only entries (e.g., {example}) that may be Entrez "
        "Gene IDs. This command requires HGNC gene symbols."
    ),
    "whitespace": (
        "Gene symbol contains whitespace: '{example}'. Check input formatting."
    ),
}


def _build_error_message(problems: Dict[str, List[str]]) -> str:
    """Build a human-readable error message from validation problems."""
    parts: List[str] = []
    for ptype, tokens in problems.items():
        if ptype == "mouse_symbol":
            # Mouse symbols use >50% threshold; handled in parse_geneset_tsv.
            continue
        template = _ERROR_TEMPLATES.get(ptype, "Unknown problem type: {example}")
        parts.append(template.format(example=tokens[0]))
    return " ".join(parts)


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------


@dataclass
class ParseResult:
    """Result of parsing a gene set TSV file."""

    gene_sets: Dict[str, List[str]]
    n_lines_parsed: int = 0
    n_lines_skipped: int = 0
    n_duplicates: int = 0
    warnings: List[str] = field(default_factory=list)


def parse_geneset_tsv(path: Path) -> ParseResult:
    """Parse a headerless two-column TSV into gene sets.

    Parameters
    ----------
    path : Path
        Input file. Expected format: ``gene_set_name<TAB>gene_symbol``,
        one line per membership. ``#`` comments and blank lines ignored.

    Returns
    -------
    ParseResult

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file has wrong number of columns on any data line, or
        contains non-HGNC identifiers (Ensembl, Entrez, mouse symbols).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    gene_sets: Dict[str, List[str]] = {}
    seen: Dict[str, Set[str]] = {}  # set_name → seen genes (for dedup)
    n_lines_parsed = 0
    n_lines_skipped = 0
    n_duplicates = 0
    warnings: List[str] = []

    with open(path) as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()

            # Skip blanks and comments.
            if not line or line.startswith("#"):
                n_lines_skipped += 1
                continue

            n_lines_parsed += 1
            fields = line.split("\t")

            if len(fields) != 2:
                raise ValueError(
                    f"Line {lineno}: expected 2 tab-separated columns, "
                    f"got {len(fields)}: {line!r}"
                )

            set_name, gene = fields[0].strip(), fields[1].strip()

            # Empty gene symbol → warn and skip.
            if not gene:
                warnings.append(f"Line {lineno}: empty gene symbol, skipping.")
                continue

            # Duplicate detection within the same set.
            if set_name not in seen:
                seen[set_name] = set()
                gene_sets[set_name] = []

            if gene in seen[set_name]:
                n_duplicates += 1
                warnings.append(
                    f"Duplicate gene '{gene}' in set '{set_name}', "
                    "keeping first occurrence."
                )
                continue

            seen[set_name].add(gene)
            gene_sets[set_name].append(gene)

    if not gene_sets:
        raise ValueError(f"No gene set entries found in {path}")

    # --- Validate gene identifiers across all unique symbols ---
    all_genes = list({g for genes in gene_sets.values() for g in genes})
    _, problems = validate_gene_ids(all_genes)

    # Hard errors for ensembl, entrez, whitespace.
    hard_errors = {
        k: v
        for k, v in problems.items()
        if k in ("ensembl_gene", "ensembl_transcript", "entrez", "whitespace")
    }
    if hard_errors:
        raise ValueError(_build_error_message(hard_errors))

    # Mouse symbols: only error if >50% of all unique symbols match.
    mouse_tokens = problems.get("mouse_symbol", [])
    if mouse_tokens and len(mouse_tokens) / len(all_genes) > 0.5:
        raise ValueError(
            f"Found possible mouse gene symbols (e.g., {mouse_tokens[0]}). "
            "This toolkit operates on human genes. Convert to human "
            "orthologs first (e.g., via Ensembl BioMart ortholog mapping)."
        )

    return ParseResult(
        gene_sets=gene_sets,
        n_lines_parsed=n_lines_parsed,
        n_lines_skipped=n_lines_skipped,
        n_duplicates=n_duplicates,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# HGNC validation
# ---------------------------------------------------------------------------


@dataclass
class ValidationResult:
    """Result of validating gene symbols against HGNC."""

    recognized: Set[str] = field(default_factory=set)
    aliases_resolved: Dict[str, str] = field(default_factory=dict)
    unrecognized: Set[str] = field(default_factory=set)
    gene_sets: Dict[str, List[str]] = field(default_factory=dict)


def validate_with_hgnc(
    gene_sets: Dict[str, List[str]],
    hgnc_path: str,
) -> ValidationResult:
    """Validate gene symbols against HGNC and resolve aliases.

    For each gene symbol:

    1. If it is a current HGNC-approved symbol, it is recognized.
    2. If it is a known alias or previous symbol, it is resolved to the
       canonical symbol and the gene set entry is updated.
    3. If it matches nothing, it is marked unrecognized (included as-is
       with a warning).

    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Gene sets from :func:`parse_geneset_tsv`.
    hgnc_path : str
        Path to HGNC Hail Table (``.ht``) or TSV file.

    Returns
    -------
    ValidationResult
    """
    # transient: Task 4 (#121) replaces this with a GeneCatalogStreamer parameter
    from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

    catalog = HGNCGeneCatalogStreamer.from_path(hgnc_path)
    canonical_symbols = catalog._canonical_symbols
    alias_to_canonical = catalog._alias_to_canonical

    # Classify every unique gene across all sets.
    all_genes = {g for genes in gene_sets.values() for g in genes}
    recognized: Set[str] = set()
    aliases_resolved: Dict[str, str] = {}
    unrecognized: Set[str] = set()

    for gene in all_genes:
        if gene in canonical_symbols:
            recognized.add(gene)
        elif gene in alias_to_canonical:
            canonical = alias_to_canonical[gene]
            aliases_resolved[gene] = canonical
            recognized.add(canonical)
        else:
            unrecognized.add(gene)

    # Rebuild gene sets with aliases resolved, deduplicating.
    resolved_sets: Dict[str, List[str]] = {}
    for set_name, genes in gene_sets.items():
        seen: Set[str] = set()
        resolved: List[str] = []
        for gene in genes:
            canonical = aliases_resolved.get(gene, gene)
            if canonical not in seen:
                seen.add(canonical)
                resolved.append(canonical)
            else:
                logger.warning(
                    "Alias resolution created duplicate '%s' (from '%s') "
                    "in set '%s', keeping first occurrence.",
                    canonical,
                    gene,
                    set_name,
                )
        resolved_sets[set_name] = resolved

    return ValidationResult(
        recognized=recognized,
        aliases_resolved=aliases_resolved,
        unrecognized=unrecognized,
        gene_sets=resolved_sets,
    )
