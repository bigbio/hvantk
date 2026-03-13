# Plan: `hvantk prepare-geneset` CLI Command

## Context

hvantk follows a CLI-first design. The current pipeline for running gene-set-based
analyses (PSROC, EnrichEx burden, EnrichEx overlap) requires a `GeneSetCollection`
JSON file. Every step in the workflow has a CLI command except the critical
preparation step:

```
download data     →  hvantk clinvar-downloader       ✅ CLI
build tables      →  hvantk mktable clinvar          ✅ CLI
prepare gene sets →  ???                             ❌ Python only
run analysis      →  hvantk psroc --gene-sets        ✅ CLI
```

A typical end-user scenario: a clinical researcher has 10 gene panels from their
lab as plain text files or a spreadsheet export. Today they must write a Python
script to convert these into the JSON format that `--gene-sets` expects. This
breaks the CLI-first workflow and misses the opportunity to validate gene symbols
before launching expensive Hail computations.

### Motivation

1. **Close the CLI workflow gap** — no Python detour required to go from text
   files to analysis.
2. **Validate before compute** — catch gene symbol mismatches, aliases, and
   Ensembl IDs *before* a 30-minute Spark job silently drops genes.
3. **Standardize output** — guaranteed `GeneSetCollection` JSON compatible with
   `hvantk psroc --gene-sets` and `hvantk enrichex burden -s`.

### Design Principles

- **One input format, strict** — reduces parser complexity and error surface.
- **HGNC gene symbols only** — this is what all downstream pipelines match on
  (`SYMBOL` field in burden, `GENEINFO` in PSROC, set intersection in overlap).
- **Human genes only** — all reference data (ClinVar, dbNSFP, gnomAD, HGNC) is
  human. Cross-species ortholog mapping is out of scope.
- **Clear errors for everything else** — detect and reject Ensembl IDs, mouse
  symbols, and mixed formats with actionable messages.

---

## Input Specification

### Format: Headerless Two-Column TSV

```
<gene_set_name>\t<gene_symbol>
```

One line per gene-set membership. A gene can appear in multiple gene sets.
Lines starting with `#` are comments. Blank lines are ignored.

**Example** (`panels.tsv`):

```
cardiac_panel	MYH7
cardiac_panel	TNNT2
cardiac_panel	LMNA
cardiac_panel	SCN5A
epilepsy_panel	SCN1A
epilepsy_panel	SCN2A
epilepsy_panel	KCNQ2
epilepsy_panel	STXBP1
ras_pathway	BRAF
ras_pathway	KRAS
ras_pathway	NRAS
```

**Why this format:**
- Two columns forces the user to name their panels explicitly.
- No header avoids "is row 1 data or metadata?" ambiguity.
- TSV is universally exportable from Excel, Google Sheets, R, Python.
- One gene per line (not comma-separated) avoids quoting issues and is
  `grep`/`sort`/`wc -l` friendly.

**Why only this format (not GMT, JSON, plain text directories):**
- GMT and JSON are already loadable by `--gene-sets` directly (PSROC and
  `load_gene_sets()` auto-detect by extension). This command targets users
  who do *not* have their data in those formats yet.
- Supporting multiple input formats multiplies parser code and error modes
  for marginal benefit. Users with GMT files don't need this command.

### Gene Identifier: HGNC Symbol

The command accepts **only HGNC-approved human gene symbols** (e.g., `BRCA1`,
`TP53`, `SCN1A`). This constraint exists because:

- Burden analysis defaults to `gene_field="SYMBOL"` and matches gene sets
  against the MatrixTable's VEP-annotated `SYMBOL` row field.
- PSROC filters ClinVar's `GENEINFO` field, which contains HGNC symbols.
- Overlap analysis uses pure set intersection; all existing examples use symbols.
- The HGNC alias resolution infrastructure (`gene_aliases.py`) operates on symbols.

### What We Explicitly Reject

| Pattern detected | Action | Message |
|---|---|---|
| Ensembl gene IDs (`ENSG00000...`) | **Error**, abort | `"Found Ensembl gene IDs (e.g., {example}). This command requires HGNC gene symbols. Convert first using biomart, hvantk GeneMapper, or https://www.genenames.org/tools/multi-symbol-checker/"` |
| Ensembl transcript IDs (`ENST00000...`) | **Error**, abort | Same as above, noting transcript vs gene distinction |
| Mouse-style symbols (lowercase: `Trp53`, `Scn1a`) | **Error**, abort | `"Found possible mouse gene symbols (e.g., {example}). This toolkit operates on human genes. Convert to human orthologs first (e.g., via Ensembl BioMart ortholog mapping)."` |
| Entrez numeric IDs (`7157`) | **Error**, abort | `"Found numeric-only entries (e.g., {example}) that may be Entrez Gene IDs. This command requires HGNC gene symbols."` |
| Whitespace in gene names | **Error**, abort | `"Gene symbol contains whitespace: '{token}'. Check input formatting."` |
| Empty gene symbol (blank second column) | **Warning**, skip line | `"Line {n}: empty gene symbol, skipping."` |
| Duplicate gene in same set | **Warning**, deduplicate | `"Duplicate gene '{gene}' in set '{set_name}', keeping first occurrence."` |

Detection heuristics:

```python
import re

def detect_id_type(token: str) -> str:
    """Classify a gene identifier token."""
    if re.match(r"^ENSG\d{11}(\.\d+)?$", token):
        return "ensembl_gene"
    if re.match(r"^ENST\d{11}(\.\d+)?$", token):
        return "ensembl_transcript"
    if re.match(r"^\d+$", token):
        return "entrez"
    if token[0].isupper() and any(c.islower() for c in token[1:3]) and token != token.upper():
        # Heuristic: mouse symbols are capitalized first letter only (Trp53)
        # Human symbols are all-caps or mixed with numbers (TP53, SCN1A, HLA-DRB1)
        # This catches the common case but is not perfect
        if re.match(r"^[A-Z][a-z]", token) and not re.search(r"[A-Z]{2}", token):
            return "mouse_symbol"
    if " " in token or "\t" in token:
        return "whitespace"
    return "symbol"
```

**Note on the mouse symbol heuristic**: This is intentionally conservative. Some
valid human symbols (e.g., `Oct4`) could false-positive. The detection should scan
all entries first and only raise the error if >50% match the mouse pattern, to
avoid false positives on edge cases. Individual ambiguous symbols are fine —
the HGNC validation step (with `--hgnc`) is the authoritative check.

---

## Output Specification

### Format: `GeneSetCollection` JSON

The output follows the existing `GeneSetCollection.save()` schema, ensuring
direct compatibility with:

- `hvantk psroc --gene-sets output.json`
- `hvantk enrichex burden -s output.json`
- `hvantk enrichex overlap -s output.json`
- `GeneSetCollection.load("output.json")` in Python

```json
{
  "gene_sets": {
    "cardiac_panel": {
      "name": "cardiac_panel",
      "genes": ["LMNA", "MYH7", "SCN5A", "TNNT2"],
      "source": "prepare-geneset",
      "n_genes": 4,
      "metadata": {}
    },
    "epilepsy_panel": {
      "name": "epilepsy_panel",
      "genes": ["KCNQ2", "SCN1A", "SCN2A", "STXBP1"],
      "source": "prepare-geneset",
      "n_genes": 4,
      "metadata": {}
    }
  },
  "background_genes": ["BRAF", "KCNQ2", "KRAS", "LMNA", "MYH7", "NRAS", "SCN1A", "SCN2A", "SCN5A", "STXBP1", "TNNT2"],
  "n_gene_sets": 3,
  "n_background": 11,
  "source_description": "Custom gene sets from panels.tsv",
  "metadata": {
    "created_by": "hvantk prepare-geneset",
    "input_file": "panels.tsv",
    "hgnc_validated": false
  }
}
```

When `--hgnc` is provided and aliases are resolved, `metadata` additionally
includes:

```json
{
  "metadata": {
    "created_by": "hvantk prepare-geneset",
    "input_file": "panels.tsv",
    "hgnc_validated": true,
    "hgnc_path": "/data/hgnc.ht",
    "aliases_resolved": {"FANCD1": "BRCA2", "ERCC11": "ERCC1"},
    "unrecognized_symbols": ["FAKEGENE", "NOTREAL"]
  }
}
```

### Background Genes

`background_genes` is set to the **union of all genes across all sets**. This is
the same default as `GeneSetCollection.load_gmt(background_strategy="union")`.

The user can supply an explicit background via `--background` (a one-column text
file of gene symbols). This is important for overlap analysis where the background
universe affects the Fisher's exact test denominator.

---

## CLI Interface

### Command Signature

```
hvantk prepare-geneset [OPTIONS]
```

### Options

```
Required:
  -i, --input PATH          Input TSV file (headerless, two columns:
                            gene_set_name<TAB>gene_symbol). One gene per line,
                            # comments and blank lines ignored.

  -o, --output PATH         Output path for GeneSetCollection JSON file.

Optional — Validation:
  --hgnc PATH               Path to HGNC data (Hail Table .ht or TSV) for
                            symbol validation and alias resolution. When
                            provided, each gene symbol is checked against
                            the HGNC canonical set. Aliases and previous
                            symbols are resolved to current approved symbols.
                            Unrecognized symbols are reported as warnings.

Optional — Filtering:
  --min-genes INTEGER       Exclude gene sets with fewer than this many genes
                            after validation. [default: 0, no filtering]

  --background PATH         Text file with background gene universe (one gene
                            symbol per line). If not provided, background is
                            the union of all genes across all sets.

Optional — Output:
  --overwrite               Overwrite output file if it exists.

  --export-gmt PATH         Additionally export as GMT file (for GSEA/MSigDB
                            tool compatibility).

General:
  --log-level               DEBUG|INFO|WARNING|ERROR [default: INFO]
  --help                    Show this message and exit.
```

### Usage Examples

```bash
# Basic: convert TSV panels to JSON
hvantk prepare-geneset \
  -i lab_panels.tsv \
  -o gene_sets.json

# With HGNC validation and alias resolution
hvantk prepare-geneset \
  -i lab_panels.tsv \
  -o gene_sets.json \
  --hgnc /data/tables/hgnc.ht

# With filtering and explicit background
hvantk prepare-geneset \
  -i lab_panels.tsv \
  -o gene_sets.json \
  --hgnc /data/tables/hgnc.ht \
  --min-genes 5 \
  --background protein_coding_genes.txt

# Full pipeline example (end-to-end)
hvantk prepare-geneset -i panels.tsv -o panels.json --hgnc /data/hgnc.ht
hvantk psroc --gene-sets panels.json --clinvar-ht /data/clinvar.ht \
  --dbnsfp-ht /data/dbnsfp.ht --scores "CADD_phred,REVEL_score" -o /results/
```

### Console Output

The command prints a validation summary to stderr:

```
Loading gene sets from lab_panels.tsv ...
Parsed 342 genes across 10 gene sets.

Validating against HGNC (/data/tables/hgnc.ht) ...
  cardiac_panel:    45 genes → 43 recognized, 2 aliases resolved
                    FANCD1 → BRCA2, ERCC11 → ERCC1
  epilepsy_panel:   38 genes → 37 recognized, 0 aliases, 1 unrecognized
                    ⚠ Unrecognized: FAKEGENE
  neuro_panel:      52 genes → 52 recognized ✓
  ...

Summary:
  Gene sets:           10
  Total unique genes:  312 (after alias resolution)
  Recognized:          310 (99.4%)
  Aliases resolved:    4
  Unrecognized:        2 (FAKEGENE, NOTREAL)

⚠ 2 unrecognized symbols will be included as-is. They will not match
  any gene in downstream analyses. Review or remove them from your input.

Saved to gene_sets.json
```

---

## Implementation Steps

### Step 1: Add Input Parser

**New file**: `hvantk/utils/geneset_io.py`

This module contains the pure parsing and validation logic, independent of Click.
Keeping it separate from the CLI command allows reuse from the Python API.

```python
"""
Gene set input parsing and validation for prepare-geneset.

Parses a headerless two-column TSV (gene_set_name<TAB>gene_symbol) into a
dictionary of {set_name: List[str]}, with format and identifier validation.
"""

from pathlib import Path
from typing import Dict, List, Set, Tuple
import logging
import re

logger = logging.getLogger(__name__)


# --- ID type detection ---

_ENSEMBL_GENE_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$")
_ENSEMBL_TRANSCRIPT_RE = re.compile(r"^ENST\d{11}(\.\d+)?$")
_ENTREZ_RE = re.compile(r"^\d+$")
_MOUSE_SYMBOL_RE = re.compile(r"^[A-Z][a-z][a-z0-9]+$")


def detect_id_type(token: str) -> str:
    """Classify a gene identifier token.

    Returns one of: "ensembl_gene", "ensembl_transcript", "entrez",
    "mouse_symbol", "whitespace", "symbol".
    """
    ...


def validate_gene_ids(
    genes: List[str],
) -> Tuple[List[str], Dict[str, List[str]]]:
    """Check a list of gene tokens for non-symbol identifiers.

    Returns
    -------
    valid : List[str]
        Tokens classified as "symbol".
    problems : Dict[str, List[str]]
        Mapping of problem type → list of offending tokens.
        Keys: "ensembl_gene", "ensembl_transcript", "entrez",
        "mouse_symbol", "whitespace".
    """
    ...


# --- TSV parsing ---

class ParseResult:
    """Result of parsing a gene set TSV file."""
    gene_sets: Dict[str, List[str]]   # set_name → [genes] (deduplicated, ordered)
    n_lines_parsed: int
    n_lines_skipped: int               # comments + blanks
    n_duplicates: int                   # duplicate gene-in-set entries
    warnings: List[str]                 # human-readable warning messages


def parse_geneset_tsv(path: Path) -> ParseResult:
    """Parse a headerless two-column TSV into gene sets.

    Parameters
    ----------
    path : Path
        Input file. Expected format: gene_set_name<TAB>gene_symbol,
        one line per membership. # comments and blank lines ignored.

    Returns
    -------
    ParseResult

    Raises
    ------
    FileNotFoundError
        If path does not exist.
    ValueError
        If file has wrong number of columns (not exactly 2) on any
        non-comment, non-blank line.
    """
    ...
```

**Key behaviors:**

- Parse line by line. Strip whitespace. Skip `#` comments and blank lines.
- Split on `\t`. If a line does not produce exactly 2 fields, raise
  `ValueError` with the line number and content.
- Collect `{set_name: [gene_symbol, ...]}` preserving insertion order.
- Track duplicates (same gene in same set) → warn and deduplicate.
- After parsing all lines, call `validate_gene_ids()` on the full list of
  unique gene symbols. If any `problems` dict is non-empty, raise
  `ValueError` with the categorized error message (see rejection table above).
  Use the >50% threshold for mouse symbol detection to avoid false positives.

**Tests**: `hvantk/tests/test_geneset_io.py`

- Valid two-column file → correct gene_sets dict.
- Comments and blanks → skipped.
- Wrong column count → `ValueError` with line number.
- Ensembl IDs → `ValueError` naming the offending tokens.
- Mouse symbols → `ValueError` when majority are mouse-pattern.
- Entrez IDs → `ValueError`.
- Duplicates → warning in ParseResult, deduplicated output.
- Empty gene symbol → warning, line skipped.
- Mixed valid + few mouse-like → passes (below threshold).

---

### Step 2: Add HGNC Validation Function

**File**: `hvantk/utils/geneset_io.py` (same module)

```python
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ValidationResult:
    """Result of validating gene symbols against HGNC."""
    recognized: Set[str]              # symbols found in HGNC canonical set
    aliases_resolved: Dict[str, str]  # original_symbol → canonical_symbol
    unrecognized: Set[str]            # symbols not found anywhere in HGNC
    gene_sets: Dict[str, List[str]]   # updated gene sets with aliases resolved


def validate_with_hgnc(
    gene_sets: Dict[str, List[str]],
    hgnc_path: str,
) -> ValidationResult:
    """Validate gene symbols against HGNC and resolve aliases.

    For each gene symbol:
    1. If it's a current HGNC-approved symbol → recognized.
    2. If it's a known alias or previous symbol → resolve to canonical,
       update the gene set entry, log the mapping.
    3. If it matches nothing → unrecognized (included as-is with warning).

    Uses the existing ``_load_hgnc_symbol_maps()`` from
    ``hvantk.utils.gene_aliases`` to load the HGNC canonical set,
    alias-to-canonical map, and canonical-to-aliases map. This avoids
    duplicating HGNC loading logic.

    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Gene sets from parse_geneset_tsv().
    hgnc_path : str
        Path to HGNC Hail Table (.ht) or TSV file.

    Returns
    -------
    ValidationResult
    """
    ...
```

**Key behaviors:**

- Reuse `_load_hgnc_symbol_maps()` from `hvantk/utils/gene_aliases.py` —
  this already supports both `.ht` and `.tsv` formats and builds the
  `canonical_symbols`, `alias_to_canonical`, and `canonical_to_aliases` maps.
- Iterate over all unique genes across all sets.
- For aliases: replace the gene in every set where it appears.
- Handle the edge case where resolving an alias creates a duplicate within a
  set (e.g., user listed both `FANCD1` and `BRCA2` → after resolution both
  become `BRCA2`). Deduplicate and warn.
- The function does NOT remove unrecognized symbols. It includes them as-is
  and reports them so the user can decide.

**Tests**: `hvantk/tests/test_geneset_io.py` (extend)

- All canonical symbols → `unrecognized` is empty.
- Known alias → resolved, appears in `aliases_resolved` map.
- Unknown symbol → appears in `unrecognized`.
- Alias that duplicates existing canonical in same set → deduplicated with warning.
- Both `.ht` and `.tsv` HGNC paths work (mock or use test fixtures).

---

### Step 3: Add CLI Command

**New file**: `hvantk/commands/prepare_geneset_cli.py`

```python
"""
CLI command for preparing gene set collections from plain TSV files.

Converts a headerless two-column TSV (gene_set_name<TAB>gene_symbol) into
a GeneSetCollection JSON file, with optional HGNC symbol validation.
"""

import logging
import click

logger = logging.getLogger(__name__)


@click.command(name="prepare-geneset")
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(exists=True),
    required=True,
    help="Input TSV file (headerless, two columns: "
    "gene_set_name<TAB>gene_symbol).",
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    required=True,
    help="Output path for GeneSetCollection JSON.",
)
@click.option(
    "--hgnc",
    type=click.Path(exists=True),
    default=None,
    help="Path to HGNC data (.ht or .tsv) for symbol validation "
    "and alias resolution.",
)
@click.option(
    "--min-genes",
    type=click.IntRange(min=0),
    default=0,
    help="Exclude gene sets with fewer than this many genes. "
    "[default: 0, no filtering]",
)
@click.option(
    "--background",
    type=click.Path(exists=True),
    default=None,
    help="Text file with background gene universe (one symbol per "
    "line). Default: union of all genes across sets.",
)
@click.option(
    "--export-gmt",
    type=click.Path(),
    default=None,
    help="Additionally export as GMT file.",
)
@click.option("--overwrite", is_flag=True, help="Overwrite output if exists.")
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    default="INFO",
)
def prepare_geneset_cmd(
    input_path,
    output,
    hgnc,
    min_genes,
    background,
    export_gmt,
    overwrite,
    log_level,
):
    """Prepare a GeneSetCollection JSON from a plain TSV file.

    Reads a headerless two-column TSV (gene_set_name<TAB>gene_symbol) and
    produces a GeneSetCollection JSON file compatible with:

    \b
      hvantk psroc --gene-sets OUTPUT
      hvantk enrichex burden -s OUTPUT
      hvantk enrichex overlap -s OUTPUT

    Optionally validates symbols against HGNC, resolving aliases and
    reporting unrecognized entries.

    \b
    Example:
      hvantk prepare-geneset -i panels.tsv -o panels.json --hgnc /data/hgnc.ht
    """
    ...
```

**Implementation flow inside the command function:**

```
1. parse_geneset_tsv(input_path)
   → abort on ValueError (wrong columns, bad ID types)

2. if --hgnc:
       validate_with_hgnc(gene_sets, hgnc_path)
       → print per-set validation summary
       → replace gene_sets with resolved version

3. if --min-genes > 0:
       filter out small sets, log which were removed

4. if --background:
       load_gene_set(path=background)
   else:
       union of all genes

5. Build GeneSetCollection via load_gene_sets_from_dict()
   → set source_description and metadata

6. collection.save(output)
   if --export-gmt: collection.save_gmt(export_gmt)

7. Print summary to stderr
```

**Tests**: `hvantk/tests/test_prepare_geneset_cli.py`

- Use `click.testing.CliRunner` to invoke the command.
- Valid TSV → JSON output matches expected structure.
- With `--hgnc` → aliases resolved in output JSON.
- With `--min-genes` → small sets excluded.
- With `--background` → background_genes in output matches file.
- Missing input → error.
- Bad format (3 columns) → error with line number.
- Ensembl IDs in input → error with guidance message.
- `--overwrite` flag behavior (refuse without, allow with).
- `--export-gmt` → GMT file written alongside JSON.

---

### Step 4: Wire into Main CLI

**File**: `hvantk/hvantk.py`

Add the import and `cli.add_command()` call:

```python
from hvantk.commands.prepare_geneset_cli import prepare_geneset_cmd

# Add after clingen_genesets_cmd (line 87), logically grouped with
# gene set preparation commands:
cli.add_command(prepare_geneset_cmd)  # Custom gene set preparation
```

---

### Step 5: Add Test Data Fixtures

**New files** in `hvantk/tests/testdata/prepare_geneset/`:

```
valid_panels.tsv          — 3 panels, 15 genes, all valid HGNC symbols
with_aliases.tsv          — includes FANCD1 (alias of BRCA2), ERCC11
with_ensembl_ids.tsv      — includes ENSG00000141510 (should fail)
with_mouse_symbols.tsv    — includes Trp53, Scn1a (should fail)
with_entrez_ids.tsv       — includes 7157, 672 (should fail)
with_duplicates.tsv       — same gene in same set twice
with_comments.tsv         — # comment lines and blank lines
wrong_columns.tsv         — 3-column line mixed in
background_genes.txt      — one-column gene list for --background
```

These are small, static text files (no Hail dependency). The HGNC validation
tests can use a minimal mock or the existing HGNC test fixtures if available.

---

### Step 6: Update Documentation

#### `docs/tools/psroc.md`

In the "Preparing Gene Set Collections" section, add a new subsection
**"From Custom Gene Panels (Lab / Spreadsheet)"** before the existing
"From ClinGen" subsection:

```markdown
### From Custom Gene Panels

If you have gene panels as text files or spreadsheet exports, use
`hvantk prepare-geneset` to convert them:

1. Format your panels as a headerless two-column TSV
   (gene_set_name<TAB>gene_symbol), one gene per line:

   ```
   cardiac	MYH7
   cardiac	TNNT2
   cardiac	LMNA
   epilepsy	SCN1A
   epilepsy	SCN2A
   ```

2. Convert to gene set JSON:

   ```bash
   # Basic conversion
   hvantk prepare-geneset -i panels.tsv -o panels.json

   # With HGNC validation (recommended)
   hvantk prepare-geneset -i panels.tsv -o panels.json --hgnc /data/hgnc.ht
   ```

3. Run PSROC:

   ```bash
   hvantk psroc --gene-sets panels.json \
     --clinvar-ht /data/clinvar.ht \
     --dbnsfp-ht /data/dbnsfp.ht \
     --scores "CADD_phred,REVEL_score" \
     -o /results/psroc
   ```

**Gene identifier requirements:**
- Use HGNC-approved human gene symbols (e.g., BRCA1, TP53, SCN1A).
- Ensembl IDs, mouse symbols, and Entrez IDs are not accepted.
- Use `--hgnc` to automatically resolve outdated aliases (e.g., FANCD1 → BRCA2).
```

#### `docs/tools/enrichex.md`

Add the same subsection in the "Gene Set Format" section, replacing the
generic "Gene Identifiers" note with a concrete reference to `prepare-geneset`.

#### `docs/library/usage.md`

Add a new section **"6) Prepare Custom Gene Sets"** between the existing
sections, showing the TSV → JSON → analysis flow.

---

## File Summary

| Action | File | Description |
|---|---|---|
| **New** | `hvantk/utils/geneset_io.py` | Parsing, ID detection, HGNC validation |
| **New** | `hvantk/commands/prepare_geneset_cli.py` | Click CLI command |
| **New** | `hvantk/tests/test_geneset_io.py` | Unit tests for parser + validation |
| **New** | `hvantk/tests/test_prepare_geneset_cli.py` | CLI integration tests |
| **New** | `hvantk/tests/testdata/prepare_geneset/*.tsv` | Test fixtures (7-8 small files) |
| **Edit** | `hvantk/hvantk.py` | Wire `prepare_geneset_cmd` into main CLI |
| **Edit** | `docs/tools/psroc.md` | Add "From Custom Gene Panels" subsection |
| **Edit** | `docs/tools/enrichex.md` | Add "From Custom Gene Panels" subsection |
| **Edit** | `docs/library/usage.md` | Add "Prepare Custom Gene Sets" section |

## Dependencies

No new dependencies. All building blocks exist:

- `gene_aliases._load_hgnc_symbol_maps()` — HGNC loading and alias resolution
- `gene_sets.GeneSetCollection` / `load_gene_sets_from_dict()` — output format
- `gene_sets.load_gene_set()` — background file loading
- `click` — CLI framework (already used everywhere)

## Out of Scope

- **Ensembl ID auto-conversion**: Detecting Ensembl IDs is in scope (as an
  error). Converting them is not — it would require HGNC/Ensembl mapping
  which adds complexity and a mandatory `--hgnc` dependency. Users should
  convert beforehand.
- **Mouse ortholog mapping**: Out of scope. Human-only toolkit.
- **GMT/JSON input**: These formats are already consumable by `--gene-sets`
  directly. This command targets users who don't have their data in those
  formats.
- **Interactive mode / web UI**: CLI-only.
- **Modifying existing `gene_sets.py`**: The new module `geneset_io.py`
  complements it. No changes to existing data structures needed.

## Implementation Order

1. **Step 1** (parser) + **Step 5** (test data) — can be done together, no
   Hail dependency, pure Python. Write tests first (TDD).
2. **Step 2** (HGNC validation) — depends on Step 1 output. Uses existing
   `gene_aliases` infrastructure.
3. **Step 3** (CLI) + **Step 4** (wire) — integrates Steps 1-2 with Click.
4. **Step 6** (docs) — after the command works end-to-end.

Estimated touch points: 4 new files, 4 edited files. No changes to existing
business logic. The feature is purely additive.
