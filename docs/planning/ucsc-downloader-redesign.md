# UCSC Downloader Redesign

**Date:** 2026-03-11
**Status:** Phases 0-2 implemented; Phase 3 planned
**Modules:** `hvantk/datasets/ucsc_cell_datasets.py`, `hvantk/commands/ucsc_downloader.py`

## Problem Statement

The UCSC Cell Browser downloader cannot download data from **collection datasets**
(56% of the catalog — 149 out of 267 entries). Collections like `hoc` have no
expression matrix at the top level; their data lives in child datasets at
subdirectory paths (e.g., `hoc/all-heart`).

**Root causes:**
1. The local catalog JSON (`cells_ucsc_datasets.json`) only stores top-level
   entries — no children for collections.
2. The downloader constructs URLs as `base_url/{name}/{filename}`, but child
   dataset names include parent paths (`hoc/all-heart`), which are rejected by
   the slash validation.
3. There is no mechanism to distinguish downloadable leaf datasets from
   non-downloadable collections before attempting the download.

**UCSC API observations** (from `cells.ucsc.edu`):
- `cells.ucsc.edu/dataset.json` — root catalog listing all top-level entries.
- `cells.ucsc.edu/{name}/dataset.json` — per-dataset metadata; collections
  include a `datasets` array with children (full paths like `hoc/all-heart`);
  leaf datasets include `fileVersions` with actual file paths/sizes.
- Expression files follow the convention:
  `cells.ucsc.edu/{full_path}/exprMatrix.tsv.gz` and
  `cells.ucsc.edu/{full_path}/meta.tsv`.

---

## Design Principles

1. **Offline-first**: The local JSON catalog is the primary source. Network
   calls happen only for collection resolution and explicit refresh.
2. **Graceful degradation**: If UCSC API is unreachable, the tool still works
   with cached data (leaf datasets always downloadable, collections show
   "run --refresh to resolve children").
3. **Backward compatible**: Existing commands like
   `--dataset adultPancreas` continue to work unchanged.
4. **Minimal network**: Fetch only what's needed — one HTTP request per
   collection resolution, not a full catalog crawl.

---

## Phase 0 — Quick Unblock (Minimal Fix) ✅ Implemented

> Goal: Allow power users to download collection children immediately.

### 0.1 Relax slash validation

**File:** `hvantk/commands/ucsc_downloader.py`

Current validation rejects any dataset name containing `/`. Change to:
- Allow forward slashes in dataset names (they're part of UCSC's path scheme).
- Keep rejecting `..` (path traversal) and `\` (Windows paths).
- Keep rejecting whitespace.

```python
# Before
invalid = (
    (".." in dataset)
    or any(ch in dataset for ch in ["/", "\\"])
    or any(ch.isspace() for ch in dataset)
)

# After
invalid = (
    (".." in dataset)
    or ("\\" in dataset)
    or any(ch.isspace() for ch in dataset)
)
```

### 0.2 Warn on collection names

When the user provides a name that matches a known collection (has
`isCollection=True` in the catalog), print a warning before the URL check fails:

```
Warning: "hoc" is a collection with 10 child datasets.
Try a child dataset path instead (e.g., "hoc/all-heart").
Use --list_datasets --search hoc to find children.
```

This requires loading the catalog and checking `isCollection` early in the
command flow, before URL construction.

### 0.3 Fix target directory for slashed names

Current code: `target_dir = os.path.join(output_dir, dataset)`.
For `dataset="hoc/all-heart"`, this creates `data/ucsc/hoc/all-heart/` which is
fine — the nested structure mirrors UCSC's. No change needed, but add a test to
confirm this works.

### Tests

- Test that `--dataset hoc/all-heart` passes validation.
- Test that `--dataset ../etc/passwd` is rejected.
- Test that `--dataset hoc` (collection) prints a warning.
- Test target directory path for slashed dataset names.

---

## Phase 1 — Collection Resolution (partially implemented)

> Goal: When a user provides a collection name, resolve its children
> dynamically and let them pick one.

### 1.1 Add `fetch_children()` to `UCSCDataset` ✅ Implemented

**File:** `hvantk/datasets/ucsc_cell_datasets.py`

```python
def fetch_children(self) -> list["UCSCDataset"]:
    """Fetch child datasets from UCSC API for a collection.

    Returns empty list if not a collection or if fetch fails.
    """
    if not self.isCollection:
        return []

    import requests
    url = f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/dataset.json"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    data = resp.json()

    children = []
    for child in data.get("datasets", []):
        children.append(UCSCDataset(
            shortLabel=child.get("shortLabel", ""),
            name=child["name"],
            md5=child.get("md5", ""),
            sampleCount=child.get("sampleCount"),
            isCollection=child.get("isCollection", False),
            datasetCount=child.get("datasetCount"),
            body_parts=child.get("body_parts", []),
            organisms=child.get("organisms", []),
            diseases=child.get("diseases", []),
        ))
    return children
```

### 1.2 Add `children` cache field to `UCSCDataset` ✅ Implemented

Add optional field:

```python
children: Optional[List["UCSCDataset"]] = field(default=None, repr=False)
```

`fetch_children()` populates this. When serialized back to JSON
(Phase 3), children are included so subsequent runs don't need network.

### 1.3 Interactive collection resolution in downloader

**File:** `hvantk/commands/ucsc_downloader.py`

When `--dataset` matches a collection:

```
"hoc" is a collection with 10 child datasets.
Fetching children from UCSC...

  1. hoc/all-heart         (142,946 cells)
  2. hoc/blood             (12,345 cells)
  3. hoc/cardiomyocyte     (8,901 cells)
  ...

Select a dataset [1-10], or 0 to cancel:
```

Use `click.prompt()` with `type=click.IntRange(0, n)` for input.

If `--dataset hoc/all-heart` is given directly (already a leaf path), skip
resolution and download immediately.

### 1.4 Non-interactive mode

Add `--non-interactive` flag (or detect piped stdin). In non-interactive mode,
collection names print the child list and exit with code 1 instead of
prompting.

### Tests

- Mock `requests.get` for `hoc/dataset.json` → return fixture with children.
- Test interactive resolution (mock `click.prompt`).
- Test non-interactive mode exits with child list.
- Test nested collections (collection whose child is also a collection).
- Test network failure graceful fallback.

---

## Phase 2 — Search, Filter, and Rich Display (mostly implemented)

> Goal: Make `--list_datasets` useful for discovering datasets.

### 2.1 Add `--search` filter ✅ Implemented

**File:** `hvantk/commands/ucsc_downloader.py`

```python
@click.option(
    "--search",
    type=str,
    default=None,
    help="Filter dataset names and labels by search term (case-insensitive).",
)
```

Search matches against `name`, `shortLabel`, `body_parts`, `organisms`, and
`diseases`. Example:

```bash
$ hvantk ucsc-downloader --list_datasets --search heart
Available datasets (filtered by "heart"):
  heart-cell-atlas         (collection, 10 datasets)
  hoc                      (collection, 10 datasets) [body: heart]
  mouse-cardiac            (collection, 9 datasets)  [body: heart]
  mouse-dev-heart          (collection, 4 datasets)  [body: heart]
  cardiogenesis-atac       (collection, 2 datasets)  [body: heart]
```

### 2.2 Add `search()` method to `UCSCDataSetCollection` ✅ Implemented

**File:** `hvantk/datasets/ucsc_cell_datasets.py`

```python
def search(self, query: str) -> "UCSCDataSetCollection":
    """Filter datasets by case-insensitive query across name, label, and facets."""
    q = query.lower()
    matches = []
    for ds in self.datasets:
        searchable = " ".join([
            ds.name,
            ds.shortLabel,
            " ".join(ds.body_parts or []),
            " ".join(ds.organisms or []),
            " ".join(ds.diseases or []),
        ]).lower()
        if q in searchable:
            matches.append(ds)
    return UCSCDataSetCollection(datasets=matches)
```

### 2.3 Rich display format ✅ Implemented

Improve `_print_dataset_names()` to show:
- Collection vs leaf distinction (icon or label).
- Sample count for leaf datasets.
- Key facets (organism, body part) when they match the search.
- Children (indented) when `--search` matches a collection and children are
  cached.

```
Available datasets (267 total, filtered by "heart"):

  heart-cell-atlas  (collection, 10 datasets)  [Human, heart]
    ├── heart-cell-atlas/fetal-heart      (45,000 cells)
    ├── heart-cell-atlas/adult-heart      (32,000 cells)
    └── ... 8 more (use --dataset heart-cell-atlas to browse)

  hoc  (collection, 10 datasets)  [Human, heart]
    ├── hoc/all-heart                     (142,946 cells)
    ├── hoc/cardiomyocyte                 (8,901 cells)
    └── ... 8 more

  mouse-dev-heart  (collection, 4 datasets)  [Mouse, heart]

5 datasets found.
```

### 2.4 Add `--organism` filter

Optional convenience filter alongside `--search`:

```bash
$ hvantk ucsc-downloader --list_datasets --organism human --search pancreas
```

### Tests

- Test `search()` finds by name, shortLabel, body_parts.
- Test case-insensitive matching.
- Test empty results.
- Test `--search` CLI integration (mock catalog).
- Test display format for collections vs leaves.

---

## Phase 3 — Catalog Refresh and Caching

> Goal: Keep the local catalog up-to-date without manual JSON file management.

### 3.1 Add `--refresh` flag to downloader

```bash
$ hvantk ucsc-downloader --list_datasets --refresh
Fetching catalog from UCSC Cell Browser...
Updated: 283 datasets (was 267).
Resolving collections... 149 collections found.
Catalog saved to ~/.hvantk/ucsc_catalog.json
```

### 3.2 Catalog fetch and merge

**File:** `hvantk/datasets/ucsc_cell_datasets.py`

```python
@classmethod
def from_remote(cls, base_url: str = UCSC_CELL_BROWSER_BASE_URL) -> "UCSCDataSetCollection":
    """Fetch the latest catalog from UCSC Cell Browser API."""
    import requests
    resp = requests.get(f"{base_url}/dataset.json", timeout=30)
    resp.raise_for_status()
    data = resp.json()
    datasets = [UCSCDataset(**ds) for ds in data.get("datasets", [])]
    return cls(datasets=datasets)
```

### 3.3 User-local cache

Store the enriched catalog (with resolved children) at
`~/.hvantk/ucsc_catalog.json`. Loading priority:

1. User cache (`~/.hvantk/ucsc_catalog.json`) if it exists and is < 30 days old.
2. Bundled catalog (`hvantk/resources/cells_ucsc_datasets.json`) as fallback.
3. Remote fetch with `--refresh`.

Add `UCSC_CATALOG_CACHE_PATH` and `UCSC_CATALOG_MAX_AGE_DAYS` to constants.

### 3.4 Bulk children resolution

When `--refresh` is used, resolve children for all collections in parallel
(bounded concurrency with `concurrent.futures.ThreadPoolExecutor`):

```python
def resolve_all_children(self, max_workers: int = 5) -> None:
    """Fetch children for all collection datasets in parallel."""
```

This is the most network-intensive operation. Only runs on explicit
`--refresh`, never automatically.

### 3.5 Serialize enriched catalog

Add `to_json()` method to `UCSCDataSetCollection` that writes the full catalog
including children:

```json
{
  "datasets": [
    {
      "name": "hoc",
      "shortLabel": "Heart of Cells",
      "isCollection": true,
      "datasetCount": 10,
      "children": [
        {"name": "hoc/all-heart", "shortLabel": "All Heart", "sampleCount": 142946},
        {"name": "hoc/blood", "shortLabel": "Blood", "sampleCount": 12345}
      ]
    }
  ],
  "metadata": {
    "fetched_at": "2026-03-11T18:00:00Z",
    "source": "https://cells.ucsc.edu/dataset.json",
    "total_datasets": 283,
    "total_leaf_datasets": 450
  }
}
```

### Tests

- Test `from_remote()` with mocked HTTP.
- Test cache loading priority (user cache > bundled > remote).
- Test cache age check (> 30 days triggers warning).
- Test `to_json()` round-trip (write → read → same data).
- Test parallel children resolution with mock.

---

## Implementation Order

| Step | Phase | Effort | Dependencies | Impact | Status |
|------|-------|--------|--------------|--------|--------|
| 1    | 0.1   | Small  | None         | Unblocks power users immediately | ✅ Done |
| 2    | 0.2   | Small  | None         | Better error messages | ✅ Done |
| 3    | 2.1-2.2 | Small | None       | Search/filter for discovery | ✅ Done |
| 4    | 1.1-1.2 | Small | 0.1        | `fetch_children()` + cache field | ✅ Done |
| 5    | 2.3   | Small  | 2.1, 1.1    | Rich display with children | ✅ Done |
| 6    | 1.3   | Medium | 1.1         | Interactive collection resolution in downloader | Planned |
| 7    | 1.4   | Small  | 1.3         | Non-interactive/scripting support | Planned |
| 8    | 3.1-3.5 | Medium | 1.1       | Catalog refresh and caching | Planned |
| 9    | 2.4   | Small  | 2.1         | Organism filter (polish) | Planned |

**Remaining work:** Steps 6-9 (interactive resolution, non-interactive mode,
catalog caching, organism filter).

---

## Files Modified

| File | Phases | Changes |
|------|--------|---------|
| `hvantk/commands/ucsc_downloader.py` | 0, 1, 2 | Slash validation, collection warning, search, interactive resolution |
| `hvantk/datasets/ucsc_cell_datasets.py` | 1, 2, 3 | `fetch_children()`, `search()`, `from_remote()`, `to_json()`, `children` field |
| `hvantk/core/constants.py` | 3 | Cache path and max age constants |
| `hvantk/tests/test_ucsc_downloader.py` | 0, 1, 2 | New tests for validation, collection resolution, search |
| `hvantk/tests/test_ucsc_dataset_collection.py` | 1, 2, 3 | Tests for `fetch_children()`, `search()`, caching |

---

## Open Questions

1. **Nested collections**: Some children may themselves be collections (e.g.,
   `heart-cell-atlas` → `heart-cell-atlas/fetal` → further children). How deep
   should resolution go? Proposal: one level only in Phase 1; recursive in
   Phase 3 with `--refresh`.

2. **Bundled JSON update frequency**: Should we periodically regenerate the
   bundled `cells_ucsc_datasets.json` as part of releases? Or rely entirely on
   `--refresh` for up-to-date data?

3. **Dataset validation before download**: Phase 1 fetches `dataset.json` to
   list children. Should we also check `fileVersions` to confirm the leaf
   dataset actually has expression/metadata files before downloading?
