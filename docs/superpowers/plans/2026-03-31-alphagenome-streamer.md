# AlphaGenome Variant Prediction Streamer — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an AlphaGenome API-driven streamer that produces per-modality Hail Tables from variant predictions.

**Architecture:** `AlphaGenomeStreamer` extends `HailDataStreamer`, reads variants from a Hail Table or TSV, groups them into genomic intervals (adaptive or fixed), calls the AlphaGenome `predict_variant()` API with rate limiting and checkpoint-based resumption, and assembles per-modality Hail Tables as output. A YAML config file controls API key, ontology terms, output types, and interval strategy.

**Tech Stack:** Python, Hail, `alphagenome` SDK (optional dep), PyYAML, Click CLI

**Spec:** `docs/superpowers/specs/2026-03-31-alphagenome-streamer-design.md`

---

## File Structure

| File | Responsibility |
|------|---------------|
| `hvantk/data/alphagenome_streamer.py` (new) | Config loading, `compute_intervals()`, `AlphaGenomeStreamer` class, checkpoint logic, rate-limited API calls |
| `hvantk/tables/table_builders.py` (modify) | Add `create_alphagenome_tb()` builder function |
| `hvantk/tables/registry.py` (modify) | Register `"alphagenome"` in `TABLE_BUILDERS` |
| `hvantk/commands/make_table_cli.py` (modify) | Add `mktable alphagenome` CLI command |
| `hvantk/core/constants.py` (modify) | Add AlphaGenome default constants |
| `hvantk/data/__init__.py` (modify) | Export `AlphaGenomeStreamer` |
| `hvantk/tests/test_alphagenome_streamer.py` (new) | Unit tests with mocked API |
| `hvantk/tests/testdata/alphagenome_config.yaml` (new) | Test config fixture |

---

### Task 1: Config Loading and Validation

**Files:**
- Create: `hvantk/tests/testdata/alphagenome_config.yaml`
- Create: `hvantk/tests/test_alphagenome_streamer.py`
- Create: `hvantk/data/alphagenome_streamer.py`

- [ ] **Step 1: Create test config fixture**

```yaml
# hvantk/tests/testdata/alphagenome_config.yaml
api:
  key: "test-api-key-123"
  max_retries: 3
  retry_backoff: 2.0
  request_timeout: 120

ontology:
  terms:
    - "UBERON:0001157"
    - "UBERON:0000955"
  output_types:
    - RNA_SEQ
    - CHROMATIN

intervals:
  default_size: 1048576
  adaptive: true
  adaptive_max_size: 1048576
  density_window: 50000
```

- [ ] **Step 2: Write failing tests for config loading**

```python
# hvantk/tests/test_alphagenome_streamer.py
import os
import pytest
import yaml

TESTDATA_DIR = os.path.join(os.path.dirname(__file__), "testdata")
CONFIG_PATH = os.path.join(TESTDATA_DIR, "alphagenome_config.yaml")


class TestLoadConfig:
    def test_load_valid_config(self):
        from hvantk.data.alphagenome_streamer import load_config

        config = load_config(CONFIG_PATH)
        assert config["api"]["key"] == "test-api-key-123"
        assert config["api"]["max_retries"] == 3
        assert config["ontology"]["terms"] == ["UBERON:0001157", "UBERON:0000955"]
        assert config["ontology"]["output_types"] == ["RNA_SEQ", "CHROMATIN"]
        assert config["intervals"]["adaptive"] is True

    def test_load_config_file_not_found(self):
        from hvantk.data.alphagenome_streamer import load_config

        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path.yaml")

    def test_load_config_missing_api_section(self, tmp_path):
        from hvantk.data.alphagenome_streamer import load_config

        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text(yaml.dump({"ontology": {"terms": []}}))
        with pytest.raises(ValueError, match="api"):
            load_config(str(bad_config))

    def test_load_config_missing_ontology_section(self, tmp_path):
        from hvantk.data.alphagenome_streamer import load_config

        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text(yaml.dump({"api": {"key": "x"}, "intervals": {}}))
        with pytest.raises(ValueError, match="ontology"):
            load_config(str(bad_config))

    def test_resolve_api_key_from_env(self, tmp_path, monkeypatch):
        from hvantk.data.alphagenome_streamer import load_config

        cfg = {
            "api": {"key": None, "max_retries": 3, "retry_backoff": 2.0, "request_timeout": 120},
            "ontology": {"terms": ["UBERON:0001157"], "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 1048576, "adaptive": True,
                          "adaptive_max_size": 1048576, "density_window": 50000},
        }
        cfg_path = tmp_path / "cfg.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        monkeypatch.setenv("ALPHAGENOME_API_KEY", "env-key-456")
        config = load_config(str(cfg_path))
        assert config["api"]["key"] == "env-key-456"

    def test_resolve_api_key_missing_everywhere(self, tmp_path, monkeypatch):
        from hvantk.data.alphagenome_streamer import load_config

        cfg = {
            "api": {"key": None, "max_retries": 3, "retry_backoff": 2.0, "request_timeout": 120},
            "ontology": {"terms": ["UBERON:0001157"], "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 1048576, "adaptive": True,
                          "adaptive_max_size": 1048576, "density_window": 50000},
        }
        cfg_path = tmp_path / "cfg.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
        with pytest.raises(ValueError, match="API key"):
            load_config(str(cfg_path))
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestLoadConfig -v`
Expected: FAIL with `ModuleNotFoundError` or `ImportError` (module doesn't exist yet)

- [ ] **Step 4: Implement config loading**

```python
# hvantk/data/alphagenome_streamer.py
"""AlphaGenome variant prediction streamer.

Streams variant positions through the AlphaGenome API and produces
per-modality Hail Tables with full multimodal predictions.
"""

import json
import logging
import os
import random
import time
from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional, Tuple

import yaml

from hvantk.data.data_streamer import HailDataStreamer

logger = logging.getLogger(__name__)

_REQUIRED_SECTIONS = ("api", "ontology")
_DEFAULT_INTERVALS = {
    "default_size": 1_048_576,
    "adaptive": True,
    "adaptive_max_size": 1_048_576,
    "density_window": 50_000,
}


def load_config(config_path: str) -> Dict[str, Any]:
    """Load and validate an AlphaGenome YAML config file.

    Auth resolution order: config api.key > ALPHAGENOME_API_KEY env var > error.

    Parameters
    ----------
    config_path : str
        Path to YAML config file.

    Returns
    -------
    dict
        Validated config dict with resolved API key and interval defaults.

    Raises
    ------
    FileNotFoundError
        If config_path does not exist.
    ValueError
        If required sections are missing or API key cannot be resolved.
    """
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    for section in _REQUIRED_SECTIONS:
        if section not in config or config[section] is None:
            raise ValueError(
                f"Config missing required section: '{section}'. "
                f"Required sections: {_REQUIRED_SECTIONS}"
            )

    # Resolve API key: config > env var > error
    api_key = config["api"].get("key")
    if not api_key:
        api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise ValueError(
            "API key not found. Set 'api.key' in config or "
            "ALPHAGENOME_API_KEY environment variable."
        )
    config["api"]["key"] = api_key

    # Apply interval defaults
    if "intervals" not in config or config["intervals"] is None:
        config["intervals"] = dict(_DEFAULT_INTERVALS)
    else:
        for k, v in _DEFAULT_INTERVALS.items():
            config["intervals"].setdefault(k, v)

    return config
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestLoadConfig -v`
Expected: All 6 tests PASS

- [ ] **Step 6: Commit**

```bash
git add hvantk/data/alphagenome_streamer.py hvantk/tests/test_alphagenome_streamer.py hvantk/tests/testdata/alphagenome_config.yaml
git commit -m "feat(alphagenome): add config loading and validation"
```

---

### Task 2: Interval Computation — `compute_intervals()`

**Files:**
- Modify: `hvantk/data/alphagenome_streamer.py`
- Modify: `hvantk/tests/test_alphagenome_streamer.py`

- [ ] **Step 1: Write failing tests for interval computation**

Append to `hvantk/tests/test_alphagenome_streamer.py`:

```python
@dataclass
class SimpleVariant:
    chrom: str
    pos: int
    ref: str
    alt: str


class TestComputeIntervals:
    """Test compute_intervals() as a pure function."""

    def _make_config(self, adaptive=True, default_size=1_048_576,
                     adaptive_max_size=1_048_576, density_window=50_000):
        return {
            "intervals": {
                "adaptive": adaptive,
                "default_size": default_size,
                "adaptive_max_size": adaptive_max_size,
                "density_window": density_window,
            }
        }

    def test_single_variant_fixed_mode(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [SimpleVariant("chr1", 500_000, "A", "T")]
        config = self._make_config(adaptive=False, default_size=100_000)
        result = compute_intervals(variants, config)
        assert len(result) == 1
        interval, grouped_variants = result[0]
        assert interval.chrom == "chr1"
        assert interval.start <= 500_000
        assert interval.end >= 500_000
        assert interval.end - interval.start == 100_000
        assert len(grouped_variants) == 1

    def test_two_close_variants_adaptive_grouped(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr1", 120_000, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        # Should be grouped into 1 interval
        assert len(result) == 1
        assert len(result[0][1]) == 2

    def test_two_distant_variants_adaptive_separate(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr1", 2_000_000, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        assert len(result) == 2
        assert len(result[0][1]) == 1
        assert len(result[1][1]) == 1

    def test_different_chromosomes_always_separate(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr2", 100_010, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        assert len(result) == 2

    def test_large_cluster_split_at_max_size(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        # Variants spanning > adaptive_max_size should be split
        variants = [SimpleVariant("chr1", i * 10_000, "A", "T") for i in range(200)]
        config = self._make_config(
            adaptive=True, density_window=50_000, adaptive_max_size=500_000
        )
        result = compute_intervals(variants, config)
        # With variants spanning 0..1_990_000 and max_size=500_000, need multiple intervals
        assert len(result) > 1
        for interval, group in result:
            assert interval.end - interval.start <= 500_000

    def test_interval_centered_on_variant(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [SimpleVariant("chr5", 1_000_000, "A", "T")]
        config = self._make_config(adaptive=False, default_size=200_000)
        result = compute_intervals(variants, config)
        interval = result[0][0]
        # Should be centered: midpoint of interval ~= variant position
        midpoint = (interval.start + interval.end) // 2
        assert abs(midpoint - 1_000_000) <= 1  # allow rounding

    def test_empty_variant_list(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        result = compute_intervals([], self._make_config())
        assert result == []
```

Add import at top of test file:

```python
from dataclasses import dataclass
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestComputeIntervals -v`
Expected: FAIL with `ImportError` (compute_intervals not defined)

- [ ] **Step 3: Implement data classes and `compute_intervals()`**

Add to `hvantk/data/alphagenome_streamer.py` after the `load_config` function:

```python
@dataclass
class VariantRecord:
    """A single variant position for AlphaGenome prediction."""
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass
class GenomicInterval:
    """A genomic interval for an AlphaGenome API call."""
    chrom: str
    start: int
    end: int


def compute_intervals(
    variants: List[Any],
    config: Dict[str, Any],
) -> List[Tuple[GenomicInterval, List[Any]]]:
    """Group variants into genomic intervals for API calls.

    In adaptive mode, nearby variants on the same chromosome are grouped
    into shared intervals to minimize API calls. In fixed mode, each
    variant gets its own interval.

    Parameters
    ----------
    variants : list
        Variant objects with .chrom, .pos, .ref, .alt attributes.
    config : dict
        Config dict with 'intervals' section.

    Returns
    -------
    list of (GenomicInterval, list of variants)
        Each tuple is an interval and the variants it contains.
    """
    if not variants:
        return []

    interval_cfg = config["intervals"]
    default_size = interval_cfg["default_size"]
    adaptive = interval_cfg.get("adaptive", True)

    if not adaptive:
        return _compute_fixed_intervals(variants, default_size)

    adaptive_max_size = interval_cfg.get("adaptive_max_size", default_size)
    density_window = interval_cfg.get("density_window", 50_000)
    return _compute_adaptive_intervals(
        variants, default_size, adaptive_max_size, density_window
    )


def _center_interval(pos: int, size: int) -> Tuple[int, int]:
    """Return (start, end) centered on pos with given size."""
    half = size // 2
    start = max(0, pos - half)
    end = start + size
    return start, end


def _compute_fixed_intervals(
    variants: List[Any], default_size: int
) -> List[Tuple[GenomicInterval, List[Any]]]:
    """One interval per variant, centered on variant position."""
    result = []
    for v in variants:
        start, end = _center_interval(v.pos, default_size)
        interval = GenomicInterval(chrom=v.chrom, start=start, end=end)
        result.append((interval, [v]))
    return result


def _compute_adaptive_intervals(
    variants: List[Any],
    default_size: int,
    max_size: int,
    density_window: int,
) -> List[Tuple[GenomicInterval, List[Any]]]:
    """Group nearby variants into shared intervals."""
    # Sort by (chrom, pos)
    sorted_variants = sorted(variants, key=lambda v: (v.chrom, v.pos))

    # Group consecutive variants on the same chrom within density_window
    groups: List[List[Any]] = []
    current_group: List[Any] = [sorted_variants[0]]

    for v in sorted_variants[1:]:
        prev = current_group[-1]
        if v.chrom == prev.chrom and (v.pos - prev.pos) <= density_window:
            current_group.append(v)
        else:
            groups.append(current_group)
            current_group = [v]
    groups.append(current_group)

    # Convert groups to intervals, splitting if they exceed max_size
    result = []
    for group in groups:
        span_start = group[0].pos
        span_end = group[-1].pos
        span = span_end - span_start

        if span <= max_size:
            # Fits in one interval — center on group midpoint
            midpoint = (span_start + span_end) // 2
            size = max(default_size, span + density_window)
            size = min(size, max_size)
            start, end = _center_interval(midpoint, size)
            interval = GenomicInterval(chrom=group[0].chrom, start=start, end=end)
            result.append((interval, group))
        else:
            # Split into sub-groups that fit within max_size
            sub_group: List[Any] = [group[0]]
            for v in group[1:]:
                if (v.pos - sub_group[0].pos) <= max_size:
                    sub_group.append(v)
                else:
                    _emit_subgroup(sub_group, default_size, max_size, density_window, result)
                    sub_group = [v]
            _emit_subgroup(sub_group, default_size, max_size, density_window, result)

    return result


def _emit_subgroup(
    sub_group: List[Any],
    default_size: int,
    max_size: int,
    density_window: int,
    result: List[Tuple[GenomicInterval, List[Any]]],
) -> None:
    """Create an interval for a sub-group and append to result."""
    sg_start = sub_group[0].pos
    sg_end = sub_group[-1].pos
    midpoint = (sg_start + sg_end) // 2
    size = max(default_size, (sg_end - sg_start) + density_window)
    size = min(size, max_size)
    start, end = _center_interval(midpoint, size)
    interval = GenomicInterval(chrom=sub_group[0].chrom, start=start, end=end)
    result.append((interval, sub_group))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestComputeIntervals -v`
Expected: All 7 tests PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/data/alphagenome_streamer.py hvantk/tests/test_alphagenome_streamer.py
git commit -m "feat(alphagenome): add compute_intervals with adaptive grouping"
```

---

### Task 3: Checkpoint Manager

**Files:**
- Modify: `hvantk/data/alphagenome_streamer.py`
- Modify: `hvantk/tests/test_alphagenome_streamer.py`

- [ ] **Step 1: Write failing tests for checkpoint logic**

Append to `hvantk/tests/test_alphagenome_streamer.py`:

```python
class TestCheckpointManager:
    def test_fresh_start_no_checkpoints(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        mgr = CheckpointManager(str(tmp_path / "out"))
        assert mgr.completed_intervals == []
        assert mgr.failed_variants == []

    def test_save_and_reload_state(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        mgr.mark_interval_complete("chr1:100000-200000")
        mgr.record_failed_variant("chr1", 150000, "A", "T", "timeout")
        mgr.save_state()

        # Reload from disk
        mgr2 = CheckpointManager(out_dir)
        assert "chr1:100000-200000" in mgr2.completed_intervals
        assert len(mgr2.failed_variants) == 1
        assert mgr2.failed_variants[0]["chrom"] == "chr1"

    def test_save_batch_data(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        batch_data = {"variant_1": {"rna_seq": [1.0, 2.0]}}
        mgr.save_batch(0, batch_data)

        batch_path = os.path.join(out_dir, "_checkpoints", "batch_000.json")
        assert os.path.isfile(batch_path)
        with open(batch_path) as f:
            loaded = json.load(f)
        assert loaded == batch_data

    def test_is_interval_complete(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        mgr = CheckpointManager(str(tmp_path / "out"))
        mgr.mark_interval_complete("chr1:100-200")
        assert mgr.is_interval_complete("chr1:100-200") is True
        assert mgr.is_interval_complete("chr2:100-200") is False

    def test_clear_checkpoints(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        mgr.mark_interval_complete("chr1:100-200")
        mgr.save_batch(0, {"data": True})
        mgr.save_state()
        mgr.clear()
        assert mgr.completed_intervals == []
        assert not os.path.isfile(
            os.path.join(out_dir, "_checkpoints", "state.json")
        )
```

Add `import json` at top of test file.

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestCheckpointManager -v`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement CheckpointManager**

Add to `hvantk/data/alphagenome_streamer.py` after `_emit_subgroup`:

```python
class CheckpointManager:
    """Manages checkpoint state for resumable AlphaGenome API runs.

    Tracks which intervals have been processed and which variants failed,
    persisting state to a _checkpoints directory under the output path.
    """

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self._checkpoint_dir = os.path.join(output_dir, "_checkpoints")
        self._state_path = os.path.join(self._checkpoint_dir, "state.json")
        self.completed_intervals: List[str] = []
        self.failed_variants: List[Dict[str, Any]] = []
        self._load_existing_state()

    def _load_existing_state(self) -> None:
        """Load state from disk if a previous checkpoint exists."""
        if os.path.isfile(self._state_path):
            with open(self._state_path) as f:
                state = json.load(f)
            self.completed_intervals = state.get("completed_intervals", [])
            self.failed_variants = state.get("failed_variants", [])
            logger.info(
                f"Resumed checkpoint: {len(self.completed_intervals)} intervals complete, "
                f"{len(self.failed_variants)} failed variants"
            )

    def save_state(self) -> None:
        """Persist current state to disk."""
        os.makedirs(self._checkpoint_dir, exist_ok=True)
        state = {
            "completed_intervals": self.completed_intervals,
            "failed_variants": self.failed_variants,
        }
        with open(self._state_path, "w") as f:
            json.dump(state, f, indent=2)

    def save_batch(self, batch_index: int, data: Any) -> None:
        """Save raw batch data to a numbered JSON file."""
        os.makedirs(self._checkpoint_dir, exist_ok=True)
        batch_path = os.path.join(
            self._checkpoint_dir, f"batch_{batch_index:03d}.json"
        )
        with open(batch_path, "w") as f:
            json.dump(data, f)

    def mark_interval_complete(self, interval_key: str) -> None:
        """Mark an interval as successfully processed."""
        self.completed_intervals.append(interval_key)

    def is_interval_complete(self, interval_key: str) -> bool:
        """Check whether an interval has already been processed."""
        return interval_key in self.completed_intervals

    def record_failed_variant(
        self, chrom: str, pos: int, ref: str, alt: str, reason: str
    ) -> None:
        """Record a variant that could not be processed."""
        self.failed_variants.append(
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "reason": reason}
        )

    def clear(self) -> None:
        """Remove all checkpoint state from memory and disk."""
        self.completed_intervals = []
        self.failed_variants = []
        if os.path.isdir(self._checkpoint_dir):
            import shutil
            shutil.rmtree(self._checkpoint_dir)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestCheckpointManager -v`
Expected: All 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/data/alphagenome_streamer.py hvantk/tests/test_alphagenome_streamer.py
git commit -m "feat(alphagenome): add CheckpointManager for resumable runs"
```

---

### Task 4: Rate-Limited API Caller

**Files:**
- Modify: `hvantk/data/alphagenome_streamer.py`
- Modify: `hvantk/tests/test_alphagenome_streamer.py`

- [ ] **Step 1: Write failing tests for rate-limited calling**

Append to `hvantk/tests/test_alphagenome_streamer.py`:

```python
from unittest.mock import MagicMock, patch


class TestRateLimitedCaller:
    def _make_api_config(self, max_retries=3, retry_backoff=0.01, request_timeout=5):
        return {
            "api": {
                "key": "test-key",
                "max_retries": max_retries,
                "retry_backoff": retry_backoff,
                "request_timeout": request_timeout,
            }
        }

    def test_successful_call(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        mock_model.predict_variant.return_value = {"rna_seq": [1.0]}
        caller = RateLimitedCaller(mock_model, self._make_api_config())

        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=["UBERON:0001157"], output_types=["RNA_SEQ"],
        )
        assert result == {"rna_seq": [1.0]}
        mock_model.predict_variant.assert_called_once()

    def test_retry_on_transient_error(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        # Fail twice, succeed third time
        mock_model.predict_variant.side_effect = [
            Exception("503 Server Error"),
            Exception("503 Server Error"),
            {"rna_seq": [1.0]},
        ]
        caller = RateLimitedCaller(mock_model, self._make_api_config())
        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=[], output_types=[],
        )
        assert result == {"rna_seq": [1.0]}
        assert mock_model.predict_variant.call_count == 3

    def test_max_retries_exceeded_returns_none(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        mock_model.predict_variant.side_effect = Exception("503 Server Error")
        caller = RateLimitedCaller(
            mock_model, self._make_api_config(max_retries=2)
        )
        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=[], output_types=[],
        )
        assert result is None
        assert mock_model.predict_variant.call_count == 2

    def test_cooldown_after_consecutive_rate_limits(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        caller = RateLimitedCaller(mock_model, self._make_api_config())
        # Simulate 3 consecutive rate-limit hits
        caller._consecutive_rate_limits = 3
        assert caller._should_cooldown() is True
        caller._consecutive_rate_limits = 2
        assert caller._should_cooldown() is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestRateLimitedCaller -v`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement RateLimitedCaller**

Add to `hvantk/data/alphagenome_streamer.py` after `CheckpointManager`:

```python
class RateLimitedCaller:
    """Wraps AlphaGenome model with retry logic and adaptive throttling.

    Retries transient errors (429, 5xx) with exponential backoff + jitter.
    Enters cooldown after consecutive rate-limit hits.
    """

    _COOLDOWN_THRESHOLD = 3
    _COOLDOWN_SECONDS = 60.0
    _BASE_DELAY = 0.5

    def __init__(self, model: Any, config: Dict[str, Any]):
        self._model = model
        api_cfg = config["api"]
        self._max_retries = api_cfg.get("max_retries", 3)
        self._retry_backoff = api_cfg.get("retry_backoff", 2.0)
        self._request_timeout = api_cfg.get("request_timeout", 120)
        self._consecutive_rate_limits = 0
        self._total_calls = 0
        self._total_failures = 0

    def _should_cooldown(self) -> bool:
        return self._consecutive_rate_limits >= self._COOLDOWN_THRESHOLD

    def _is_rate_limit_error(self, exc: Exception) -> bool:
        msg = str(exc).lower()
        return "429" in msg or "rate limit" in msg

    def _is_transient_error(self, exc: Exception) -> bool:
        msg = str(exc).lower()
        return (
            self._is_rate_limit_error(exc)
            or "500" in msg
            or "502" in msg
            or "503" in msg
            or "504" in msg
            or "server error" in msg
            or "timeout" in msg
            or "connection" in msg
        )

    def call_predict_variant(
        self,
        interval: Any,
        variant: Any,
        ontology_terms: List[str],
        output_types: List[str],
    ) -> Optional[Any]:
        """Call model.predict_variant with retry and backoff.

        Returns the API response on success, or None if all retries fail.
        """
        if self._should_cooldown():
            logger.warning(
                f"Cooldown: {self._COOLDOWN_THRESHOLD} consecutive rate limits. "
                f"Pausing {self._COOLDOWN_SECONDS}s."
            )
            time.sleep(self._COOLDOWN_SECONDS)
            self._consecutive_rate_limits = 0

        for attempt in range(1, self._max_retries + 1):
            try:
                self._total_calls += 1
                result = self._model.predict_variant(
                    interval=interval,
                    variant=variant,
                    ontology_terms=ontology_terms,
                    requested_outputs=output_types,
                )
                self._consecutive_rate_limits = 0
                time.sleep(self._BASE_DELAY)
                return result

            except Exception as exc:
                if self._is_rate_limit_error(exc):
                    self._consecutive_rate_limits += 1

                if self._is_transient_error(exc) and attempt < self._max_retries:
                    delay = self._retry_backoff ** attempt * (1 + random.random() * 0.5)
                    logger.warning(
                        f"Transient error (attempt {attempt}/{self._max_retries}): "
                        f"{exc}. Retrying in {delay:.1f}s."
                    )
                    time.sleep(delay)
                else:
                    self._total_failures += 1
                    logger.error(
                        f"Failed after {attempt} attempt(s): {exc}"
                    )
                    return None

    @property
    def stats(self) -> Dict[str, int]:
        return {
            "total_calls": self._total_calls,
            "total_failures": self._total_failures,
        }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestRateLimitedCaller -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/data/alphagenome_streamer.py hvantk/tests/test_alphagenome_streamer.py
git commit -m "feat(alphagenome): add RateLimitedCaller with retry and backoff"
```

---

### Task 5: AlphaGenomeStreamer Class

**Files:**
- Modify: `hvantk/data/alphagenome_streamer.py`
- Modify: `hvantk/tests/test_alphagenome_streamer.py`

- [ ] **Step 1: Write failing tests for the streamer**

Append to `hvantk/tests/test_alphagenome_streamer.py`:

```python
class TestAlphaGenomeStreamer:
    def _make_variant_tsv(self, tmp_path, variants):
        """Write a TSV with chrom/pos/ref/alt columns."""
        tsv_path = tmp_path / "variants.tsv"
        lines = ["chrom\tpos\tref\talt\n"]
        for v in variants:
            lines.append(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\n")
        tsv_path.write_text("".join(lines))
        return str(tsv_path)

    def _make_config_file(self, tmp_path, overrides=None):
        cfg = {
            "api": {"key": "test-key", "max_retries": 2,
                    "retry_backoff": 0.01, "request_timeout": 5},
            "ontology": {"terms": ["UBERON:0001157"],
                         "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 100_000, "adaptive": False,
                          "adaptive_max_size": 100_000, "density_window": 50_000},
        }
        if overrides:
            for section, vals in overrides.items():
                cfg.setdefault(section, {}).update(vals)
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        return str(cfg_path)

    @patch("hvantk.data.alphagenome_streamer._create_dna_client")
    def test_setup_loads_tsv_input(self, mock_create_client, tmp_path):
        from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

        mock_create_client.return_value = MagicMock()
        tsv = self._make_variant_tsv(tmp_path, [("chr1", 500000, "A", "T")])
        cfg = self._make_config_file(tmp_path)
        out = str(tmp_path / "output")

        streamer = AlphaGenomeStreamer(
            input_path=tsv, output_dir=out, config_path=cfg,
        )
        streamer.setup()
        assert len(streamer._variants) == 1
        assert streamer._variants[0].chrom == "chr1"
        streamer.teardown()

    @patch("hvantk.data.alphagenome_streamer._create_dna_client")
    def test_stream_calls_api_per_variant(self, mock_create_client, tmp_path):
        from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

        mock_model = MagicMock()
        mock_result = MagicMock()
        mock_result.reference = MagicMock()
        mock_result.alternate = MagicMock()
        mock_model.predict_variant.return_value = mock_result
        mock_create_client.return_value = mock_model

        tsv = self._make_variant_tsv(
            tmp_path,
            [("chr1", 500000, "A", "T"), ("chr2", 600000, "G", "C")],
        )
        cfg = self._make_config_file(tmp_path)
        out = str(tmp_path / "output")

        streamer = AlphaGenomeStreamer(
            input_path=tsv, output_dir=out, config_path=cfg,
        )
        streamer.setup()
        batches = list(streamer.stream())
        # With adaptive=False, each variant is its own interval
        assert mock_model.predict_variant.call_count == 2
        streamer.teardown()

    @patch("hvantk.data.alphagenome_streamer._create_dna_client")
    def test_no_resume_clears_checkpoints(self, mock_create_client, tmp_path):
        from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

        mock_create_client.return_value = MagicMock()
        tsv = self._make_variant_tsv(tmp_path, [("chr1", 500000, "A", "T")])
        cfg = self._make_config_file(tmp_path)
        out = str(tmp_path / "output")

        # Create a fake checkpoint
        ckpt_dir = os.path.join(out, "_checkpoints")
        os.makedirs(ckpt_dir, exist_ok=True)
        state_path = os.path.join(ckpt_dir, "state.json")
        with open(state_path, "w") as f:
            json.dump({"completed_intervals": ["chr1:450000-550000"], "failed_variants": []}, f)

        streamer = AlphaGenomeStreamer(
            input_path=tsv, output_dir=out, config_path=cfg, no_resume=True,
        )
        streamer.setup()
        assert not os.path.isfile(state_path)
        streamer.teardown()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestAlphaGenomeStreamer -v`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement AlphaGenomeStreamer**

Add to `hvantk/data/alphagenome_streamer.py` after `RateLimitedCaller`:

```python
import csv as csv_module


def _create_dna_client(api_key: str) -> Any:
    """Create an AlphaGenome DNA client. Isolated for mocking."""
    try:
        from alphagenome.models import dna_client
    except ImportError:
        raise ImportError(
            "alphagenome package not installed. "
            "Install with: pip install alphagenome"
        )
    return dna_client.create(api_key)


def _load_variants_from_tsv(tsv_path: str) -> List[VariantRecord]:
    """Load variants from a TSV file with chrom/pos/ref/alt columns."""
    variants = []
    with open(tsv_path) as f:
        reader = csv_module.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append(VariantRecord(
                chrom=row["chrom"],
                pos=int(row["pos"]),
                ref=row["ref"],
                alt=row["alt"],
            ))
    return variants


def _interval_key(interval: GenomicInterval) -> str:
    """Create a string key for checkpoint tracking."""
    return f"{interval.chrom}:{interval.start}-{interval.end}"


class AlphaGenomeStreamer(HailDataStreamer):
    """Streams variant predictions from the AlphaGenome API.

    Reads variants from a Hail Table or TSV, groups them into genomic
    intervals, calls the AlphaGenome predict_variant API with rate
    limiting and checkpoint-based resumption, and assembles per-modality
    Hail Tables as output.

    Parameters
    ----------
    input_path : str
        Path to a Hail Table (.ht) or TSV file with chrom/pos/ref/alt columns.
    output_dir : str
        Directory for output Hail Tables and checkpoints.
    config_path : str
        Path to AlphaGenome YAML config file.
    no_resume : bool
        If True, clear existing checkpoints and start fresh.
    chunk_size : int
        Number of intervals to process per batch before checkpointing.
    """

    def __init__(
        self,
        input_path: str,
        output_dir: str,
        config_path: str,
        no_resume: bool = False,
        chunk_size: int = 50,
    ):
        super().__init__(name="AlphaGenomeStreamer", chunk_size=chunk_size,
                         init_hail=input_path.endswith(".ht"))
        self.input_path = input_path
        self.output_dir = output_dir
        self.config_path = config_path
        self.no_resume = no_resume

        self._config: Dict[str, Any] = {}
        self._model: Any = None
        self._caller: Optional[RateLimitedCaller] = None
        self._variants: List[VariantRecord] = []
        self._interval_groups: List[Tuple[GenomicInterval, List[VariantRecord]]] = []
        self._checkpoint: Optional[CheckpointManager] = None
        self._start_time: float = 0.0

    def setup(self) -> None:
        """Initialize Hail, load config, authenticate, load variants, compute intervals."""
        super().setup()
        self._start_time = time.time()

        # Load config
        self._config = load_config(self.config_path)
        self.logger.info(f"Config loaded from {self.config_path}")

        # Authenticate
        self._model = _create_dna_client(self._config["api"]["key"])
        self._caller = RateLimitedCaller(self._model, self._config)
        self.logger.info("AlphaGenome client authenticated")

        # Handle checkpoints
        if self.no_resume:
            # Clear any existing checkpoints
            mgr = CheckpointManager(self.output_dir)
            mgr.clear()
            self.logger.info("Cleared existing checkpoints (--no-resume)")
        self._checkpoint = CheckpointManager(self.output_dir)

        # Load variants
        if self.input_path.endswith(".ht"):
            self._variants = self._load_variants_from_hail_table()
        else:
            self._variants = _load_variants_from_tsv(self.input_path)
        self.logger.info(f"Loaded {len(self._variants)} variants from {self.input_path}")

        # Compute intervals
        self._interval_groups = compute_intervals(self._variants, self._config)
        self.logger.info(f"Computed {len(self._interval_groups)} intervals")

        # Filter already-completed intervals
        pending = [
            (iv, vs) for iv, vs in self._interval_groups
            if not self._checkpoint.is_interval_complete(_interval_key(iv))
        ]
        skipped = len(self._interval_groups) - len(pending)
        if skipped > 0:
            self.logger.info(f"Resuming: skipping {skipped} completed intervals")
        self._interval_groups = pending

    def _load_variants_from_hail_table(self) -> List[VariantRecord]:
        """Load variants from a Hail Table keyed by (locus, alleles)."""
        import hail as hl

        ht = hl.read_table(self.input_path)
        rows = ht.select(
            chrom=ht.locus.contig,
            pos=ht.locus.position,
            ref=ht.alleles[0],
            alt=ht.alleles[1],
        ).collect()
        return [
            VariantRecord(chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt)
            for r in rows
        ]

    def stream(self) -> Iterator[Dict[str, Any]]:
        """Stream predictions batch by batch, checkpointing after each."""
        try:
            from alphagenome.data import genome as ag_genome
            from alphagenome.models import dna_client as ag_client
        except ImportError:
            raise ImportError(
                "alphagenome package not installed. "
                "Install with: pip install alphagenome"
            )

        ontology_terms = self._config["ontology"].get("terms", [])
        output_types_raw = self._config["ontology"].get("output_types", [])
        output_types = [getattr(ag_client.OutputType, ot) for ot in output_types_raw]

        total = len(self._interval_groups)
        processed_variants = 0

        for batch_idx in range(0, total, self.chunk_size):
            batch = self._interval_groups[batch_idx:batch_idx + self.chunk_size]
            batch_results: Dict[str, List[Dict[str, Any]]] = {}

            for interval, variants in batch:
                ag_interval = ag_genome.Interval(
                    chromosome=interval.chrom,
                    start=interval.start,
                    end=interval.end,
                )
                for v in variants:
                    ag_variant = ag_genome.Variant(
                        chromosome=v.chrom,
                        position=v.pos,
                        reference_bases=v.ref,
                        alternate_bases=v.alt,
                    )
                    result = self._caller.call_predict_variant(
                        interval=ag_interval,
                        variant=ag_variant,
                        ontology_terms=ontology_terms,
                        output_types=output_types,
                    )
                    variant_key = f"{v.chrom}:{v.pos}:{v.ref}>{v.alt}"
                    if result is not None:
                        batch_results[variant_key] = result
                        processed_variants += 1
                    else:
                        self._checkpoint.record_failed_variant(
                            v.chrom, v.pos, v.ref, v.alt, "max retries exceeded"
                        )

                self._checkpoint.mark_interval_complete(_interval_key(interval))

            # Checkpoint this batch
            batch_num = batch_idx // self.chunk_size
            self._checkpoint.save_batch(batch_num, {
                k: str(v) for k, v in batch_results.items()
            })
            self._checkpoint.save_state()

            elapsed = time.time() - self._start_time
            self.logger.info(
                f"Batch {batch_num}: {len(batch)} intervals, "
                f"{processed_variants} variants processed, "
                f"elapsed {elapsed:.0f}s"
            )

            yield batch_results

    def process_chunk(self, chunk: Any) -> Any:
        """Identity — processing happens in stream()."""
        return chunk

    def teardown(self) -> None:
        """Log final summary."""
        elapsed = time.time() - self._start_time
        stats = self._caller.stats if self._caller else {}
        failed_count = len(self._checkpoint.failed_variants) if self._checkpoint else 0
        self.logger.info(
            f"AlphaGenome run complete. "
            f"API calls: {stats.get('total_calls', 0)}, "
            f"failures: {stats.get('total_failures', 0)}, "
            f"failed variants: {failed_count}, "
            f"runtime: {elapsed:.0f}s"
        )
        super().teardown()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py::TestAlphaGenomeStreamer -v`
Expected: All 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/data/alphagenome_streamer.py hvantk/tests/test_alphagenome_streamer.py
git commit -m "feat(alphagenome): add AlphaGenomeStreamer class"
```

---

### Task 6: Constants, Exports, and Builder Integration

**Files:**
- Modify: `hvantk/core/constants.py`
- Modify: `hvantk/data/__init__.py`
- Modify: `hvantk/tables/table_builders.py`
- Modify: `hvantk/tables/registry.py`

- [ ] **Step 1: Add AlphaGenome constants**

Add to `hvantk/core/constants.py` before the `__all__` list:

```python
# AlphaGenome defaults
ALPHAGENOME_DEFAULT_INTERVAL_SIZE = 1_048_576  # 1Mbp
ALPHAGENOME_DEFAULT_DENSITY_WINDOW = 50_000    # 50kb
ALPHAGENOME_DEFAULT_RETRY_BACKOFF = 2.0
ALPHAGENOME_DEFAULT_MAX_RETRIES = 3
ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT = 120
```

Add to the `__all__` list:

```python
    "ALPHAGENOME_DEFAULT_INTERVAL_SIZE",
    "ALPHAGENOME_DEFAULT_DENSITY_WINDOW",
    "ALPHAGENOME_DEFAULT_RETRY_BACKOFF",
    "ALPHAGENOME_DEFAULT_MAX_RETRIES",
    "ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT",
```

- [ ] **Step 2: Export AlphaGenomeStreamer from data package**

Add to `hvantk/data/__init__.py`:

Import line:
```python
from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer
```

Add `"AlphaGenomeStreamer"` to the `__all__` list.

- [ ] **Step 3: Add builder function to table_builders.py**

Add at the end of `hvantk/tables/table_builders.py`:

```python
def create_alphagenome_tb(
    input_path: str,
    output_path: str,
    config_path: str,
    no_resume: bool = False,
    overwrite: bool = False,
) -> None:
    """Build per-modality Hail Tables from AlphaGenome variant predictions.

    Runs the AlphaGenomeStreamer to call the API for each variant, then
    writes per-modality Hail Tables to output_path (a directory).

    Parameters
    ----------
    input_path : str
        Path to Hail Table (.ht) or TSV with chrom/pos/ref/alt columns.
    output_path : str
        Output directory for per-modality Hail Tables.
    config_path : str
        Path to AlphaGenome YAML config file.
    no_resume : bool
        If True, discard existing checkpoints and restart.
    overwrite : bool
        If True, overwrite existing output tables.
    """
    from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

    if overwrite and os.path.isdir(output_path):
        import shutil
        shutil.rmtree(output_path)

    streamer = AlphaGenomeStreamer(
        input_path=input_path,
        output_dir=output_path,
        config_path=config_path,
        no_resume=no_resume,
    )
    streamer.setup()
    try:
        for _batch in streamer.stream():
            pass  # checkpointing handled internally
    finally:
        streamer.teardown()
```

- [ ] **Step 4: Register in registry**

Add to `hvantk/tables/registry.py` inside the `TABLE_BUILDERS` dict:

```python
    "alphagenome": create_table_adapter(
        "hvantk.tables.table_builders", "create_alphagenome_tb"
    ),
```

- [ ] **Step 5: Commit**

```bash
git add hvantk/core/constants.py hvantk/data/__init__.py hvantk/tables/table_builders.py hvantk/tables/registry.py
git commit -m "feat(alphagenome): add constants, exports, builder, and registry entry"
```

---

### Task 7: CLI Command

**Files:**
- Modify: `hvantk/commands/make_table_cli.py`

- [ ] **Step 1: Add the CLI command**

Add to `hvantk/commands/make_table_cli.py` after the existing commands:

```python
@mktable_group.command("alphagenome")
@click.option(
    "--input", "input_path", required=True, type=str,
    help="Path to Hail Table (.ht) or TSV with chrom/pos/ref/alt columns",
)
@click.option(
    "--output-dir", required=True, type=str,
    help="Output directory for per-modality Hail Tables",
)
@click.option(
    "--config", "config_path", required=True, type=str,
    help="Path to AlphaGenome YAML config file",
)
@click.option(
    "--no-resume", is_flag=True,
    help="Discard existing checkpoints and restart from scratch",
)
@_overwrite_opt
def mktable_alphagenome(input_path, output_dir, config_path, no_resume, overwrite):
    """Build per-modality Hail Tables from AlphaGenome variant predictions."""
    from hvantk.tables.table_builders import create_alphagenome_tb

    logger.info("Building AlphaGenome prediction tables")
    create_alphagenome_tb(
        input_path=input_path,
        output_path=output_dir,
        config_path=config_path,
        no_resume=no_resume,
        overwrite=overwrite,
    )
    click.echo(f"AlphaGenome tables created at {output_dir}")
```

- [ ] **Step 2: Verify CLI loads without error**

Run: `python -c "from hvantk.commands.make_table_cli import mktable_group; print([c.name for c in mktable_group.commands.values()])"`
Expected: Output includes `'alphagenome'` in the list

- [ ] **Step 3: Commit**

```bash
git add hvantk/commands/make_table_cli.py
git commit -m "feat(alphagenome): add mktable alphagenome CLI command"
```

---

### Task 8: Final Integration Test and Cleanup

**Files:**
- Modify: `hvantk/tests/test_alphagenome_streamer.py`

- [ ] **Step 1: Add an end-to-end integration test (mocked API)**

Append to `hvantk/tests/test_alphagenome_streamer.py`:

```python
class TestEndToEnd:
    """End-to-end test with mocked AlphaGenome API."""

    @patch("hvantk.data.alphagenome_streamer._create_dna_client")
    def test_full_pipeline_tsv_input(self, mock_create_client, tmp_path):
        from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

        # Setup mock
        mock_model = MagicMock()
        mock_result = MagicMock()
        mock_result.reference = MagicMock()
        mock_result.alternate = MagicMock()
        mock_model.predict_variant.return_value = mock_result
        mock_create_client.return_value = mock_model

        # Write variant TSV
        tsv_path = tmp_path / "variants.tsv"
        tsv_path.write_text(
            "chrom\tpos\tref\talt\n"
            "chr1\t100000\tA\tT\n"
            "chr1\t120000\tG\tC\n"
            "chr2\t500000\tC\tG\n"
        )

        # Write config
        cfg = {
            "api": {"key": "test-key", "max_retries": 2,
                    "retry_backoff": 0.01, "request_timeout": 5},
            "ontology": {"terms": ["UBERON:0001157"],
                         "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 100_000, "adaptive": True,
                          "adaptive_max_size": 100_000, "density_window": 50_000},
        }
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text(yaml.dump(cfg))

        out_dir = str(tmp_path / "output")

        # Run streamer
        streamer = AlphaGenomeStreamer(
            input_path=str(tsv_path),
            output_dir=out_dir,
            config_path=str(cfg_path),
        )
        streamer.setup()
        batches = list(streamer.stream())
        streamer.teardown()

        # Verify
        assert mock_model.predict_variant.call_count == 3
        assert len(batches) >= 1
        # Checkpoint state should exist
        state_path = os.path.join(out_dir, "_checkpoints", "state.json")
        assert os.path.isfile(state_path)
        with open(state_path) as f:
            state = json.load(f)
        assert len(state["completed_intervals"]) > 0

    @patch("hvantk.data.alphagenome_streamer._create_dna_client")
    def test_resume_skips_completed(self, mock_create_client, tmp_path):
        from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

        mock_model = MagicMock()
        mock_result = MagicMock()
        mock_result.reference = MagicMock()
        mock_result.alternate = MagicMock()
        mock_model.predict_variant.return_value = mock_result
        mock_create_client.return_value = mock_model

        tsv_path = tmp_path / "variants.tsv"
        tsv_path.write_text(
            "chrom\tpos\tref\talt\n"
            "chr1\t100000\tA\tT\n"
            "chr2\t500000\tC\tG\n"
        )
        cfg = {
            "api": {"key": "test-key", "max_retries": 2,
                    "retry_backoff": 0.01, "request_timeout": 5},
            "ontology": {"terms": ["UBERON:0001157"],
                         "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 100_000, "adaptive": False,
                          "adaptive_max_size": 100_000, "density_window": 50_000},
        }
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        out_dir = str(tmp_path / "output")

        # Run first time
        streamer = AlphaGenomeStreamer(
            input_path=str(tsv_path), output_dir=out_dir,
            config_path=str(cfg_path),
        )
        streamer.setup()
        list(streamer.stream())
        streamer.teardown()
        first_call_count = mock_model.predict_variant.call_count

        # Run again — should skip all intervals
        mock_model.predict_variant.reset_mock()
        streamer2 = AlphaGenomeStreamer(
            input_path=str(tsv_path), output_dir=out_dir,
            config_path=str(cfg_path),
        )
        streamer2.setup()
        batches = list(streamer2.stream())
        streamer2.teardown()

        assert mock_model.predict_variant.call_count == 0
        assert len(batches) == 0
```

- [ ] **Step 2: Run the full test suite**

Run: `pytest hvantk/tests/test_alphagenome_streamer.py -v`
Expected: All tests PASS

- [ ] **Step 3: Commit**

```bash
git add hvantk/tests/test_alphagenome_streamer.py
git commit -m "test(alphagenome): add end-to-end integration tests"
```

- [ ] **Step 4: Run full project test suite (fast tests only)**

Run: `pytest -q`
Expected: All existing tests still PASS, new tests PASS
