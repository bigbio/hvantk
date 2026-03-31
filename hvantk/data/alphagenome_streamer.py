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
    sorted_variants = sorted(variants, key=lambda v: (v.chrom, v.pos))

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

    result: List[Tuple[GenomicInterval, List[Any]]] = []
    for group in groups:
        span_start = group[0].pos
        span_end = group[-1].pos
        span = span_end - span_start

        if span <= max_size:
            midpoint = (span_start + span_end) // 2
            size = max(default_size, span + density_window)
            size = min(size, max_size)
            start, end = _center_interval(midpoint, size)
            interval = GenomicInterval(chrom=group[0].chrom, start=start, end=end)
            result.append((interval, group))
        else:
            sub_group: List[Any] = [group[0]]
            for v in group[1:]:
                if (v.pos - sub_group[0].pos) <= max_size:
                    sub_group.append(v)
                else:
                    _emit_subgroup(
                        sub_group, default_size, max_size, density_window, result
                    )
                    sub_group = [v]
            _emit_subgroup(
                sub_group, default_size, max_size, density_window, result
            )

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


class CheckpointManager:
    """Manages checkpoint state for resumable AlphaGenome API runs."""

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self._checkpoint_dir = os.path.join(output_dir, "_checkpoints")
        self._state_path = os.path.join(self._checkpoint_dir, "state.json")
        self.completed_intervals: List[str] = []
        self.failed_variants: List[Dict[str, Any]] = []
        self._load_existing_state()

    def _load_existing_state(self) -> None:
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
        os.makedirs(self._checkpoint_dir, exist_ok=True)
        state = {
            "completed_intervals": self.completed_intervals,
            "failed_variants": self.failed_variants,
        }
        with open(self._state_path, "w") as f:
            json.dump(state, f, indent=2)

    def save_batch(self, batch_index: int, data: Any) -> None:
        os.makedirs(self._checkpoint_dir, exist_ok=True)
        batch_path = os.path.join(
            self._checkpoint_dir, f"batch_{batch_index:03d}.json"
        )
        with open(batch_path, "w") as f:
            json.dump(data, f)

    def mark_interval_complete(self, interval_key: str) -> None:
        self.completed_intervals.append(interval_key)

    def is_interval_complete(self, interval_key: str) -> bool:
        return interval_key in self.completed_intervals

    def record_failed_variant(
        self, chrom: str, pos: int, ref: str, alt: str, reason: str
    ) -> None:
        self.failed_variants.append(
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "reason": reason}
        )

    def clear(self) -> None:
        self.completed_intervals = []
        self.failed_variants = []
        if os.path.isdir(self._checkpoint_dir):
            import shutil

            shutil.rmtree(self._checkpoint_dir)
