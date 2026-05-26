"""AlphaGenome variant prediction streamer.

Streams variant positions through the AlphaGenome API and produces
multimodal variant effect predictions. Outputs are currently written as
consolidated JSON (predictions.json); per-modality Hail Table assembly
will be added once the AlphaGenome SDK response structure is validated.
"""

import csv as csv_module
import json
import logging
import os
import random
import re
import time
from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional, Tuple

import yaml

from hvantk.core.constants import (
    ALPHAGENOME_DEFAULT_INTERVAL_SIZE,
    ALPHAGENOME_DEFAULT_DENSITY_WINDOW,
    ALPHAGENOME_DEFAULT_MAX_RETRIES,
    ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT,
    ALPHAGENOME_DEFAULT_RETRY_BACKOFF,
)
from hvantk.core.utils.streaming import HailDataStreamer

logger = logging.getLogger(__name__)

_REQUIRED_SECTIONS = ("api", "ontology")
_DEFAULT_INTERVALS = {
    "default_size": ALPHAGENOME_DEFAULT_INTERVAL_SIZE,
    "adaptive": True,
    "adaptive_max_size": ALPHAGENOME_DEFAULT_INTERVAL_SIZE,
    "density_window": ALPHAGENOME_DEFAULT_DENSITY_WINDOW,
}
_CHR_X_ORDER = 23
_CHR_Y_ORDER = 24
_CHR_M_ORDER = 25


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

    if not isinstance(config, dict):
        raise ValueError(
            f"Config file is empty or malformed (expected a YAML mapping, "
            f"got {type(config).__name__}): {config_path}"
        )

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
    def _chrom_sort_key(chrom: str) -> Tuple[int, Any]:
        c = chrom.lower()
        if c.startswith("chr"):
            c = c[3:]
        if c.isdigit():
            return (0, int(c))
        if c == "x":
            return (1, _CHR_X_ORDER)
        if c == "y":
            return (1, _CHR_Y_ORDER)
        if c in {"m", "mt"}:
            return (1, _CHR_M_ORDER)
        m = re.match(r"^(\d+)(.*)$", c)
        if m:
            return (2, int(m.group(1)), m.group(2))
        return (3, c)

    sorted_variants = sorted(variants, key=lambda v: (_chrom_sort_key(v.chrom), v.pos))

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
        self.completed_intervals: set = set()
        self.failed_variants: List[Dict[str, Any]] = []
        self._load_existing_state()

    def _load_existing_state(self) -> None:
        if os.path.isfile(self._state_path):
            with open(self._state_path) as f:
                state = json.load(f)
            self.completed_intervals = set(state.get("completed_intervals", []))
            self.failed_variants = state.get("failed_variants", [])
            logger.info(
                f"Resumed checkpoint: {len(self.completed_intervals)} intervals complete, "
                f"{len(self.failed_variants)} failed variants"
            )

    def save_state(self) -> None:
        os.makedirs(self._checkpoint_dir, exist_ok=True)
        state = {
            "completed_intervals": sorted(self.completed_intervals),
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
        self.completed_intervals.add(interval_key)

    def is_interval_complete(self, interval_key: str) -> bool:
        return interval_key in self.completed_intervals

    def record_failed_variant(
        self, chrom: str, pos: int, ref: str, alt: str, reason: str
    ) -> None:
        self.failed_variants.append(
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "reason": reason}
        )

    def clear(self) -> None:
        self.completed_intervals = set()
        self.failed_variants = []
        if os.path.isdir(self._checkpoint_dir):
            import shutil

            shutil.rmtree(self._checkpoint_dir)


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
        self._max_retries = api_cfg.get("max_retries", ALPHAGENOME_DEFAULT_MAX_RETRIES)
        self._retry_backoff = api_cfg.get("retry_backoff", ALPHAGENOME_DEFAULT_RETRY_BACKOFF)
        self._request_timeout = api_cfg.get("request_timeout", ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT)
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
                    delay = self._retry_backoff ** attempt * (
                        1 + random.random() * 0.5
                    )
                    logger.warning(
                        f"Transient error (attempt {attempt}/{self._max_retries}): "
                        f"{exc}. Retrying in {delay:.1f}s."
                    )
                    time.sleep(delay)
                else:
                    self._total_failures += 1
                    logger.error(f"Failed after {attempt} attempt(s): {exc}")
                    return None

    @property
    def stats(self) -> Dict[str, int]:
        return {
            "total_calls": self._total_calls,
            "total_failures": self._total_failures,
        }


def _import_alphagenome() -> Tuple[Any, Any]:
    """Import alphagenome modules. Isolated for mocking."""
    try:
        from alphagenome.data import genome as ag_genome
        from alphagenome.models import dna_client as ag_client
    except ImportError:
        raise ImportError(
            "alphagenome package not installed. "
            "Install with: pip install alphagenome"
        )
    return ag_genome, ag_client


def _create_dna_client(api_key: str, timeout: Optional[float] = None) -> Any:
    """Create an AlphaGenome DNA client. Isolated for mocking."""
    _, ag_client = _import_alphagenome()
    kwargs: Dict[str, Any] = {}
    if timeout is not None:
        kwargs["timeout"] = timeout
    return ag_client.create(api_key, **kwargs)


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


def _serialize_value(obj: Any) -> Any:
    """Recursively convert an SDK object to a JSON-serializable structure.

    Handles numpy arrays, pandas DataFrames, dataclass-like SDK objects
    (TrackData, Interval, Output), and plain Python types.
    """
    if obj is None:
        return None
    # Plain JSON types
    if isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, (list, tuple)):
        return [_serialize_value(item) for item in obj]
    if isinstance(obj, dict):
        return {str(k): _serialize_value(v) for k, v in obj.items()}
    if obj.__class__.__module__.startswith("unittest.mock"):
        return str(obj)
    # numpy array-like → nested list
    tolist = getattr(obj, "tolist", None)
    if callable(tolist):
        try:
            return _serialize_value(tolist())
        except Exception:
            pass
    # pandas DataFrame-like → list of dicts
    to_dict = getattr(obj, "to_dict", None)
    iterrows = getattr(obj, "iterrows", None)
    if callable(to_dict) and callable(iterrows):
        try:
            return _serialize_value(to_dict(orient="records"))
        except Exception:
            pass
    # Dataclass / SDK objects — recurse into public attributes
    if hasattr(obj, "__dict__"):
        return {
            k: _serialize_value(v)
            for k, v in vars(obj).items()
            if not k.startswith("_") and not callable(v)
        }
    # Fallback
    return str(obj)


def _serialize_prediction(result: Any) -> Dict[str, Any]:
    """Serialize an AlphaGenome VariantOutput for JSON checkpointing.

    Recursively converts the reference/alternate Output objects —
    including nested TrackData (numpy arrays, pandas metadata) and
    Interval objects — into JSON-safe dicts.
    """
    try:
        if hasattr(result, "reference") and hasattr(result, "alternate"):
            return {
                "reference": _serialize_value(result.reference),
                "alternate": _serialize_value(result.alternate),
            }
        return {"raw": str(result)}
    except Exception:
        return {"raw": str(result)}


class AlphaGenomeStreamer(HailDataStreamer):
    """Streams variant predictions from the AlphaGenome API.

    Reads variants from a Hail Table or TSV, groups them into genomic
    intervals, calls the AlphaGenome predict_variant API with rate
    limiting and checkpoint-based resumption, and writes consolidated
    JSON predictions. Per-modality Hail Table assembly is deferred
    until the SDK response structure is validated.

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
        super().setup()
        self._start_time = time.time()

        self._config = load_config(self.config_path)
        self.logger.info(f"Config loaded from {self.config_path}")

        self._model = _create_dna_client(
            self._config["api"]["key"],
            timeout=self._config["api"].get("request_timeout"),
        )
        self._caller = RateLimitedCaller(self._model, self._config)
        self.logger.info("AlphaGenome client authenticated")

        if self.no_resume:
            mgr = CheckpointManager(self.output_dir)
            mgr.clear()
            self.logger.info("Cleared existing checkpoints (--no-resume)")
        self._checkpoint = CheckpointManager(self.output_dir)

        if self.input_path.endswith(".ht"):
            self._variants = self._load_variants_from_hail_table()
        else:
            self._variants = _load_variants_from_tsv(self.input_path)
        self.logger.info(f"Loaded {len(self._variants)} variants from {self.input_path}")

        self._interval_groups = compute_intervals(self._variants, self._config)
        self.logger.info(f"Computed {len(self._interval_groups)} intervals")

        pending = [
            (iv, vs) for iv, vs in self._interval_groups
            if not self._checkpoint.is_interval_complete(_interval_key(iv))
        ]
        skipped = len(self._interval_groups) - len(pending)
        if skipped > 0:
            self.logger.info(f"Resuming: skipping {skipped} completed intervals")
        self._interval_groups = pending

    _MAX_VARIANTS_COLLECT = 100_000

    def _load_variants_from_hail_table(self) -> List[VariantRecord]:
        """Load variants from a Hail Table keyed by (locus, alleles).

        Uses .collect() which materializes all rows on the driver.
        Raises ValueError if the table exceeds _MAX_VARIANTS_COLLECT rows
        to prevent driver OOM on unexpectedly large inputs.
        """
        import hail as hl

        ht = hl.read_table(self.input_path)
        n_variants = ht.count()
        if n_variants > self._MAX_VARIANTS_COLLECT:
            raise ValueError(
                f"Input table has {n_variants} variants, exceeding the "
                f"maximum of {self._MAX_VARIANTS_COLLECT} for API-bound "
                f"prediction. Use a smaller variant set or export to TSV."
            )
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
        ag_genome, ag_client = _import_alphagenome()

        ontology_terms = self._config["ontology"].get("terms", [])
        output_types_raw = self._config["ontology"].get("output_types", [])

        # Validate output types against SDK before starting long-running job
        output_types = []
        for ot in output_types_raw:
            if not hasattr(ag_client.OutputType, ot):
                allowed = [n for n in dir(ag_client.OutputType) if not n.startswith("_")]
                raise ValueError(
                    f"Invalid AlphaGenome OutputType in config: '{ot}'. "
                    f"Allowed values: {', '.join(sorted(allowed))}"
                )
            output_types.append(getattr(ag_client.OutputType, ot))

        total = len(self._interval_groups)
        processed_variants = 0

        for batch_idx in range(0, total, self.chunk_size):
            batch = self._interval_groups[batch_idx:batch_idx + self.chunk_size]
            batch_results: Dict[str, Any] = {}

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

            batch_num = batch_idx // self.chunk_size
            self._checkpoint.save_batch(batch_num, {
                k: _serialize_prediction(v) for k, v in batch_results.items()
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
        """Log summary and assemble per-modality outputs from checkpoints."""
        elapsed = time.time() - self._start_time
        stats = self._caller.stats if self._caller else {}
        failed_count = len(self._checkpoint.failed_variants) if self._checkpoint else 0

        self._assemble_outputs()

        self.logger.info(
            f"AlphaGenome run complete. "
            f"API calls: {stats.get('total_calls', 0)}, "
            f"failures: {stats.get('total_failures', 0)}, "
            f"failed variants: {failed_count}, "
            f"runtime: {elapsed:.0f}s"
        )
        super().teardown()

    def _assemble_outputs(self) -> None:
        """Assemble predictions from batch checkpoint files.

        Merges all batch_NNN.json files into a consolidated predictions.json.
        Full assembly into per-modality Hail Tables requires the AlphaGenome
        SDK to inspect prediction structure and is deferred until the SDK
        is available for validation.
        """
        if not self._checkpoint:
            return

        checkpoint_dir = self._checkpoint._checkpoint_dir
        if not os.path.isdir(checkpoint_dir):
            return

        batch_files = sorted(
            f for f in os.listdir(checkpoint_dir)
            if f.startswith("batch_") and f.endswith(".json")
        )
        if not batch_files:
            return

        all_predictions: Dict[str, Any] = {}
        for batch_file in batch_files:
            batch_path = os.path.join(checkpoint_dir, batch_file)
            with open(batch_path) as f:
                batch_data = json.load(f)
            all_predictions.update(batch_data)

        merged_path = os.path.join(self.output_dir, "predictions.json")
        os.makedirs(self.output_dir, exist_ok=True)
        with open(merged_path, "w") as f:
            json.dump(all_predictions, f, indent=2)
        self.logger.info(
            f"Assembled {len(all_predictions)} variant predictions to {merged_path}"
        )
