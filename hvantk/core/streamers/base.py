"""Low-level chunked-IO streaming primitives.

This module pre-dates the Phase A artifact contract and is intentionally
NOT part of the artifact pipeline. Streamers handle the byte-level
download → on-disk landing of raw upstream files (e.g. ClinVar VCF chunks,
HGNC TSV chunks) before any parse/build step runs. The artifact pipeline
(``hvantk/core/io`` + ``AnnotationTable``/``ExpressionMatrix``/``GeneSet``)
starts where streamers end: at well-formed on-disk files.

For algorithm-side persistence with provenance, use ``hvantk.core.io.save``
(artifact API) or ``hvantk.core.io.save_native`` (native passthrough), NOT
the streamer persistence methods.

Streamer subclasses (``HailDataStreamer``, etc.) reason in terms of Hail
tables, JSON, text, and bytes because that's the shape of raw upstream
data, not because the streamer layer is platform-architectural. If a
new upstream source needs a streamer, it adds a subclass here; if the
source produces ready-to-build files, no streamer is needed.
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Iterator, Optional, List
import hail as hl
from hvantk.core.utils.hail_context import init_hail, hail_initialized
import os
import json

logger = logging.getLogger(__name__)


class DataStreamer(ABC):
    """
    Abstract base class for data streamers that process genomic data in chunks.
    """

    def __init__(self, name: str, chunk_size: int = 10000):
        self.name = name
        self.chunk_size = chunk_size
        self.logger = logging.getLogger(f"{__name__}.{name}")

    @abstractmethod
    def stream(self) -> Iterator[Any]:
        """
        Stream data in chunks. Must be implemented by subclasses.

        Yields:
            Chunks of processed data
        """
        pass

    @abstractmethod
    def process_chunk(self, chunk: Any) -> Any:
        """
        Process a single chunk of data. Must be implemented by subclasses.

        Args:
            chunk: A chunk of data to process

        Returns:
            Processed chunk
        """
        pass

    def setup(self) -> None:
        """
        Setup method called before streaming starts. Override if needed.
        """
        self.logger.info(f"Setting up {self.name} streamer")

    def teardown(self) -> None:
        """
        Cleanup method called after streaming ends. Override if needed.
        """
        self.logger.info(f"Tearing down {self.name} streamer")


class HailDataStreamer(DataStreamer):
    """
    Base class for Hail-based data streamers.
    """

    def __init__(self, name: str, chunk_size: int = 10000, init_hail: bool = True):
        super().__init__(name, chunk_size)
        self.init_hail = init_hail
        # Track whether THIS streamer triggered initialization (informational only)
        self._streamer_initialized_hail = False

    def setup(self) -> None:
        """Initialize Hail if needed"""
        super().setup()
        if self.init_hail:
            was_already = hail_initialized()
            init_hail()  # idempotent global initializer
            if not was_already and hail_initialized():
                self._streamer_initialized_hail = True
                self.logger.info("Hail initialized (by streamer)")
            else:
                self.logger.debug("Hail already initialized; streamer proceeding")

    def teardown(self) -> None:
        """No-op for global Hail lifecycle (do not stop shared Hail context)."""
        super().teardown()
        # Intentionally NOT calling hl.stop() here to avoid shutting down a shared
        # global session that other streamers or user code may still need. Users
        # can call hvantk.core.context.shutdown_hail() explicitly if desired.


class StreamProcessor:
    """
    Orchestrates multiple data streamers in a pipeline.
    """

    def __init__(self, name: str):
        self.name = name
        self.streamers: List[DataStreamer] = []
        self.logger = logging.getLogger(f"{__name__}.{name}")

    def add_streamer(self, streamer: DataStreamer) -> "StreamProcessor":
        """Add a streamer to the pipeline"""
        self.streamers.append(streamer)
        # Avoid attribute errors with mocks or lightweight objects lacking a `name`
        _sname = getattr(
            streamer, "name", getattr(streamer, "__name__", streamer.__class__.__name__)
        )
        self.logger.info(f"Added streamer: {_sname}")
        return self

    def process(self, output_path: Optional[str] = None) -> Any:
        """
        Execute the streaming pipeline.

        Args:
            output_path: Optional path to save results

        Returns:
            Final processed result
        """
        self.logger.info(
            f"Starting {self.name} pipeline with {len(self.streamers)} streamers"
        )

        # Setup all streamers
        for streamer in self.streamers:
            streamer.setup()

        try:
            result = None

            # Process through each streamer in sequence
            for i, streamer in enumerate(self.streamers):
                _sname = getattr(
                    streamer,
                    "name",
                    getattr(streamer, "__name__", streamer.__class__.__name__),
                )
                self.logger.info(
                    f"Processing with streamer {i+1}/{len(self.streamers)}: {_sname}"
                )

                if i == 0:
                    # First streamer processes raw data
                    result = list(streamer.stream())
                else:
                    # Subsequent streamers process output from previous streamer
                    processed_chunks = []
                    incoming_chunks = (
                        result if isinstance(result, (list, tuple)) else [result]
                    )
                    for chunk in incoming_chunks:
                        if hasattr(streamer, "set_input"):
                            setup_ok = streamer.set_input(chunk)
                            if setup_ok:
                                # set_input succeeded; stream produces zero or more outputs
                                processed_chunks.extend(streamer.stream())
                            else:
                                # Fallback to single-chunk processing
                                processed_chunks.append(streamer.process_chunk(chunk))
                        else:
                            processed_chunks.append(streamer.process_chunk(chunk))
                    result = processed_chunks

            if output_path and result:
                self._save_result(result, output_path)

            return result

        finally:
            # Teardown all streamers
            for streamer in reversed(self.streamers):
                streamer.teardown()

    def _save_result(self, result: Any, output_path: str) -> None:
        """Persist the final pipeline result to disk.

        Supported result types:
          - str -> UTF-8 text file
          - bytes / bytearray -> binary file
          - dict / list (JSON serializable) -> pretty-printed JSON file
          - hail.Table -> checkpoint (.ht) (if output_path does not end with .ht, it is used as given)
          - list of hail.Table -> union then checkpoint

        For any other type, raise NotImplementedError to force subclasses to
        implement a custom serialization strategy.
        """
        if not output_path or not isinstance(output_path, str):
            raise ValueError("output_path must be a non-empty string")

        # Ensure parent directory exists
        parent = os.path.dirname(output_path) or "."
        try:
            os.makedirs(parent, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Failed creating parent directory '{parent}': {e}")
            raise

        try:
            # Hail Table or list[Hail Table]
            if isinstance(result, hl.Table):
                self.logger.info(f"Saving Hail Table to {output_path}")
                result.checkpoint(output_path, overwrite=True)
                return
            if (
                isinstance(result, list)
                and result
                and all(isinstance(r, hl.Table) for r in result)
            ):
                self.logger.info(
                    f"Unioning {len(result)} Hail Tables and saving to {output_path}"
                )
                combined = result[0]
                for tb in result[1:]:
                    combined = combined.union(tb)
                combined.checkpoint(output_path, overwrite=True)
                return

            # Simple Python types
            if isinstance(result, str):
                self.logger.info(f"Writing text result to {output_path}")
                with open(output_path, "w", encoding="utf-8") as fh:
                    fh.write(result)
                return
            if isinstance(result, (bytes, bytearray)):
                self.logger.info(f"Writing binary result to {output_path}")
                with open(output_path, "wb") as fh:
                    fh.write(result)
                return
            if isinstance(result, (dict, list)):
                self.logger.info(f"Writing JSON result to {output_path}")
                with open(output_path, "w", encoding="utf-8") as fh:
                    json.dump(result, fh, indent=2, ensure_ascii=False)
                return

            # Unsupported type -> delegate responsibility
            msg = (
                "_save_result does not know how to persist object of type "
                f"{type(result).__name__}; subclasses must override _save_result"
            )
            self.logger.error(msg)
            raise NotImplementedError(msg)
        except Exception as e:
            self.logger.error(f"Failed saving result to {output_path}: {e}")
            raise
