# Data Streamer Base Classes
# Provides a streaming architecture for processing large genomic datasets

import logging
from abc import ABC, abstractmethod
from typing import Any, Iterator, Optional, List
import hail as hl

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
        self._hail_initialized = False

    def setup(self) -> None:
        """Initialize Hail if needed"""
        super().setup()
        if self.init_hail and not self._hail_initialized:
            hl.init()
            self._hail_initialized = True
            self.logger.info("Hail initialized")

    def teardown(self) -> None:
        """Stop Hail if we initialized it"""
        super().teardown()
        if self._hail_initialized:
            hl.stop()
            self.logger.info("Hail stopped")


class StreamProcessor:
    """
    Orchestrates multiple data streamers in a pipeline.
    """

    def __init__(self, name: str):
        self.name = name
        self.streamers: List[DataStreamer] = []
        self.logger = logging.getLogger(f"{__name__}.{name}")

    def add_streamer(self, streamer: DataStreamer) -> 'StreamProcessor':
        """Add a streamer to the pipeline"""
        self.streamers.append(streamer)
        # Avoid attribute errors with mocks or lightweight objects lacking a `name`
        _sname = getattr(streamer, "name", getattr(streamer, "__name__", streamer.__class__.__name__))
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
        self.logger.info(f"Starting {self.name} pipeline with {len(self.streamers)} streamers")

        # Setup all streamers
        for streamer in self.streamers:
            streamer.setup()

        try:
            result = None

            # Process through each streamer in sequence
            for i, streamer in enumerate(self.streamers):
                _sname = getattr(streamer, "name", getattr(streamer, "__name__", streamer.__class__.__name__))
                self.logger.info(f"Processing with streamer {i+1}/{len(self.streamers)}: {_sname}")

                if i == 0:
                    # First streamer processes raw data
                    result = list(streamer.stream())
                else:
                    # Subsequent streamers process output from previous streamer
                    processed_chunks = []
                    for chunk in result if isinstance(result, (list, tuple)) else [result]:
                        processed_chunks.extend(streamer.stream() if hasattr(streamer, 'set_input')
                                              and streamer.set_input(chunk) else [streamer.process_chunk(chunk)])
                    result = processed_chunks

            if output_path and result:
                self._save_result(result, output_path)

            return result

        finally:
            # Teardown all streamers
            for streamer in reversed(self.streamers):
                streamer.teardown()

    def _save_result(self, result: Any, output_path: str) -> None:
        """Save the final result to the specified path"""
        self.logger.info(f"Saving results to {output_path}")
        # Implementation depends on result type - will be handled by specific streamers
