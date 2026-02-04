"""hvantk.data package.

Data management utilities including:
- Base streamer classes (DataStreamer, HailDataStreamer, StreamProcessor)
- Data source streamers (ClinvarDataStreamer, ClinGenStreamer)
- Gene ID mapping utility (GeneMapper)
"""

from hvantk.data.clingen_streamer import ClinGenStreamer
from hvantk.data.clinvar_streamer import (
    ClinvarDataStreamer,
    ClinvarTrainingSetProcessor,
    create_clinvar_training_set_streamer,
)
from hvantk.data.data_streamer import (
    DataStreamer,
    HailDataStreamer,
    StreamProcessor,
)
from hvantk.data.gene_mapper import GeneMapper

__all__ = [
    # Base classes
    "DataStreamer",
    "HailDataStreamer",
    "StreamProcessor",
    # Data source streamers
    "ClinGenStreamer",
    "ClinvarDataStreamer",
    "ClinvarTrainingSetProcessor",
    "create_clinvar_training_set_streamer",
    # Gene ID mapping
    "GeneMapper",
]
