"""hvantk.data package.

Data management utilities including:
- Base streamer classes (DataStreamer, HailDataStreamer, StreamProcessor)
- Gene-disease validity base (GeneDiseaseValidityStreamer)
- Data source streamers (ClinvarDataStreamer, ClinGenStreamer, GenCCStreamer)
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
from hvantk.data.gene_disease_streamer import GeneDiseaseValidityStreamer
from hvantk.data.gene_mapper import GeneMapper
from hvantk.data.gencc_streamer import GenCCStreamer

__all__ = [
    # Base classes
    "DataStreamer",
    "HailDataStreamer",
    "StreamProcessor",
    "GeneDiseaseValidityStreamer",
    # Data source streamers
    "ClinGenStreamer",
    "GenCCStreamer",
    "ClinvarDataStreamer",
    "ClinvarTrainingSetProcessor",
    "create_clinvar_training_set_streamer",
    # Gene ID mapping
    "GeneMapper",
]
