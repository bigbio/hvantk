# local/rerank_engine/catalog/profile.py
from dataclasses import dataclass, field
from typing import Optional, Callable
from hvantk.algorithms.rerank.config import PriorSpec, CohortSpec

@dataclass
class DiseaseProfile:
    name: str
    # inputs / universe
    prior: Optional[PriorSpec] = None
    cohort: Optional[CohortSpec] = None
    # precomputed feature source: callable -> (matrix[gene+feature cols], {family: [cols]})
    precomputed: Optional[Callable] = None
    # tissue / cell-type context (drives expression/eqtl/ptm axes)
    tissue: Optional[str] = None
    cell_types: list = field(default_factory=list)
    dev_window: Optional[str] = None
    # labels
    disease_terms: list = field(default_factory=list)
    label_classifications: list = field(default_factory=lambda: ["Definitive", "Strong", "Moderate"])
    extra_positive_genes: set = field(default_factory=set)
    # per-cohort PTM feature parquet (built by chd_ptm_features.py prep step)
    ptm_features_path: Optional[str] = None
    # knobs
    min_label_coverage: float = 0.5
    extra_flagged_genes: list = field(default_factory=list)
