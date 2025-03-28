from typing import List, Optional
from dataclasses import dataclass, field
import json

from hvantk.utils.constants import (UCSC_CELL_BROWSER_BASE_URL,
                                    EXPRESSION_MATRIX_FILE_NAME,
                                    METADATA_FILE_NAME)
from hvantk.utils.file_utils import download_file

@dataclass
class DatasetFacets:
    body_parts: List[str]
    organisms: List[str]
    projects: List[str]
    diseases: List[str]
    life_stages: List[str]
    domains: List[str]
    assays: List[str]
    sources: List[str]

@dataclass
class UCSCDataset:
    shortLabel: str
    name: str
    md5: str
    hasFiles: Optional[List[str]] = field(default_factory=list)
    body_parts: Optional[List[str]] = field(default_factory=list)
    organisms: Optional[List[str]] = field(default_factory=list)
    tags: Optional[List[str]] = field(default_factory=list)
    projects: Optional[List[str]] = field(default_factory=list)
    diseases: Optional[List[str]] = field(default_factory=list)
    life_stages: Optional[List[str]] = field(default_factory=list)
    domains: Optional[List[str]] = field(default_factory=list)
    sources: Optional[List[str]] = field(default_factory=list)
    assays: Optional[List[str]] = field(default_factory=list)
    facets: Optional[DatasetFacets] = None
    sampleCount: Optional[int] = None
    isCollection: Optional[bool] = False
    collectionCount: Optional[int] = None
    datasetCount: Optional[int] = None

    def __post_init__(self):
        if self.facets is None:
            self.facets = DatasetFacets(
                body_parts=self.body_parts,
                organisms=self.organisms,
                projects=self.projects,
                diseases=self.diseases,
                life_stages=self.life_stages,
                domains=self.domains,
                assays=self.assays,
                sources=self.sources
            )

    def summary(self) -> str:
        return (f"Dataset {self.name} ({self.shortLabel}): \n"
                f"{len(self.body_parts)} body parts, \n"
                f"{len(self.organisms)} organisms, \n"
                f"{self.sampleCount} samples")

    def download_expression_matrix(self, out_dir: str):
        url_download = f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/{EXPRESSION_MATRIX_FILE_NAME}"
        download_file(url=url_download, out_dir=out_dir, file_name=EXPRESSION_MATRIX_FILE_NAME)


    def download_metadata(self, out_dir: str):
        url_download = f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/{METADATA_FILE_NAME}"
        download_file(url=url_download, out_dir=out_dir, file_name=METADATA_FILE_NAME)


@dataclass
class UCSCDataSetCollection:
    shortLabel: str
    abstract: str
    inDir: str
    name: str
    datasets: List[UCSCDataset]

    @classmethod
    def from_json(cls, json_path: str) -> 'UCSCDataSetCollection':
        try:
            with open(json_path, 'r') as file:
                data = json.load(file)
            
            required_keys = ['shortLabel', 'abstract', 'inDir', 'name', 'datasets']
            missing_keys = [key for key in required_keys if key not in data]
            if missing_keys:
                raise ValueError(f"Missing required keys in JSON: {', '.join(missing_keys)}")
                
            datasets = [UCSCDataset(**ds) for ds in data['datasets']]
            return cls(
                shortLabel=data['shortLabel'],
                abstract=data['abstract'],
                inDir=data['inDir'],
                name=data['name'],
                datasets=datasets
            )
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON format: {str(e)}")
        except OSError as e:
            raise ValueError(f"Could not read file {json_path}: {str(e)}")
    def get_dataset_by_name(self, dataset_name: str) -> Optional[UCSCDataset]:
        for dataset in self.datasets:
            if dataset.name == dataset_name:
                return dataset
        return None

    def total_samples(self) -> int:
        return sum(dataset.sampleCount or 0 for dataset in self.datasets)

    def list_dataset_names(self) -> List[str]:
        return [dataset.name for dataset in self.datasets]
