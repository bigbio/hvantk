"""
Protocol definitions for hvantk builders, streamers, and downloaders.

These protocols define the contracts that different components must follow,
enabling better type checking and clearer interfaces for extensibility.
"""

from typing import Protocol, Any, Dict, Callable, Union
from pathlib import Path
import hail as hl


class Builder(Protocol):
    """
    Protocol for data builders that create Hail Tables or MatrixTables from raw inputs.

    Builders transform raw data files (VCF, TSV, BED, etc.) into Hail's native
    Table or MatrixTable format with standardized schemas and keys.

    Example:
        class ClinVarBuilder:
            def build(self, input_path: str, **params: Any) -> hl.Table:
                ht = hl.import_vcf(input_path, reference_genome=params.get('reference_genome'))
                return ht.key_by('locus', 'alleles')

            def validate_schema(self, ht: hl.Table) -> bool:
                return 'locus' in ht.row and 'alleles' in ht.row

            def get_metadata(self) -> Dict[str, Any]:
                return {'type': 'variant', 'key': ['locus', 'alleles']}
    """

    def build(self, input_path: str, **params: Any) -> Union[hl.Table, hl.MatrixTable]:
        """
        Build a Hail Table or MatrixTable from raw input.

        Parameters
        ----------
        input_path : str
            Path to the raw input file (VCF, TSV, BED, etc.)
        **params : Any
            Additional parameters specific to the builder
            (e.g., reference_genome, fields to select, filtering options)

        Returns
        -------
        Union[hl.Table, hl.MatrixTable]
            The constructed Hail Table or MatrixTable with appropriate keys

        Raises
        ------
        ValueError
            If input_path is invalid or params are incorrect
        FileNotFoundError
            If input_path does not exist
        """
        ...

    def validate_schema(self, data: Union[hl.Table, hl.MatrixTable]) -> bool:
        """
        Validate that the output schema matches expectations.

        Parameters
        ----------
        data : Union[hl.Table, hl.MatrixTable]
            The Hail Table or MatrixTable to validate

        Returns
        -------
        bool
            True if schema is valid, False otherwise
        """
        ...

    def get_metadata(self) -> Dict[str, Any]:
        """
        Return metadata about this builder.

        Returns
        -------
        Dict[str, Any]
            Metadata dictionary containing:
            - type: str - 'variant', 'gene', 'protein', 'expression', etc.
            - key: List[str] - List of key field names
            - schema_version: str - Schema version (optional)
            - description: str - Human-readable description (optional)
        """
        ...


class Downloader(Protocol):
    """
    Protocol for data downloaders that fetch external datasets.

    Downloaders handle fetching data from remote sources (HTTP, FTP, S3, etc.)
    with checksum verification and metadata tracking.

    Example:
        class UCSCDownloader:
            def download(self, dataset_id: str, output_dir: Path, **params: Any) -> Path:
                url = f"https://cells.ucsc.edu/{dataset_id}"
                output_path = output_dir / f"{dataset_id}.tsv.gz"
                # ... download logic ...
                return output_path

            def verify_checksum(self, file_path: Path, expected: str) -> bool:
                import hashlib
                with open(file_path, 'rb') as f:
                    actual = hashlib.sha256(f.read()).hexdigest()
                return actual == expected

            def get_metadata(self, dataset_id: str) -> Dict[str, Any]:
                return {'source': 'UCSC Cell Browser', 'version': '2024-01'}
    """

    def download(self, dataset_id: str, output_dir: Path, **params: Any) -> Path:
        """
        Download dataset and return path to downloaded file(s).

        Parameters
        ----------
        dataset_id : str
            Identifier for the dataset to download
        output_dir : Path
            Directory where downloaded files should be saved
        **params : Any
            Download parameters (e.g., version, format, credentials)

        Returns
        -------
        Path
            Path to the downloaded file or directory

        Raises
        ------
        ValueError
            If dataset_id is invalid
        ConnectionError
            If download fails
        """
        ...

    def verify_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """
        Verify file integrity using checksum.

        Parameters
        ----------
        file_path : Path
            Path to the file to verify
        expected_checksum : str
            Expected checksum in format "algorithm:hash" (e.g., "sha256:abc123...")

        Returns
        -------
        bool
            True if checksum matches, False otherwise

        Raises
        ------
        FileNotFoundError
            If file_path does not exist
        ValueError
            If checksum format is invalid
        """
        ...

    def get_metadata(self, dataset_id: str) -> Dict[str, Any]:
        """
        Return metadata about the dataset.

        Parameters
        ----------
        dataset_id : str
            Identifier for the dataset

        Returns
        -------
        Dict[str, Any]
            Metadata dictionary containing:
            - source: str - Data source name
            - version: str - Dataset version
            - url: str - Download URL (optional)
            - checksum: str - Expected checksum (optional)
            - license: str - Data license (optional)
            - description: str - Dataset description (optional)
        """
        ...


# Type aliases for convenience
TableOrMatrix = Union[hl.Table, hl.MatrixTable]
BuilderFunction = Callable[..., Union[hl.Table, hl.MatrixTable]]
