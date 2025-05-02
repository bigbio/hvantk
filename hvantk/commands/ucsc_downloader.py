import click
import os

from hvantk.datasets.ucsc_cell_datasets import UCSCDataSetCollection
from hvantk.settings import CONTEXT_SETTINGS
from hvantk.utils.file_utils import download_file, url_exists
from hvantk.utils.constants import (
    UCSC_CELL_BROWSER_BASE_URL,
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
    UCSC_JSON_FILE_PATH,
)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """A package for gene and variant annotation."""
    pass


def _print_dataset_names():
    """
    Prints the names of all available UCSC datasets in a formatted list.

    :raises FileNotFoundError: If the JSON file containing UCSC dataset
                               information is not found at the specified
                               path.
    :raises JSONDecodeError: If the JSON file contains invalid JSON format.
    """
    try:
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        dataset_names = collection.list_dataset_names()
        if not dataset_names:
            click.echo("No datasets available.")
        else:
            click.echo("Available datasets:")
            for name in dataset_names:
                click.echo(f"- {name}")
    except ValueError as e:
        click.echo(f"Error: {e}")


@click.command("ucsc-downloader", short_help="Download UCSC Cell Browser data")
@click.option(
    "--dataset",
    default="adultPancreas",
    type=str,
    required=True,
    help="Dataset to download from UCSC Cell Browser.",
)
@click.option(
    "--output-dir",
    default="data/ucsc",
    type=str,
    required=False,
    help="Output directory for downloaded files.",
)
@click.option(
    "--base_url",
    default=UCSC_CELL_BROWSER_BASE_URL,
    type=str,
    required=False,
    help="URL to download the data from.",
)
@click.option(
    "--list_datasets",
    is_flag=True,
    required=False,
    help="List all available datasets in the UCSC collection.",
)
@click.pass_context
def ucsc_downloader(ctx, dataset, output_dir, base_url, list_datasets):
    """
    Download expression matrix and metadata from the UCSC Cell Browser for a specified dataset (e.g., adultPancreas).

    Usage examples:
        # Download a specific dataset
        hvantk ucsc-downloader --dataset <dataset_name> --output-dir <output_directory>

        # Download a specific dataset with a custom base URL
        hvantk ucsc-downloader --dataset <dataset_name> --output-dir <output_directory> --base_url <custom_base_url>

        # List all available datasets
        hvantk ucsc-downloader --list_datasets
    """

    if list_datasets:
        _print_dataset_names()
        return

    # Construct the URLs for the expression matrix and metadata files
    url_expression_matrix = f"{base_url}/{dataset}/{EXPRESSION_MATRIX_FILE_NAME}"
    url_metadata = f"{base_url}/{dataset}/{METADATA_FILE_NAME}"

    # check that both urls exist
    if not url_exists(url_expression_matrix):
        click.echo(f"Error: Expression matrix URL does not exist: {url_expression_matrix}")
        return
    if not url_exists(url_metadata):
        click.echo(f"Error: Metadata URL does not exist: {url_metadata}")
        return

    # Validate output directory
    if os.path.exists(output_dir) and not os.path.isdir(output_dir):
        click.echo(f"Error: Output path exists but is not a directory: {output_dir}")
        return

    # make dataset-specific directory
    dataset_dir = f"{output_dir}/{dataset}"

    # Check if dataset directory already exists and has content
    if os.path.exists(dataset_dir) and os.listdir(dataset_dir):
        if not click.confirm(f"Directory {dataset_dir} already exists and has content. Overwrite?"):
            click.echo("Download canceled.")
            return

    # Download the expression matrix file
    download_file(url_expression_matrix, dataset_dir, EXPRESSION_MATRIX_FILE_NAME)
    # Download the metadata file
    download_file(url_metadata, dataset_dir, METADATA_FILE_NAME)

    click.echo(f"Data downloaded to {dataset_dir}")

if __name__ == "__main__":
    cli()
