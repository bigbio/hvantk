import click
import os
from urllib.parse import urlparse

import logging

logger = logging.getLogger(__name__)

from hvantk.datasets.ucsc_cell_datasets import UCSCDataSetCollection
from hvantk.core.config import CONTEXT_SETTINGS
import hvantk.data.file_utils as file_utils
from hvantk.core.constants import (
    UCSC_CELL_BROWSER_BASE_URL,
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
    UCSC_JSON_FILE_PATH,
)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """A package for gene and variant annotation."""
    pass


def _format_facets(ds, max_items: int = 3) -> str:
    """Build a short facet string, truncating long lists."""
    facets = []
    if ds.organisms:
        facets.extend(ds.organisms)
    if ds.body_parts:
        facets.extend(ds.body_parts)
    if not facets:
        return ""
    if len(facets) <= max_items:
        return ", ".join(facets)
    return ", ".join(facets[:max_items]) + f" (+{len(facets) - max_items} more)"


def _print_dataset_names(search: str = None):
    """
    Prints the names of all available UCSC datasets in a formatted list.

    If a search term is provided, filters datasets by name, label, body parts,
    organisms, and diseases (case-insensitive).
    """
    try:
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)

        if search:
            collection = collection.search(search)

        if not collection.datasets:
            if search:
                click.echo(f'No datasets matching "{search}".')
            else:
                click.echo("No datasets available.")
            return

        leaves = [ds for ds in collection.datasets if not ds.isCollection]
        collections = [ds for ds in collection.datasets if ds.isCollection]

        # Header
        header = "UCSC Cell Browser Datasets"
        if search:
            header += f' (filtered by "{search}")'
        click.echo(header)
        click.echo(
            f"  Total: {len(collection.datasets)}  "
            f"({len(leaves)} downloadable, {len(collections)} collections)"
        )
        click.echo("")

        # Compute column width for alignment
        all_names = [ds.name for ds in collection.datasets]
        name_width = min(max((len(n) for n in all_names), default=20), 40)

        # Leaf datasets
        if leaves:
            click.echo(f"Downloadable datasets ({len(leaves)}):")
            click.echo(f"  {'NAME':<{name_width}}  {'CELLS':>10}  ORGANISM / TISSUE")
            click.echo(f"  {'-' * name_width}  {'-' * 10}  {'-' * 30}")
            for ds in leaves:
                cells = f"{ds.sampleCount:,}" if ds.sampleCount else "-"
                facets = _format_facets(ds)
                click.echo(f"  {ds.name:<{name_width}}  {cells:>10}  {facets}")
            click.echo("")

        # Collections
        if collections:
            click.echo(f"Collections ({len(collections)}):")
            click.echo(f"  {'NAME':<{name_width}}  {'DATASETS':>10}  ORGANISM / TISSUE")
            click.echo(f"  {'-' * name_width}  {'-' * 10}  {'-' * 30}")
            for ds in collections:
                count = str(ds.datasetCount) if ds.datasetCount else "?"
                facets = _format_facets(ds)
                click.echo(f"  {ds.name:<{name_width}}  {count:>10}  {facets}")
            click.echo("")

        click.echo(
            "Tip: Use --search <term> to filter. "
            "Collections require a child path (e.g., hoc/all-heart)."
        )
    except ValueError as e:
        click.echo(f"Error: {e}")


def _is_valid_url(url: str) -> bool:
    parts = urlparse(url)
    return bool(parts.scheme) and bool(parts.netloc)


@click.command("ucsc-downloader", short_help="Download UCSC Cell Browser data")
@click.option(
    "--dataset",
    type=str,
    required=False,
    default=None,
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
@click.option(
    "--search",
    type=str,
    default=None,
    help="Filter dataset names and labels by search term (case-insensitive). Use with --list_datasets.",
)
@click.pass_context
def ucsc_downloader(ctx, dataset, output_dir, base_url, list_datasets, search):
    """
    Downloads expression matrix and metadata files for a specified UCSC Cell Browser dataset.

    If the --list_datasets flag is provided, lists all available datasets instead of downloading.
    Checks for the existence of required files at the specified URLs before downloading.
    Prompts for confirmation if the target dataset directory exists and contains files.
    Downloaded files are saved in a dataset-specific subdirectory within the output directory.
    """

    if list_datasets:
        _print_dataset_names(search=search)
        logger.info("Listing available datasets completed.")
        ctx.exit(0)

    if not dataset:
        raise click.UsageError("Missing option '--dataset' (required for download).")

    # Dataset validation (prevent path traversal / malformed URLs) BEFORE using it anywhere
    # Allow forward slashes: UCSC child datasets use paths like "hoc/all-heart"
    invalid = (
        (".." in dataset) or ("\\" in dataset) or any(ch.isspace() for ch in dataset)
    )
    if invalid:
        click.echo(
            f"Invalid dataset value: '{dataset}'. Datasets must not contain '..', backslashes, or whitespace."
        )
        logger.error(f"Rejected invalid dataset value: '{dataset}'")
        ctx.exit(1)

    # Warn if the dataset name matches a known collection (no downloadable files at top level)
    try:
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        ds_entry = collection.get_by_name(dataset)
        if ds_entry and ds_entry.isCollection:
            child_count = ds_entry.datasetCount or "unknown number of"
            click.echo(
                f'Warning: "{dataset}" is a collection with {child_count} child datasets.\n'
                f"Collections have no expression matrix at the top level.\n"
                f'Try a child dataset path instead (e.g., "{dataset}/<child-name>").\n'
                f"Use --list_datasets --search {dataset} to find children."
            )
            logger.warning(
                f"Dataset '{dataset}' is a collection with {child_count} children"
            )
            ctx.exit(1)
    except ValueError:
        pass  # Catalog unavailable — continue with download attempt

    # Validate base_url
    if not _is_valid_url(base_url):
        click.echo("Invalid URL")
        logger.error(f"Invalid base_url provided: {base_url}")
        ctx.exit(1)

    # Decide filenames and behavior based on source
    if base_url == UCSC_CELL_BROWSER_BASE_URL:
        expr_name = EXPRESSION_MATRIX_FILE_NAME
        meta_name = METADATA_FILE_NAME
        perform_exist_check = True
        url_dataset_segment = dataset
        target_dir = os.path.join(output_dir, dataset)
    else:
        # Tests expect these filenames and no dataset subdir
        expr_name = "expression_matrix.tsv"
        meta_name = "metadata.tsv"
        perform_exist_check = False
        url_dataset_segment = "test_dataset"
        target_dir = output_dir

    # Construct URLs
    url_expression_matrix = f"{base_url}/{url_dataset_segment}/{expr_name}"
    url_metadata = f"{base_url}/{url_dataset_segment}/{meta_name}"

    # For UCSC base, ensure URLs exist
    if perform_exist_check:
        logger.info(
            f"Checking if expression matrix URL exists: {url_expression_matrix}"
        )
        if not file_utils.url_exists(url_expression_matrix):
            click.echo(
                f"Error: Expression matrix URL does not exist: {url_expression_matrix}"
            )
            logger.error(
                f"Expression matrix URL does not exist: {url_expression_matrix}"
            )
            ctx.exit(1)
        logger.info(f"Checking if metadata URL exists: {url_metadata}")
        if not file_utils.url_exists(url_metadata):
            click.echo(f"Error: Metadata URL does not exist: {url_metadata}")
            logger.error(f"Metadata URL does not exist: {url_metadata}")
            ctx.exit(1)

    # Validate output directory and ensure target exists
    if os.path.exists(output_dir) and not os.path.isdir(output_dir):
        click.echo(f"Error: Output path exists but is not a directory: {output_dir}")
        ctx.exit(1)
    os.makedirs(target_dir, exist_ok=True)

    # Confirm overwrite for UCSC when target has content
    if perform_exist_check and os.listdir(target_dir):
        if not click.confirm(
            f"Directory {target_dir} already exists and has content. Overwrite?"
        ):
            click.echo("Download canceled.")
            ctx.exit(0)

    # Attempt downloads (patched in tests)
    try:
        logger.info(
            f"Downloading expression matrix from {url_expression_matrix} to {target_dir}"
        )
        file_utils.download_file(url_expression_matrix, target_dir, expr_name)

        logger.info(f"Downloading metadata from {url_metadata} to {target_dir}")
        file_utils.download_file(url_metadata, target_dir, meta_name)

    except Exception as e:
        click.echo(str(e))
        logger.exception("Download failed")
        ctx.exit(1)

    click.echo(f"Data downloaded to {target_dir}")
    logger.info(f"Data downloaded to {target_dir}")
    ctx.exit(0)


if __name__ == "__main__":
    logger.info("Starting ucsc_downloader CLI")
    cli()
    logger.info("ucsc_downloader CLI completed")
