"""
CLI commands for expression AnnData analysis.

Grouped under ``hvantk expression``:

- ``hvantk expression describe``: inspect metadata fields in an expression .h5ad
- ``hvantk expression summarize``: collapse an expression .h5ad into a gene-level
  summary DataFrame saved as Parquet
- ``hvantk expression markers``: extract marker genes using scanpy rank_genes_groups
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group("expression", context_settings=CONTEXT_SETTINGS)
def expression_group():
    """Expression AnnData analysis commands."""
    pass


@expression_group.command("describe")
@click.option(
    "-m",
    "--matrix-table",
    "matrix_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to an expression AnnData file (.h5ad).",
)
def describe_expression_cmd(matrix_path):
    """Inspect observation metadata fields in an expression AnnData.

    Prints available grouping variables, their types, and level counts.
    Uses pre-computed column_summary when available (instant, no compute).

    \b
    Example:
      hvantk expression describe -m data/heart_sc.h5ad
    """
    from hvantk.core.anndata_utils import load_anndata
    from hvantk.utils.matrix_utils import describe_expression_ad

    adata = load_anndata(matrix_path)
    info = describe_expression_ad(adata)

    click.echo(f"\nExpression AnnData: {matrix_path}")
    click.echo(f"  Observations (cells/samples): {info['n_obs']:,}")
    click.echo(f"  Variables (genes):            {info['n_vars']:,}")
    click.echo("")
    if info["fields"]:
        click.echo("Metadata fields:")
        for f in info["fields"]:
            if f["dtype"] == "categorical":
                click.echo(f"  {f['name']:<40s}  categorical  ({f.get('n_unique', '?')} levels)")
            else:
                click.echo(
                    f"  {f['name']:<40s}  numeric      "
                    f"[{f.get('min', '?')}, {f.get('max', '?')}]"
                )
    else:
        click.echo("  (no metadata fields found)")


@expression_group.command("summarize")
@click.option(
    "-m",
    "--matrix-table",
    "matrix_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to an expression AnnData file (.h5ad).",
)
@click.option(
    "--group-by",
    multiple=True,
    required=True,
    help="Observation metadata field(s) to group by. Repeat for multi-field grouping.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter observations: FIELD=VALUE (repeatable). "
    "Example: --filter-by time_point=9wpc --filter-by region=LV",
)
@click.option(
    "--min-cells",
    type=int,
    default=50,
    show_default=True,
    help="Minimum cells per group.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for the summary (.parquet).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def summarize_expression_cmd(
    matrix_path,
    group_by,
    filter_by,
    min_cells,
    output,
    overwrite,
):
    """Collapse an expression AnnData into a gene-level summary DataFrame.

    \b
    Examples:

      # Single-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          -o data/heart_celltype_summary.parquet

      # Multi-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type --group-by region \\
          -o data/heart_celltype_region_summary.parquet

      # With pre-filtering
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          --filter-by time_point=9wpc --filter-by region=LV \\
          --min-cells 50 \\
          -o data/heart_9wpc_LV_summary.parquet
    """
    from pathlib import Path

    from hvantk.core.anndata_utils import load_anndata
    from hvantk.utils.matrix_utils import summarize_expression_ad

    # Parse filter_by from "field=value" strings
    filters = None
    if filter_by:
        filters = {}
        for item in filter_by:
            if "=" not in item:
                raise click.BadParameter(
                    f"Expected FIELD=VALUE format, got: '{item}'",
                    param_hint="--filter-by",
                )
            key, value = item.split("=", 1)
            filters[key.strip()] = value.strip()

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        raise SystemExit(1)

    adata = load_anndata(matrix_path)

    summary_df = summarize_expression_ad(
        adata,
        group_by=list(group_by),
        filter_by=filters,
        min_cells_per_group=min_cells,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_parquet(str(output_path))

    n_genes = summary_df["gene_id"].nunique() if len(summary_df) > 0 else 0
    group_labels = sorted(summary_df["group"].unique().tolist()) if len(summary_df) > 0 else []

    click.echo(f"\nSummary table written to: {output}")
    click.echo(f"  Genes:  {n_genes:,}")
    click.echo(f"  Groups: {len(group_labels)} (from {list(group_by)})")
    click.echo("")
    if group_labels:
        click.echo(f"Group labels ({len(group_labels)}):")
        for label in group_labels[:10]:
            # Show cell count for this group
            grp_rows = summary_df[summary_df["group"] == label]
            n_cells = int(grp_rows["n_cells"].iloc[0]) if len(grp_rows) > 0 else 0
            click.echo(f"  {label:<40s}  {n_cells:>6,} cells")
        if len(group_labels) > 10:
            click.echo(f"  ... and {len(group_labels) - 10} more")


@expression_group.command("markers")
@click.option(
    "-m",
    "--matrix-table",
    "matrix_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to expression AnnData (.h5ad).",
)
@click.option(
    "--group-by",
    required=True,
    help="Observation metadata field to group by for differential expression.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter observations: FIELD=VALUE (repeatable).",
)
@click.option(
    "--method",
    type=click.Choice(["wilcoxon", "t-test", "t-test_overestim_var", "logreg"]),
    default="wilcoxon",
    show_default=True,
    help="Scoring method (passed to scanpy.tl.rank_genes_groups).",
)
@click.option(
    "--top-n",
    type=int,
    default=200,
    show_default=True,
    help="Maximum markers per group.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt).",
)
@click.option(
    "--results-tsv",
    type=click.Path(),
    default=None,
    help="Save full results table to TSV.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def markers_cmd(
    matrix_path,
    group_by,
    filter_by,
    method,
    top_n,
    output,
    results_tsv,
    overwrite,
):
    """Extract top marker genes per group from an expression AnnData.

    Uses scanpy.tl.rank_genes_groups for differential expression testing.

    \b
    Examples:

      # Wilcoxon markers
      hvantk expression markers \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          --top-n 200 \\
          -o gene_sets/heart_markers.json

      # t-test with pre-filtering
      hvantk expression markers \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          --method t-test \\
          --filter-by region=LV \\
          -o gene_sets/heart_LV_markers.json
    """
    from pathlib import Path

    import scanpy as sc

    from hvantk.core.anndata_utils import load_anndata
    from hvantk.utils.matrix_utils import filter_by_metadata_ad
    from hvantk.utils.gene_sets import GeneSet, GeneSetCollection

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        raise SystemExit(1)

    adata = load_anndata(matrix_path)

    # Apply pre-filters
    if filter_by:
        filters = {}
        for item in filter_by:
            if "=" not in item:
                raise click.BadParameter(
                    f"Expected FIELD=VALUE format, got: '{item}'",
                    param_hint="--filter-by",
                )
            key, value = item.split("=", 1)
            filters[key.strip()] = value.strip()
        adata = filter_by_metadata_ad(adata, filters)

    # Run differential expression
    sc.tl.rank_genes_groups(adata, groupby=group_by, method=method, n_genes=top_n)
    result_df = sc.get.rank_genes_groups_df(adata, group=None)

    # Save full results if requested
    if results_tsv:
        tsv_path = Path(results_tsv)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(str(tsv_path), sep="\t", index=False)
        click.echo(f"Full results table: {tsv_path} ({len(result_df):,} rows)")

    # Build gene set collection from top markers per group
    gene_sets_dict = {}
    for grp, grp_df in result_df.groupby("group"):
        top_genes = grp_df.head(top_n)["names"].tolist()
        if top_genes:
            gs = GeneSet(name=str(grp), genes=set(top_genes))
            gene_sets_dict[str(grp)] = gs

    collection = GeneSetCollection(gene_sets=gene_sets_dict, background_genes=set())

    if not collection.gene_sets:
        click.echo("Error: No marker gene sets produced.", err=True)
        raise SystemExit(1)

    # Report
    click.echo(f"Extracted markers for {len(collection)} groups (method={method})")
    for gs in sorted(collection, key=lambda g: -g.n_genes):
        click.echo(f"  {gs.name}: {gs.n_genes} markers")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".gmt":
        collection.save_gmt(str(output_path))
    else:
        collection.save(str(output_path))

    click.echo(f"\nSaved to: {output_path}")
