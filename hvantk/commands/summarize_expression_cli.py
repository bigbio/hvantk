"""
CLI commands for expression MatrixTable analysis.

Grouped under ``hvantk expression``:

- ``hvantk expression describe``: inspect metadata fields in an expression MT
- ``hvantk expression summarize``: collapse an expression MT into a gene-level
  summary Table grouped by metadata fields
- ``hvantk expression markers``: extract marker genes from a summary Table
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group("expression", context_settings=CONTEXT_SETTINGS)
def expression_group():
    """Expression MatrixTable analysis commands."""
    pass


@expression_group.command("describe")
@click.option(
    "-m",
    "--matrix-table",
    type=click.Path(),
    required=True,
    help="Path to an expression MatrixTable (.mt).",
)
def describe_expression_cmd(matrix_table):
    """Inspect column metadata fields in an expression MatrixTable.

    Prints available grouping variables, their types, and level counts.
    Uses pre-computed column_summary when available (instant, no Spark job).

    \b
    Example:
      hvantk expression describe -m data/heart_sc.mt
    """
    from hvantk.utils.matrix_utils import describe_expression_mt

    describe_expression_mt(matrix_table)


@expression_group.command("summarize")
@click.option(
    "-m",
    "--matrix-table",
    type=click.Path(),
    required=True,
    help="Path to an expression MatrixTable (.mt).",
)
@click.option(
    "--group-by",
    multiple=True,
    required=True,
    help="Column metadata field(s) to group by. Repeat for multi-field grouping.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter columns: FIELD=VALUE (repeatable). "
    "Example: --filter-by time_point=9wpc --filter-by region=LV",
)
@click.option(
    "--expr-field",
    default="x",
    show_default=True,
    help="Entry field containing expression values.",
)
@click.option(
    "--gene-id-field",
    default="GeneID",
    show_default=True,
    help="Row field for gene IDs.",
)
@click.option(
    "--gene-name-field",
    default="Gene Name",
    show_default=True,
    help="Row field for gene names (use '' to omit).",
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
    help="Output path for the summary Hail Table (.ht).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def summarize_expression_cmd(
    matrix_table,
    group_by,
    filter_by,
    expr_field,
    gene_id_field,
    gene_name_field,
    min_cells,
    output,
    overwrite,
):
    """Collapse an expression MatrixTable into a gene-level summary Table.

    \b
    Examples:

      # Single-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type \\
          -o data/heart_celltype_summary.ht

      # Multi-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type --group-by region \\
          -o data/heart_celltype_region_summary.ht

      # With pre-filtering
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type \\
          --filter-by time_point=9wpc --filter-by region=LV \\
          --min-cells 50 \\
          -o data/heart_9wpc_LV_summary.ht
    """
    import hail as hl
    from hvantk.utils.matrix_utils import summarize_expression

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

    gene_name = gene_name_field if gene_name_field else None

    mt = hl.read_matrix_table(matrix_table)

    tb = summarize_expression(
        mt,
        group_by=list(group_by),
        filter_by=filters,
        expr_field=expr_field,
        gene_id_field=gene_id_field,
        gene_name_field=gene_name,
        min_cells_per_group=min_cells,
        output_path=output,
        overwrite=overwrite,
    )

    n_genes = tb.count()
    # Collect group labels and sample sizes from one row
    sample_stats = hl.eval(tb.take(1)[0].stats) if n_genes > 0 else {}
    group_labels = sorted(sample_stats.keys())

    click.echo(f"\nSummary table written to: {output}")
    click.echo(f"  Genes:  {n_genes:,}")
    click.echo(f"  Groups: {len(group_labels)} (from {list(group_by)})")
    click.echo("")
    click.echo("Schema:")
    click.echo("  Key:    gene_id (str)")
    click.echo("  Fields: stats -> dict<group_label, struct{mean, fraction_expressed, n_cells}>")
    click.echo("")
    if group_labels:
        click.echo(f"Group labels ({len(group_labels)}):")
        # Show first few with cell counts
        for label in group_labels[:10]:
            info = sample_stats[label]
            click.echo(f"  {label:<40s}  {info.n_cells:>6,} cells")
        if len(group_labels) > 10:
            click.echo(f"  ... and {len(group_labels) - 10} more")


@expression_group.command("markers")
@click.option(
    "-s",
    "--summary",
    type=click.Path(),
    required=True,
    help="Path to an expression summary Hail Table (.ht) from 'expression summarize'.",
)
@click.option(
    "--method",
    type=click.Choice(["fold_change", "specificity"]),
    default="fold_change",
    show_default=True,
    help="Marker scoring method.",
)
@click.option(
    "--top-n",
    type=int,
    default=200,
    show_default=True,
    help="Maximum markers per group.",
)
@click.option(
    "--min-fold-change",
    type=float,
    default=1.5,
    show_default=True,
    help="Minimum fold change to qualify as marker.",
)
@click.option(
    "--min-fraction-expressed",
    type=float,
    default=0.1,
    show_default=True,
    help="Minimum fraction of cells expressing a gene in the group.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def markers_cmd(
    summary,
    method,
    top_n,
    min_fold_change,
    min_fraction_expressed,
    output,
    overwrite,
):
    """Extract top marker genes per group from an expression summary Table.

    \b
    Examples:

      hvantk expression markers \\
          -s data/heart_celltype_summary.ht \\
          --method fold_change \\
          --top-n 200 \\
          -o gene_sets/heart_cell_types.json

      # GMT output for GSEA compatibility
      hvantk expression markers \\
          -s data/heart_celltype_summary.ht \\
          --top-n 200 \\
          -o gene_sets/heart_cell_types.gmt
    """
    from pathlib import Path

    from hvantk.utils.gene_sets import extract_marker_gene_sets

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        raise SystemExit(1)

    collection = extract_marker_gene_sets(
        summary=summary,
        n_markers=top_n,
        min_fold_change=min_fold_change,
        min_fraction_expressed=min_fraction_expressed,
        method=method,
    )

    if not collection.gene_sets:
        click.echo("Error: No marker gene sets produced.", err=True)
        raise SystemExit(1)

    # Report
    click.echo(
        f"Extracted markers for {len(collection)} groups (method={method})"
    )
    for gs in sorted(collection, key=lambda g: -g.n_genes):
        click.echo(f"  {gs.name}: {gs.n_genes} markers")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".gmt":
        collection.save_gmt(str(output_path))
    else:
        collection.save(str(output_path))

    click.echo(f"\nSaved to: {output_path}")
