"""Load a cohort's label gene set.

Labels are a pointer to a ``GeneSetCollection`` JSON -- the format
``hvantk genesets prepare`` and the clingen/gencc/cosmic extractors already emit --
never a column in the cohort table. A collection also records where it came from
(``metadata.created_by`` / ``input_file`` / ``hgnc_validated`` / ``aliases_resolved``),
which is exactly the provenance an informally-curated clinical panel most needs.

Note this uses ``hvantk.core.utils.gene_sets.GeneSet`` (a plain named set inside a
collection), NOT ``hvantk.core.models.GeneSet`` (the typed artifact). The two are
unrelated classes that share a name.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def load_label_genes(labels) -> set[str]:
    """Return the label-positive gene identifiers named by ``labels``.

    Parameters
    ----------
    labels : hvantk.algorithms.cohort.spec.CohortLabels

    Returns
    -------
    set[str]
        The selected set's members, as written in the collection. These are
        source identifiers (``genesets prepare`` guarantees HGNC-style symbols);
        mapping them onto the spine is the caller's job, with its own rate report.

    Raises
    ------
    ValueError
        If the file holds no gene sets, if ``set_name`` is omitted for a collection
        with more than one set, if ``set_name`` names a set that is not present, or
        if the selected set is empty.
    """
    from hvantk.core.utils.gene_sets import GeneSetCollection

    collection = GeneSetCollection.load(labels.gene_set)
    names = list(collection.names())

    # GeneSetCollection.from_dict reads every field with .get(), so a JSON document
    # of the wrong shape loads as an empty collection instead of raising. Silently
    # empty labels are far worse than a hard failure.
    if not names:
        raise ValueError(
            f"{labels.gene_set} contains no gene sets; expected a "
            "GeneSetCollection JSON (as written by `hvantk genesets prepare`)"
        )

    if labels.set_name is None:
        if len(names) != 1:
            raise ValueError(
                f"{labels.gene_set} holds {len(names)} gene sets "
                f"({', '.join(sorted(names))}); set labels.set_name to choose one"
            )
        selected = names[0]
    else:
        if labels.set_name not in names:
            raise ValueError(
                f"gene set {labels.set_name!r} is not in {labels.gene_set}; "
                f"available: {', '.join(sorted(names))}"
            )
        selected = labels.set_name

    genes = set(collection.get(selected).genes)
    if not genes:
        raise ValueError(f"gene set {selected!r} in {labels.gene_set} is empty")

    logger.info("labels: %d genes from set %r", len(genes), selected)
    return genes
