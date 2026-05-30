# Generic training-set builder.
#
# Composes a variant data source (a VariantTableStreamer) with a sequence of
# algorithm-level Annotators and an outcome specification (label or score) into
# a single Hail Table. Lives in tools/ because only tools/ may compose a
# skills/ data source with algorithms/ annotators.
#
# This module is deliberately source-agnostic: it imports ONLY from core/ and
# algorithms/, never from skills/. Source-specific wiring (e.g. ClinVar label
# derivation and default annotators) lives in sibling modules such as
# tools/training_sets/clinvar.py.

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional

import hail as hl

from hvantk.algorithms.annotation.annotator import Annotator
from hvantk.core.streamers.variant_table import VariantTableStreamer

logger = logging.getLogger(__name__)


@dataclass
class OutcomeSpec:
    """Describes the supervised target column of a training set.

    Attributes:
        column: Name of the outcome column in the final table.
        kind: ``"label"`` (categorical target, e.g. TP/TN) or ``"score"``
            (continuous target). This is descriptive metadata (label vs
            continuous score) recorded for logging/provenance; it does not
            branch the logic in :meth:`TrainingSetBuilder.build`.
        derive: Optional whole-table transform that adds ``column`` to the
            input table (e.g. ClinVar TP/TN derivation). If ``None`` the column
            is expected to already exist on the base table.
    """

    column: str
    kind: Literal["label", "score"]
    derive: Optional[Callable[[hl.Table], hl.Table]] = None


class TrainingSetBuilder:
    """Build a training set from a base variant table + annotators + outcome.

    The builder is parameterized over a :class:`VariantTableStreamer` (the data
    source), a list of :class:`Annotator` steps, and an :class:`OutcomeSpec`.
    It performs no gene-filtering: any gating on a gene set belongs in the
    outcome's ``derive`` function (matching the legacy ClinVar semantics).
    """

    def __init__(
        self,
        base: VariantTableStreamer,
        annotators: List[Annotator],
        outcome: OutcomeSpec,
        *,
        post_process: Optional[Callable[[hl.Table], hl.Table]] = None,
        name: str = "TrainingSetBuilder",
    ) -> None:
        self.base = base
        self.annotators = annotators
        self.outcome = outcome
        self.post_process = post_process
        self.name = name
        self.logger = logging.getLogger(f"{__name__}.{name}")

    def build(
        self, output_path: Optional[str] = None, export_tsv: bool = False
    ) -> hl.Table:
        """Run the full build and return the annotated Hail Table."""
        self.logger.info(
            "Building training set '%s': outcome column '%s' (kind=%s), "
            "%d annotator(s)",
            self.name,
            self.outcome.column,
            self.outcome.kind,
            len(self.annotators),
        )

        # 1. Base variant table (no gene-filtering in the generic builder).
        ht = self.base.to_hail()

        # 2. Derive the outcome column if a derivation was supplied.
        if self.outcome.derive is not None:
            ht = self.outcome.derive(ht)

        # 3. Require the outcome column and drop rows missing it.
        if self.outcome.column not in ht.row:
            raise ValueError(
                f"Outcome column '{self.outcome.column}' not present in the "
                f"training table after derivation. Available row fields: "
                f"{list(ht.row)}"
            )
        ht = ht.filter(hl.is_defined(ht[self.outcome.column]))

        # 4. Apply each annotator in sequence.
        for annotator in self.annotators:
            annotator.setup()
            try:
                ht = annotator.process_chunk(ht)
            finally:
                annotator.teardown()

        # 5. Optional post-processing (e.g. composite features).
        if self.post_process is not None:
            ht = self.post_process(ht)

        # 6. Optionally checkpoint + export.
        if output_path:
            self.logger.info("Writing training set to %s", output_path)
            ht = ht.checkpoint(output_path, overwrite=True)
            if export_tsv:
                ht.export(output_path + ".tsv")

        # 7. Return the lazy Hail Table.
        return ht
