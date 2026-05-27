"""Phase N: tests for provenance chaining through @algorithm decorator."""
from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd

from hvantk.core.models._expr import col
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.backends import Backend, algorithm
from hvantk.core.models.provenance import Provenance


def _prov(name: str = "test") -> Provenance:
    return Provenance(
        plugin=name,
        dataset=f"{name}:rows",
        plugin_version="0.0",
        source_fingerprint=f"sha256:{name}-abc",
        schema_id=f"{name}-rows-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


def _df() -> pd.DataFrame:
    return pd.DataFrame({"gene": ["BRCA1", "BRCA2"], "score": [0.7, 0.4]})


def test_single_artifact_input_chains_to_output():
    @algorithm(name="passthrough", backends=[Backend.PANDAS])
    def passthrough(ann: AnnotationTable) -> AnnotationTable:
        # Build a fresh artifact with empty provenance.parents — the decorator
        # should fill them in from the input.
        return AnnotationTable.from_pandas(
            ann.to_pandas(), provenance=_prov("derived")
        )

    input_ann = AnnotationTable.from_pandas(_df(), provenance=_prov("source"))
    output = passthrough(input_ann)

    # The output's provenance is the algorithm-supplied one PLUS chained parents.
    assert output.provenance.plugin == "derived"  # original schema_id stays
    assert len(output.provenance.parents) == 1
    assert output.provenance.parents[0].plugin == "source"


def test_two_artifact_inputs_chain_both_provenances():
    @algorithm(name="two_input", backends=[Backend.PANDAS])
    def two_input(a: AnnotationTable, b: AnnotationTable) -> AnnotationTable:
        merged = a.to_pandas()  # ignore b for test simplicity
        return AnnotationTable.from_pandas(merged, provenance=_prov("merged"))

    a = AnnotationTable.from_pandas(_df(), provenance=_prov("a"))
    b = AnnotationTable.from_pandas(_df(), provenance=_prov("b"))
    output = two_input(a, b)

    plugins = {p.plugin for p in output.provenance.parents}
    assert plugins == {"a", "b"}


def test_kwarg_artifact_inputs_are_chained():
    @algorithm(name="kwarg_input", backends=[Backend.PANDAS])
    def kwarg_input(*, ann: AnnotationTable) -> AnnotationTable:
        return AnnotationTable.from_pandas(
            ann.to_pandas(), provenance=_prov("derived")
        )

    src = AnnotationTable.from_pandas(_df(), provenance=_prov("source"))
    output = kwarg_input(ann=src)

    assert len(output.provenance.parents) == 1
    assert output.provenance.parents[0].plugin == "source"


def test_non_artifact_args_do_not_chain():
    @algorithm(name="scalar_args", backends=[Backend.PANDAS])
    def scalar_args(ann: AnnotationTable, threshold: float = 0.5) -> AnnotationTable:
        return AnnotationTable.from_pandas(
            ann.to_pandas(), provenance=_prov("derived")
        )

    src = AnnotationTable.from_pandas(_df(), provenance=_prov("source"))
    output = scalar_args(src, threshold=0.3)

    # Only the AnnotationTable is chained; the float threshold is ignored.
    assert len(output.provenance.parents) == 1
    assert output.provenance.parents[0].plugin == "source"


def test_non_artifact_output_does_not_chain():
    @algorithm(name="scalar_output", backends=[Backend.PANDAS])
    def scalar_output(ann: AnnotationTable) -> float:
        return ann.count() * 0.5

    src = AnnotationTable.from_pandas(_df(), provenance=_prov("source"))
    output = scalar_output(src)

    # Output is a float — no provenance attribute, no chaining attempted.
    assert isinstance(output, float)


def test_explicit_internal_chaining_is_preserved():
    """If an algorithm already sets parents internally, the decorator does not overwrite."""
    deliberate_parent = _prov("deliberate")

    @algorithm(name="self_chained", backends=[Backend.PANDAS])
    def self_chained(ann: AnnotationTable) -> AnnotationTable:
        import dataclasses
        prov = dataclasses.replace(
            _prov("derived"), parents=(deliberate_parent,)
        )
        return AnnotationTable.from_pandas(ann.to_pandas(), provenance=prov)

    src = AnnotationTable.from_pandas(_df(), provenance=_prov("source"))
    output = self_chained(src)

    # The internal chain (deliberate) is preserved; the decorator did NOT
    # overwrite with the input (source).
    assert len(output.provenance.parents) == 1
    assert output.provenance.parents[0].plugin == "deliberate"
