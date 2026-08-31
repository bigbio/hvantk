"""The ``Artifact`` union — a name for "one of the four typed artifacts".

``Artifact`` used to be a ``runtime_checkable`` Protocol here. Nothing ever
imported it and nothing dispatched on it, so it was removed as dead code. But
the *word* stayed load-bearing: docstrings, error messages and the plugin
contract all say "returns an Artifact", and out-of-tree builders reasonably
write ``-> Artifact``. This module gives that vocabulary a real symbol again,
as a union of the concrete types rather than a structural Protocol.

Why a union and not a Protocol: the four artifacts are a closed set defined in
this package, so enumerating them is both honest and more useful. ``isinstance``
works against a union on Python 3.10+, and a type checker can narrow it, which a
``runtime_checkable`` Protocol cannot do (it only checks attribute presence).

**This module is deliberately not imported by** ``hvantk.core.models.__init__``
**at import time.** Building the union requires all four classes, which pulls in
pandas (``annotation_table``) and anndata (``expression_matrix``) -- about 1.2 s.
The package resolves every public name lazily via PEP 562 precisely to avoid
that, so ``Artifact`` is registered in ``_EXPORTS`` and costs nothing until
someone actually asks for it. Annotating with it is free under
``from __future__ import annotations``, which makes annotations strings.

A plugin manifest's ``artifact_type:`` must still name a *concrete* type -- the
union would make the build-time isinstance check vacuous, so the loader rejects
it.
"""
from __future__ import annotations

from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.gene_set import GeneSet
from hvantk.core.models.variant_matrix import VariantMatrix

#: Any value a plugin builder may return. Kept in sync with the concrete
#: artifact classes by ``test_artifact_union_covers_every_artifact_class``.
Artifact = AnnotationTable | ExpressionMatrix | VariantMatrix | GeneSet

__all__ = ["Artifact"]
