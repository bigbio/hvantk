"""Optional dependencies of the PTM module, each with the extra that installs it.

Mirrors ``require_scanpy`` in ``hvantk/algorithms/expression/matrix_utils.py``: the
documented standard (docs_site/getting-started/installation.md) is that a path behind an
extra exits with an actionable message naming that extra, not a traceback. ``statsmodels``
was declared in the ``constraint`` extra by #318 but ``lmm.py`` imported it bare at module
scope, so ``hvantk ptm test --test lmm`` on an ``hvantk[ptm]`` install printed
``No module named 'statsmodels'`` and never said which extra fixes it (#362).
"""

from __future__ import annotations

_STATSMODELS_HINT = (
    "statsmodels is required for the PTM constraint tests (`hvantk ptm test`), and is "
    "not part of the base install. Install the 'constraint' extra:\n"
    "    pip install 'hvantk[constraint]'\n"
    "    poetry install --extras constraint"
)


def require_statsmodels():
    """Return ``statsmodels.formula.api``, or raise ``ImportError`` naming the extra."""
    try:
        import statsmodels.formula.api as smf
    except ModuleNotFoundError as exc:
        raise ImportError(_STATSMODELS_HINT) from exc
    return smf
