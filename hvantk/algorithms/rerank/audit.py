# hvantk/algorithms/rerank/audit.py
from abc import ABC, abstractmethod
import pandas as pd

# The columns CaseControlArchitectureAudit needs. Single source of truth for "does
# this cohort manifest carry the case/control architecture" -- shared by
# hvantk.algorithms.rerank.catalog.registry._default_audit and
# hvantk.tools.rerank.rerank_cli so the two audit-selection paths (build_config() and
# the plain YAML CLI) can never disagree about when the audit auto-wires up.
ARCHITECTURE_AUDIT_COLUMNS = frozenset({"n_case_var", "conc", "driver_af"})


def has_architecture_columns(columns) -> bool:
    """Whether ``columns`` is a superset of ``ARCHITECTURE_AUDIT_COLUMNS``.

    A cohort (or feature table) that declares all three columns is eligible for
    :class:`CaseControlArchitectureAudit` to be wired up automatically; one that
    doesn't falls back to :class:`NoAudit`.
    """
    return ARCHITECTURE_AUDIT_COLUMNS <= set(columns)


class Audit(ABC):
    """Advisory case/control architecture/QC audit. Returns a per-gene *reason*
    string ("" = passes) that FLAGS a gene for the user to review. It never
    overrides, re-ranks, or removes a gene — flagging is advisory only."""

    @abstractmethod
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        ...


class NoAudit(Audit):
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        return pd.Series("", index=unit_table.index)


class CaseControlArchitectureAudit(Audit):
    """Flags too-few-variants / recurrent control-shared driver / QC-leak common
    variant. Advisory only. Requires case/control architecture columns
    ``n_case_var``, ``conc`` and ``driver_af`` in the audit table (supplied via
    the cohort and/or feature tables). Use :class:`NoAudit` for tables that lack
    this architecture."""

    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        df = unit_table
        has_driver = ("driver_af" in df.columns) or ("driver_af_f" in df.columns)
        # "driver_af" is handled separately above (it accepts a "driver_af_f" alias),
        # so the direct-name requirement is ARCHITECTURE_AUDIT_COLUMNS minus it --
        # deriving from the constant instead of re-listing ("n_case_var", "conc") here
        # keeps this the single place the audit's own column requirement is spelled out.
        required_direct = sorted(ARCHITECTURE_AUDIT_COLUMNS - {"driver_af"})
        missing = [c for c in required_direct if c not in df.columns]
        if missing or not has_driver:
            need = missing + ([] if has_driver else ["driver_af"])
            raise ValueError(
                "CaseControlArchitectureAudit requires case/control architecture columns "
                f"{need} in the audit table (supplied via the cohort and/or feature tables). "
                "Declare all three columns (n_case_var, conc, driver_af) across one or "
                "more cohort_axes entries on the cohort manifest -- the manifest schema "
                "allows only 'prior' and 'cohort_axes', and declaring a required column "
                "as the prior does not work either: engine.rerank() excludes the prior "
                "column from the audit merge, so it is invisible to this check either "
                "way. 'hvantk rerank' always sets Config.audit explicitly, based on "
                "hvantk.algorithms.rerank.audit.has_architecture_columns(cohort.axis_columns()) "
                "-- there is no 'audit:' config key -- so a CLI user who does not want "
                "this audit should simply not declare all three columns, which falls "
                "back to NoAudit automatically. A direct Python-API caller may instead "
                "leave Config.audit unset (default NoAudit) or pass NoAudit() explicitly."
            )
        nv = pd.to_numeric(df["n_case_var"], errors="coerce")
        conc = pd.to_numeric(df["conc"], errors="coerce")
        daf = pd.to_numeric(
            df.get("driver_af", df.get("driver_af_f")), errors="coerce"
        ).fillna(0)
        reason = pd.Series("", index=df.index)
        # Precedence (last-write-wins): insufficient_data < recurrent_variant < common_driver.
        # Flagged set (reason != "") == old veto union (nv<=2)|(conc>=.6&af>5e-5)|(af>1e-3).
        reason = reason.mask(nv <= 2, "insufficient_data")
        reason = reason.mask((conc >= 0.6) & (daf > 5e-5), "recurrent_variant")
        reason = reason.mask(daf > 1e-3, "common_driver")
        return reason.fillna("")
