# hvantk/algorithms/rerank/veto.py
from abc import ABC, abstractmethod
import pandas as pd


class Veto(ABC):
    @abstractmethod
    def apply(self, unit_table: pd.DataFrame) -> pd.Series: ...


class NoOpVeto(Veto):
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        return pd.Series(False, index=unit_table.index)


class CaseControlArchitectureVeto(Veto):
    """Ported from chd_veto.py: recurrent control-shared driver / too-few variants / QC-leak common variant.

    Requires case/control architecture columns ``n_case_var``, ``conc`` and ``driver_af``
    in the veto table (supplied via the cohort and/or feature tables). Use
    :class:`NoOpVeto` for tables that lack this architecture.
    """
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        df = unit_table
        has_driver = ("driver_af" in df.columns) or ("driver_af_f" in df.columns)
        missing = [c for c in ("n_case_var", "conc") if c not in df.columns]
        if missing or not has_driver:
            need = missing + ([] if has_driver else ["driver_af"])
            raise ValueError(
                "CaseControlArchitectureVeto requires case/control architecture columns "
                f"{need} in the veto table (supplied via the cohort and/or feature tables). "
                "Use NoOpVeto (omit the cohort) for tables without this architecture.")
        nv = pd.to_numeric(df["n_case_var"], errors="coerce")
        conc = pd.to_numeric(df["conc"], errors="coerce")
        daf = pd.to_numeric(df.get("driver_af", df.get("driver_af_f")), errors="coerce").fillna(0)
        arch_fragile = (nv <= 2) | ((conc >= 0.6) & (daf > 5e-5))
        qc_flag = daf > 1e-3
        return (arch_fragile | qc_flag).fillna(False)
