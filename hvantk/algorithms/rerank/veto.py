# local/rerank_engine/veto.py
from abc import ABC, abstractmethod
import numpy as np, pandas as pd

class Veto(ABC):
    @abstractmethod
    def apply(self, unit_table: pd.DataFrame) -> pd.Series: ...

class NoOpVeto(Veto):
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        return pd.Series(False, index=unit_table.index)

class CaseControlArchitectureVeto(Veto):
    """Ported from chd_veto.py: recurrent control-shared driver / too-few variants / QC-leak common variant."""
    def apply(self, unit_table: pd.DataFrame) -> pd.Series:
        df = unit_table
        nv = pd.to_numeric(df.n_case_var, errors="coerce")
        conc = pd.to_numeric(df.conc, errors="coerce")
        daf = pd.to_numeric(df.get("driver_af", df.get("driver_af_f")), errors="coerce").fillna(0)
        arch_fragile = (nv <= 2) | ((conc >= 0.6) & (daf > 5e-5))
        qc_flag = daf > 1e-3
        return (arch_fragile | qc_flag).fillna(False)
