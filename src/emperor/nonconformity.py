"""nonconformity.py — scalar nonconformity score s_i (METHODS §2).

Lower s_i = more confident it's a true interaction. Table 5 exposes only the
AF-Multimer combined confidence `score` (ipTM axis; no pDockQ2/interface-PAE),
so the base score is:

    s_i = (1 - score_i) + w_phys * phys_penalty_i

with phys_penalty in [0,1] (0 until the structure-derived term is wired in via
the BioNeMo step; see METHODS §4). This keeps s_i purely structural/physical —
DepMap NEVER enters here (held-out firewall).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from . import config as C


def nonconformity(df: pd.DataFrame, w_phys: float = 1.0) -> np.ndarray:
    """Compute s_i from `score` (+ optional `phys_penalty` column)."""
    s = 1.0 - df["score"].to_numpy(dtype=float)
    if "phys_penalty" in df.columns:
        s = s + w_phys * df["phys_penalty"].to_numpy(dtype=float)
    return s
