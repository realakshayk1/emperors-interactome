"""covariates.py — per-pair selection covariates (PLAN_V3 T1.1).

The covariate vector used to attribute the calibration->wild shift and to fit
density ratios / conditional conformal. LOCALLY AVAILABLE covariates only:

  * score      — AF-M combined confidence (the 1-D axis WCS used).
  * iptm_mean  — interface pTM (the interface-confidence component of `score`).
  * ptm_mean   — global pTM (whole-model confidence; largely length/fold-driven).
  * iptm_ptm_gap = iptm_mean - ptm_mean — interface-vs-global confidence contrast.
  * degree     — UNION-graph degree of the pair (max of the two endpoints'
                 partner counts over ALL pairs, candidates AND decoys): a HUB /
                 selection proxy. Computed on the union (not the candidate graph)
                 so it is NOT definitionally tied to pool membership — a decoy hub
                 and a candidate hub are scored on the same footing. The wild "True"
                 candidates were pre-screened by the AP-MS+AF pipeline, so
                 high-degree hubs are still over-represented among candidates vs
                 random decoys — a covariate the 1-D `score` ratio never sees.

Covariates that need EXTERNAL data (disorder fraction via IUPred/AF pLDDT, MSA
depth/Neff, sequence length) are NOT included here; see `covariates_external()`
which requires a network grant and is a documented extension, not a dependency of
the local attribution.

Purity: these are all AF/structural/graph quantities — NO DepMap. Firewalled.
"""
from __future__ import annotations
from collections import Counter

import numpy as np
import pandas as pd

from . import config as C

LOCAL_COVARIATES = ["score", "iptm_mean", "ptm_mean", "iptm_ptm_gap", "degree"]


def _union_degree(inter: pd.DataFrame) -> dict:
    """Degree over the UNION graph (all pairs, candidates and decoys) so the
    covariate is not definitionally tied to pool membership."""
    deg = Counter()
    for a, b in zip(inter["uniprot_a"], inter["uniprot_b"]):
        deg[a] += 1
        deg[b] += 1
    return dict(deg)


def add_covariates(inter: pd.DataFrame) -> pd.DataFrame:
    """Attach the local selection-covariate columns to an interactome-shaped frame.

    `inter` must carry score, iptm_mean, ptm_mean, uniprot_a, uniprot_b.
    """
    out = inter.copy()
    out["iptm_ptm_gap"] = out["iptm_mean"].astype(float) - out["ptm_mean"].astype(float)
    deg = _union_degree(inter)
    out["degree"] = [max(deg.get(a, 0), deg.get(b, 0))
                     for a, b in zip(out["uniprot_a"], out["uniprot_b"])]
    return out


def covariate_matrix(df: pd.DataFrame, cols=None) -> np.ndarray:
    """Return a standardized covariate matrix for the given rows."""
    cols = cols or LOCAL_COVARIATES
    X = df[cols].to_numpy(dtype=float)
    mu = X.mean(0, keepdims=True)
    sd = X.std(0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def build_covariate_frame():
    """Interactome + decoys with covariates, tagged pool in {candidate, decoy}.

    Decoys are the interactome's own is_true_pair=False (random) rows, so their
    iptm_mean/ptm_mean come from the SAME table (no join needed).
    """
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    inter = add_covariates(inter)
    inter["pool"] = np.where(inter["is_true_pair"], "candidate", "decoy")
    return inter
