"""test_firewall.py — adversarial purity-firewall + no-signal tests (PLAN_V3 T0.3).

The single load-bearing project rule: DepMap co-essentiality is the held-out
REFEREE and must NEVER enter the nonconformity score, the calibration labels, or
the certification. These tests fail loudly if that firewall leaks, and a
label-shuffle negative control asserts the certified-vs-dropped enrichment
collapses to null when the structure signal is destroyed (catches a pipeline that
is leaking referee structure into selection).
"""
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from emperor import config as C
from emperor import nonconformity as nc
from emperor.conformal import conformal_pvalues, benjamini_hochberg
from emperor.nonconformity import nonconformity

ROOT = Path(C.ROOT)
DEPMAP_TOKENS = ("depmap", "gls", "coess", "co_ess", "co-ess", "essential")


# ---------- static firewall: score/label code never references DepMap ----------
@pytest.mark.parametrize("mod", [nc])
def test_nonconformity_source_has_no_depmap_token(mod):
    src = inspect.getsource(mod).lower()
    # allow the word inside the module docstring's firewall NOTE, but not in code:
    code_lines = [ln for ln in src.splitlines()
                  if not ln.strip().startswith("#") and '"""' not in ln]
    code = "\n".join(code_lines)
    # strip the docstring region heuristically: take text after the last triple-quote
    if '"""' in src:
        code = src.split('"""')[-1]
    for tok in DEPMAP_TOKENS:
        assert tok not in code, f"nonconformity code references referee token {tok!r}"


def test_nonconformity_uses_only_structural_columns():
    """nonconformity() output must be invariant to any DepMap-like column added."""
    df = pd.DataFrame({"score": np.linspace(0.1, 0.9, 20)})
    base = nonconformity(df)
    df2 = df.copy()
    df2["gls_p"] = np.random.default_rng(0).random(20)       # inject a referee column
    df2["coess"] = np.random.default_rng(1).random(20)
    after = nonconformity(df2)
    assert np.allclose(base, after), "nonconformity changed when a DepMap column was added -> LEAK"


def test_labels_parquet_carries_no_depmap_column():
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    for col in lab.columns:
        assert not any(t in col.lower() for t in DEPMAP_TOKENS), f"label column {col!r} looks DepMap-derived"


# ---------- certification is invariant to permuting the referee ----------
def test_certification_invariant_to_depmap_permutation():
    """Permuting DepMap values must not change the certified set (referee only RANKS)."""
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    s_cal = nonconformity(lab[(lab.label == 0) & (lab.split == "cal")])
    cand = inter[inter["is_true_pair"]].copy()
    p = conformal_pvalues(nonconformity(cand), s_cal)
    cert_a = benjamini_hochberg(p, C.Q)
    # The pipeline never reads DepMap to certify; certification depends only on `score`.
    # Shuffle the candidate `score`? No — that changes the input. Instead confirm the
    # certification function's inputs contain no DepMap channel: recompute with a
    # DepMap column present and permuted.
    cand2 = cand.copy()
    rng = np.random.default_rng(7)
    cand2["gls_p"] = rng.random(len(cand2))
    p2 = conformal_pvalues(nonconformity(cand2), s_cal)
    cert_b = benjamini_hochberg(p2, C.Q)
    assert np.array_equal(cert_a, cert_b), "certified set moved when DepMap column permuted -> LEAK"


# ---------- no-signal negative control ----------
def test_no_signal_control_collapses_enrichment():
    """No-signal negative control (PLAN_V3 T0.3). Keep the REAL certified/dropped
    partition, but permute the co-essentiality REFEREE across candidate pairs (break
    the edge<->co-essentiality correspondence). Under a shuffled referee the
    certified-vs-dropped gap must collapse to ~0 and the real gap must sit in the
    upper tail; if a shuffled referee reproduced the gap, the "enrichment" would be
    an artifact of the partition sizes rather than real biology.

    (Shuffling the confidence SCORE instead is not a clean null here: `dropped` is
    conditioned on the paper's high_conf flag, itself co-essentiality-enriched, so
    the score-shuffle null is asymmetric by construction. Permuting the referee is
    the correct signal-destroying control.)"""
    slice_path = C.INTERIM / "depmap_slice.npz"
    if not slice_path.exists():
        pytest.skip("depmap_slice.npz not present")
    z = np.load(slice_path, allow_pickle=True)
    genes = list(z["genes"]); gidx = {g: i for i, g in enumerate(genes)}
    gls_p, gls_sign = z["gls_p"], z["gls_sign"]

    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    s_cal = nonconformity(lab[(lab.label == 0) & (lab.split == "cal")])
    cand = inter[inter["is_true_pair"]].copy()

    # Per-edge "strong co-essential" indicator for edges with both genes in DepMap.
    strong, in_dm = [], []
    for a, b in zip(cand.gene_a, cand.gene_b):
        if a in gidx and b in gidx and a != b:
            i, j = gidx[a], gidx[b]
            strong.append(bool(gls_p[i, j] <= C.COESS_SIG_THRESHOLD and gls_sign[i, j] > 0))
            in_dm.append(True)
        else:
            strong.append(False); in_dm.append(False)
    strong = np.array(strong); in_dm = np.array(in_dm)

    p = conformal_pvalues(nonconformity(cand), s_cal)
    cert = benjamini_hochberg(p, C.Q)
    hc = cand["high_conf"].to_numpy()
    drop = hc & ~cert

    def gap(strong_vec):
        cs = strong_vec[cert & in_dm]; ds = strong_vec[drop & in_dm]
        cf = cs.mean() if cs.size else 0.0
        df_ = ds.mean() if ds.size else 0.0
        return cf - df_

    real_gap = gap(strong)

    rng = np.random.default_rng(0)
    # permute the referee only among DepMap-scored edges (preserve the base rate)
    idx_dm = np.nonzero(in_dm)[0]
    null_gaps = []
    for _ in range(500):
        perm_strong = strong.copy()
        perm_strong[idx_dm] = rng.permutation(strong[idx_dm])
        null_gaps.append(gap(perm_strong))
    null_gaps = np.array(null_gaps)
    perm_p = float((null_gaps >= real_gap).mean())

    assert real_gap > 0.1, f"real enrichment gap unexpectedly small: {real_gap}"
    assert perm_p < 0.05, f"real gap not distinguishable from shuffled-referee null (perm p={perm_p})"
    assert abs(np.mean(null_gaps)) < 0.05, f"shuffled-referee null gap not centered near 0: {np.mean(null_gaps)}"
