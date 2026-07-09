"""conformal.py — conformal p-values + Benjamini-Hochberg FDR control (METHODS §5).

For a candidate edge j with nonconformity s_j, and calibration NEGATIVES
{s_k : label=0} of size n, the (marginal) conformal p-value is

    p_j = (1 + #{k : s_k <= s_j}) / (n + 1)

Under exchangeability of a true non-interaction with the calibration negatives,
p_j is a valid p-value for H0: "j is a non-interaction" (Jin & Candès 2023;
Marandon 2023 conformal link-prediction FDR). Applying Benjamini-Hochberg to the
{p_j} at level q yields a rejected set with FDR <= q, distribution-free.

Run the audit: python -m emperor.conformal
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .nonconformity import nonconformity


def conformal_pvalues(s_test: np.ndarray, s_cal_neg: np.ndarray) -> np.ndarray:
    """Marginal conformal p-values: (1 + #{cal_neg <= s_j}) / (n+1)."""
    cal = np.sort(np.asarray(s_cal_neg, dtype=float))
    n = cal.size
    # #{k : s_cal <= s_j} via searchsorted on the sorted calibration scores
    counts = np.searchsorted(cal, np.asarray(s_test, dtype=float), side="right")
    return (1.0 + counts) / (n + 1.0)


def benjamini_hochberg(p: np.ndarray, q: float) -> np.ndarray:
    """Return a boolean rejection mask controlling FDR at q (BH step-up)."""
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, m + 1) / m)
    passed = ranked <= thresh
    reject = np.zeros(m, dtype=bool)
    if passed.any():
        kmax = np.nonzero(passed)[0].max()
        reject[order[: kmax + 1]] = True
    return reject


def certify(df_test: pd.DataFrame, s_cal_neg: np.ndarray, q: float,
            w_phys: float = 1.0) -> pd.DataFrame:
    """Certify edges in df_test at FDR level q using calibration negatives."""
    out = df_test.copy()
    out["nonconf"] = nonconformity(out, w_phys=w_phys)
    out["conf_pvalue"] = conformal_pvalues(out["nonconf"].to_numpy(), s_cal_neg)
    out[f"certified@{q}"] = benjamini_hochberg(out["conf_pvalue"].to_numpy(), q)
    return out


def _empirical_fdr(certified_mask, is_true) -> float:
    """FDR proxy on a labeled test split: among certified, fraction that are
    NOT true interactions (using the CORUM-derived label of the test pool)."""
    cert = np.asarray(certified_mask, dtype=bool)
    if cert.sum() == 0:
        return 0.0
    return float((~np.asarray(is_true, dtype=bool))[cert].mean())


def _heldout_fdr_montecarlo(lab: pd.DataFrame, n_splits: int = 200) -> dict:
    """Repeated complex-disjoint split; pooled + mean held-out FDR per q."""
    pos = lab[lab.label == 1].copy()
    neg = lab[lab.label == 0].copy()
    counts = pos["complex_id"].value_counts().to_dict()
    qs = sorted(set(C.Q_SWEEP) | {C.Q})
    pooled = {q: [0, 0] for q in qs}       # [false, total]
    per_split = {q: [] for q in qs}
    for seed in range(n_splits):
        rng = np.random.default_rng(seed)
        cxs = pos["complex_id"].dropna().unique().tolist()
        rng.shuffle(cxs)
        cal_c, cn, tn = set(), 0, 0
        for c in sorted(cxs, key=lambda c: -counts[c]):
            if cn <= tn:
                cal_c.add(c); cn += counts[c]
            else:
                tn += counts[c]
        pos_sp = pos["complex_id"].map(lambda c: "cal" if c in cal_c else "test")
        neg_sp = np.where(rng.random(len(neg)) < 0.5, "cal", "test")
        s_calneg = nonconformity(neg[neg_sp == "cal"])
        test = pd.concat([pos[pos_sp == "test"], neg[neg_sp == "test"]])
        p = conformal_pvalues(nonconformity(test), s_calneg)
        y = test["label"].to_numpy()
        for q in qs:
            rej = benjamini_hochberg(p, q)
            pooled[q][0] += int((1 - y[rej]).sum())
            pooled[q][1] += int(rej.sum())
            if rej.sum():
                per_split[q].append(float((1 - y[rej]).mean()))
    out = {}
    for q in qs:
        f, t = pooled[q]
        out[q] = dict(
            pooled_fdr=float(f / t) if t else 0.0,
            mean_split_fdr=float(np.mean(per_split[q])) if per_split[q] else 0.0,
            mean_n_certified=float(t / n_splits),
            n_splits_with_cert=len(per_split[q]),
        )
    return out


def run():
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")

    # Calibration negatives = label-0 pairs in the CAL split
    cal_neg = lab[(lab.label == 0) & (lab.split == "cal")]
    s_cal_neg = nonconformity(cal_neg)

    # Candidate pool for certification = the 1,666 real ("True") candidate edges
    cand = inter[inter["is_true_pair"]].copy()

    # Sweep q; produce certified.parquet at the primary q with all sweep columns
    result = certify(cand, s_cal_neg, C.Q)
    sweep = {}
    for q in sorted(set(C.Q_SWEEP) | {C.Q})[::-1]:
        col = f"certified@{q}"
        result[col] = benjamini_hochberg(result["conf_pvalue"].to_numpy(), q)
        sweep[q] = int(result[col].sum())

    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    dest = C.PROCESSED / "certified.parquet"
    result.to_parquet(dest, index=False)

    # --- Validity: held-out empirical FDR, Monte-Carlo over repeated
    # complex-disjoint splits (a single split gives a noisy point estimate).
    # A "false discovery" = a certified CORUM-label-0 (random-decoy) pair. Because
    # CORUM is incomplete, some certified decoys may be real interactions -> this
    # estimate is an UPPER BOUND on the true FDR (negative-set problem, LEARNINGS).
    fdr_by_q = _heldout_fdr_montecarlo(lab, n_splits=200)

    n_hc = int(cand["high_conf"].sum())
    hc_dropped = {}
    for q in sorted(set(C.Q_SWEEP) | {C.Q}):
        col = f"certified@{q}"
        dropped = int((cand["high_conf"].to_numpy() & ~result[col].to_numpy()).sum())
        hc_dropped[q] = dropped

    summary = dict(
        n_candidates=int(len(cand)), n_high_conf=n_hc,
        certified_by_q=sweep, high_conf_dropped_by_q=hc_dropped,
        heldout_fdr_by_q=fdr_by_q, q_primary=C.Q,
    )
    (C.PROCESSED / "audit_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    s = run()
    print(f"candidates={s['n_candidates']}  high-confidence(paper)={s['n_high_conf']}")
    print("certified by q:", s["certified_by_q"])
    print("high-conf edges DROPPED by q:", s["high_conf_dropped_by_q"])
    print("held-out empirical FDR by q (Monte-Carlo, UPPER bound; CORUM incomplete):")
    for q, d in s["heldout_fdr_by_q"].items():
        print(f"  q={q}: pooled_FDR={d['pooled_fdr']:.3f}  mean_split_FDR={d['mean_split_fdr']:.3f}"
              f"  mean_n_cert={d['mean_n_certified']:.0f}")


if __name__ == "__main__":
    main()
