"""shift_attribution.py — attribute the calibration->wild covariate shift (PLAN_V3 T1.1).

Pre-registered gate (DECISIONS.md 2026-07-10): produce a quantitative attribution of
how much of the calibration-decoy <-> wild-candidate divergence is explained by
`score` alone vs by the fuller local covariate vector. The gate PASSES ON PRODUCING
the attribution, whatever it shows; the number decides which WS1 fix leads (hard
negatives vs conditional/Gamma).

Method:
  * Per-covariate 2-sample KS(decoy, candidate) with the same-distribution decoy-vs-
    decoy floor (KS ~0.026 from wcs_results) as a reference.
  * Density-ratio (logistic classifier decoy=0 / candidate=1) fit on (i) `score`
    alone and (ii) the full covariate vector; report classifier AUC (a proper
    divergence proxy — AUC 0.5 = no shift), the implied importance-weight ESS
    fraction, and max weight. The AUC gap between (ii) and (i) is the share of the
    shift INVISIBLE to a 1-D `score` reweighting -> the WCS-failure diagnosis.
  * Nearest-neighbour coverage: fraction of wild candidates whose nearest decoy (in
    standardized covariate space) is within the decoy-decoy NN radius — i.e. does the
    decoy null SPAN the region the candidates occupy? Poor coverage => no reweighting
    of the existing null can fix control (need a better null): the T1.2 rationale.

Purity: covariates are AF/graph only (covariates.py) — no DepMap.
Run: python -m emperor.shift_attribution
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score

from . import config as C
from .covariates import build_covariate_frame, LOCAL_COVARIATES, covariate_matrix


def _density_ratio_auc(df, cols, seed=C.SEED):
    """5-fold CV AUC of a decoy(0)/candidate(1) classifier on `cols`, plus weight ESS."""
    y = (df["pool"].to_numpy() == "candidate").astype(int)
    X = covariate_matrix(df, cols)
    clf = LogisticRegression(max_iter=1000, C=1.0)
    proba = cross_val_predict(clf, X, y, cv=5, method="predict_proba")[:, 1]
    auc = float(roc_auc_score(y, proba))
    # importance weights w = p(candidate)/p(decoy) evaluated on the DECOY rows
    pdec = np.clip(proba[y == 0], 1e-6, 1 - 1e-6)
    w = pdec / (1 - pdec)
    w = w / w.mean()
    ess_frac = float((w.sum() ** 2) / (len(w) * (w ** 2).sum()))
    return auc, ess_frac, float(w.max())


def _nn_coverage(df, cols, seed=C.SEED):
    """Fraction of candidates within the decoy-decoy NN radius of their nearest decoy."""
    from scipy.spatial import cKDTree
    dec = df[df.pool == "decoy"]; cand = df[df.pool == "candidate"]
    Xall = covariate_matrix(df, cols)
    dmask = (df.pool == "decoy").to_numpy()
    Xdec = Xall[dmask]; Xcand = Xall[~dmask]
    tree = cKDTree(Xdec)
    # decoy-decoy NN radius distribution (exclude self -> k=2)
    dd, _ = tree.query(Xdec, k=2)
    dd = dd[:, 1]
    radius = np.quantile(dd, 0.95)          # 95th pct of decoy-decoy NN distance
    dc, _ = tree.query(Xcand, k=1)
    covered = float((dc <= radius).mean())
    return covered, float(radius), float(np.median(dc)), float(np.median(dd))


def run():
    df = build_covariate_frame()
    dec = df[df.pool == "decoy"]; cand = df[df.pool == "candidate"]

    # per-covariate KS
    ks = {}
    for c in LOCAL_COVARIATES:
        d, p = stats.ks_2samp(dec[c].to_numpy(float), cand[c].to_numpy(float))
        ks[c] = dict(ks=float(d), p=float(p))

    auc_score, ess_score, wmax_score = _density_ratio_auc(df, ["score"])
    auc_full, ess_full, wmax_full = _density_ratio_auc(df, LOCAL_COVARIATES)
    cov_score, r_s, mdc_s, mdd_s = _nn_coverage(df, ["score"])
    cov_full, r_f, mdc_f, mdd_f = _nn_coverage(df, LOCAL_COVARIATES)

    # Share of divergence invisible to a 1-D score reweighting:
    # measured as extra AUC above the score-only classifier, scaled by the total
    # divergence above chance.
    div_score = auc_score - 0.5
    div_full = auc_full - 0.5
    invisible_frac = float((div_full - div_score) / div_full) if div_full > 0 else 0.0

    out = dict(
        per_covariate_ks=ks,
        density_ratio=dict(
            score_only=dict(auc=auc_score, ess_frac=ess_score, weight_max=wmax_score),
            full_covariate=dict(auc=auc_full, ess_frac=ess_full, weight_max=wmax_full),
            divergence_invisible_to_score_frac=invisible_frac,
        ),
        nn_coverage=dict(
            score_only=dict(covered_frac=cov_score, decoy_radius=r_s,
                            median_cand_nn=mdc_s, median_decoy_nn=mdd_s),
            full_covariate=dict(covered_frac=cov_full, decoy_radius=r_f,
                                median_cand_nn=mdc_f, median_decoy_nn=mdd_f),
        ),
        interpretation=(
            "Full-covariate classifier AUC %.3f vs score-only %.3f: %.0f%% of the "
            "cal<->wild divergence is INVISIBLE to a 1-D score reweighting (the WCS "
            "failure is a null-completeness problem, not a reweighting bug). "
            "Full-covariate NN coverage of wild candidates by decoys = %.0f%% -> the "
            "decoy null %s span the candidate region."
        ) % (auc_full, auc_score, 100 * invisible_frac, 100 * cov_full,
             "does NOT" if cov_full < 0.8 else "does"),
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "shift_attribution.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print("Per-covariate KS(decoy, candidate):")
    for c, d in o["per_covariate_ks"].items():
        print(f"  {c:14s} KS={d['ks']:.3f}  p={d['p']:.2e}")
    dr = o["density_ratio"]
    print(f"\nDensity-ratio classifier AUC: score-only={dr['score_only']['auc']:.3f} "
          f"(ESS {dr['score_only']['ess_frac']:.2f}), "
          f"full={dr['full_covariate']['auc']:.3f} (ESS {dr['full_covariate']['ess_frac']:.2f})")
    print(f"  -> {100*dr['divergence_invisible_to_score_frac']:.0f}% of divergence invisible to 1-D score reweighting")
    nn = o["nn_coverage"]
    print(f"NN coverage of candidates by decoys: score-only={100*nn['score_only']['covered_frac']:.0f}%  "
          f"full-covariate={100*nn['full_covariate']['covered_frac']:.0f}%")
    print("\n" + o["interpretation"])


if __name__ == "__main__":
    main()
