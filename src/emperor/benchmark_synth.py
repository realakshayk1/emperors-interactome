"""benchmark_synth.py — semi-synthetic controlled-prevalence benchmark (PLAN_V3 T2.2).

The held-out CORUM FDR is only an UPPER bound (CORUM incomplete). Here we build a
benchmark where TRUTH IS KNOWN: draw positives from CORUM same-complex pairs and
negatives from random decoys (the two label pools the map already ships), assign
realistic per-class AF-score laws fit from the REAL data, then inject them at
controlled prevalences pi in {0.30,0.10,0.05,0.02}. Because membership is assigned
by us, the realized FDR of any selection rule is exactly computable.

We compare, at each pi:
  * conformal + BH  (should hold realized FDR <= q + finite-sample slack at all pi)
  * benchmark-tuned fixed cutoff (calibrated to q on a BALANCED set; should inflate
    as pi drops — the prevalence wedge, now on known-truth data)

Score laws: fit a simple parametric model (Beta on [0,1]) to the REAL positive
(CORUM same-complex) and negative (decoy) `score` distributions, so the synthetic
scores are realistic, not arbitrary. Purity: no DepMap. Run: python -m emperor.benchmark_synth
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .conformal import conformal_pvalues, benjamini_hochberg

PREVALENCES = [0.30, 0.10, 0.05, 0.02]
QS = [0.05, 0.10, 0.20]


def _fit_beta(x):
    x = np.clip(np.asarray(x, float), 1e-4, 1 - 1e-4)
    a, b, loc, scale = stats.beta.fit(x, floc=0, fscale=1)
    return a, b


def _score_laws():
    """Fit Beta laws to REAL positive vs negative `score` distributions."""
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    pos = lab[lab.label == 1]["score"].to_numpy()
    neg = lab[lab.label == 0]["score"].to_numpy()
    return _fit_beta(pos), _fit_beta(neg), (pos.mean(), neg.mean())


def _one_run(pos_ab, neg_ab, pi, q, n_test=4000, n_cal=1500, seed=0):
    rng = np.random.default_rng(seed)
    # calibration negatives (pure nulls) — score law = negatives
    s_cal = 1.0 - stats.beta.rvs(*neg_ab, size=n_cal, random_state=rng)
    # test pool at prevalence pi
    n_pos = rng.binomial(n_test, pi)
    n_neg = n_test - n_pos
    score_pos = stats.beta.rvs(*pos_ab, size=n_pos, random_state=rng)
    score_neg = stats.beta.rvs(*neg_ab, size=n_neg, random_state=rng)
    scores = np.concatenate([score_pos, score_neg])
    truth = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(bool)
    s_test = 1.0 - scores

    # conformal + BH
    p = conformal_pvalues(s_test, s_cal)
    rej = benjamini_hochberg(p, q)
    conf_fdr = float((~truth[rej]).mean()) if rej.sum() else 0.0
    conf_power = float(truth[rej].sum() / max(1, truth.sum()))

    # benchmark cutoff: calibrated on a BALANCED (pi=0.5) set to give FDR=q, then
    # applied at this pi (the prevalence-fragile lookup).
    nb = 2000
    bal_pos = stats.beta.rvs(*pos_ab, size=nb, random_state=rng)
    bal_neg = stats.beta.rvs(*neg_ab, size=nb, random_state=rng)
    # choose cutoff c so that among balanced predicted-positives (score>=c) FDR=q
    cand_c = np.linspace(0.1, 0.95, 200)
    best_c = 0.95
    for c in cand_c:
        tp = (bal_pos >= c).sum(); fp = (bal_neg >= c).sum()
        if tp + fp > 0 and fp / (tp + fp) <= q:
            best_c = c
            break
    sel = scores >= best_c
    bench_fdr = float((~truth[sel]).mean()) if sel.sum() else 0.0
    bench_power = float(truth[sel].sum() / max(1, truth.sum()))
    return conf_fdr, conf_power, bench_fdr, bench_power


def run(n_seeds=100):
    pos_ab, neg_ab, means = _score_laws()
    rows = []
    for pi in PREVALENCES:
        for q in QS:
            cf, cp, bf, bp = [], [], [], []
            for seed in range(n_seeds):
                a, b, c, d = _one_run(pos_ab, neg_ab, pi, q, seed=seed)
                cf.append(a); cp.append(b); bf.append(c); bp.append(d)
            rows.append(dict(
                prevalence=pi, q=q,
                conformal_fdr=float(np.mean(cf)), conformal_fdr_se=float(np.std(cf) / np.sqrt(n_seeds)),
                conformal_power=float(np.mean(cp)),
                benchmark_fdr=float(np.mean(bf)), benchmark_fdr_se=float(np.std(bf) / np.sqrt(n_seeds)),
                benchmark_power=float(np.mean(bp)),
                conformal_controls=bool(np.mean(cf) <= q + 0.03),
            ))
    out = dict(
        pos_beta=list(pos_ab), neg_beta=list(neg_ab),
        real_score_means=dict(pos=float(means[0]), neg=float(means[1])),
        prevalences=PREVALENCES, q_sweep=QS, n_seeds=n_seeds, rows=rows,
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "benchmark_synth.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"Semi-synthetic benchmark (Beta laws fit to real scores; {o['n_seeds']} seeds)")
    print(f"real score means: pos={o['real_score_means']['pos']:.3f} neg={o['real_score_means']['neg']:.3f}")
    print(f"\n{'pi':>5} {'q':>5} | {'conf FDR':>9} {'ctrl?':>5} {'pow':>5} | {'bench FDR':>9} {'pow':>5}")
    for r in o["rows"]:
        print(f"{r['prevalence']:>5} {r['q']:>5} | {r['conformal_fdr']:>9.3f} "
              f"{'YES' if r['conformal_controls'] else 'NO':>5} {r['conformal_power']:>5.2f} | "
              f"{r['benchmark_fdr']:>9.3f} {r['benchmark_power']:>5.2f}")


if __name__ == "__main__":
    main()
