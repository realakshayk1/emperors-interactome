"""prevalence.py — the prevalence-shift centerpiece (METHODS §5).

Benchmark-estimated FDR (the paper's approach) fixes a score cutoff that gives
FDR=q on a BALANCED benchmark (positive:negative ~1:1). That cutoff's *realized*
FDR climbs above q once true interactions are rare, because the benchmark's
assumed prevalence no longer matches the test pool. Conformal p-values + BH
re-derive the selection from the calibration negatives at each composition, so
the realized FDR stays <= q regardless of prevalence.

We hold positives fixed and resample negatives to sweep prevalence
pi = n_pos / (n_pos + n_neg), averaging realized FDR over many resamples.
Output: data/processed/prevalence_shift.json. Run: python -m emperor.prevalence
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .conformal import benjamini_hochberg, conformal_pvalues
from .nonconformity import nonconformity


def _benchmark_cutoff(score_pos, score_neg, q) -> float:
    """LOWEST score cutoff t whose balanced-benchmark FDR <= q (maximizes recall
    at the claimed error rate — how a benchmark-FDR method sets its threshold).
    Balanced FDR(t) = w*#{neg>=t} / (#{pos>=t} + w*#{neg>=t}), w = n_pos/n_neg
    so the benchmark is effectively 1:1."""
    w = len(score_pos) / max(len(score_neg), 1)
    cuts = np.unique(np.concatenate([score_pos, score_neg]))  # ascending
    best = None
    for t in cuts:  # from low to high; take the LOWEST t that satisfies q
        tp = (score_pos >= t).sum()
        fp = (score_neg >= t).sum() * w
        if tp + fp > 0 and fp / (tp + fp) <= q:
            best = float(t)
            break
    return best if best is not None else float(cuts.max())


def run(n_resamples: int = 300):
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    q = C.Q
    rng = np.random.default_rng(C.SEED)

    pos = lab[lab.label == 1]
    neg = lab[lab.label == 0]
    # calibration negatives (for conformal) = cal-split negatives
    cal_neg = neg[neg.split == "cal"]
    s_cal_neg = nonconformity(cal_neg)
    # benchmark cutoff set on cal split (balanced pos vs neg)
    cutoff = _benchmark_cutoff(pos[pos.split == "cal"].score.values,
                               cal_neg.score.values, q)

    # test material: test-split positives + test-split negatives (disjoint from cal)
    pos_te = pos[pos.split == "test"]
    neg_te = neg[neg.split == "test"]
    s_pos = pos_te.score.values
    s_neg_pool = neg_te.score.values
    n_pos = len(s_pos)

    prevalences = [0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01]
    bench_realized, conf_realized = [], []
    for pi in prevalences:
        n_neg = int(round(n_pos * (1 - pi) / pi))
        br, cr = [], []
        for _ in range(n_resamples):
            neg_s = rng.choice(s_neg_pool, size=n_neg, replace=True)
            s_test = np.concatenate([s_pos, neg_s])
            is_neg = np.concatenate([np.zeros(n_pos, bool), np.ones(n_neg, bool)])
            # benchmark: fixed cutoff on score
            sel_b = s_test >= cutoff
            br.append(is_neg[sel_b].mean() if sel_b.sum() else 0.0)
            # conformal: BH on conformal p-values from cal negatives
            p = conformal_pvalues(1.0 - s_test, s_cal_neg)
            sel_c = benjamini_hochberg(p, q)
            cr.append(is_neg[sel_c].mean() if sel_c.sum() else 0.0)
        bench_realized.append(float(np.mean(br)))
        conf_realized.append(float(np.mean(cr)))

    out = dict(q=q, prevalences=prevalences, benchmark_cutoff=cutoff,
               benchmark_fdr=bench_realized, conformal_fdr=conf_realized,
               n_resamples=n_resamples)
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "prevalence_shift.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"benchmark cutoff (FDR={o['q']} on balanced set): score>={o['benchmark_cutoff']:.3f}")
    print(f"{'prevalence':>10} {'benchmark realized FDR':>24} {'conformal realized FDR':>24}")
    for pi, b, c in zip(o["prevalences"], o["benchmark_fdr"], o["conformal_fdr"]):
        print(f"{pi:>10.2f} {b:>24.3f} {c:>24.3f}")


if __name__ == "__main__":
    main()
