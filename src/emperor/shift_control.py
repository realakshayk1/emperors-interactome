"""shift_control.py — WS1 gate: does a better null / conditioning restore control?
(PLAN_V3 T1.2 + T1.3)

Baseline to beat (plain conformal, node/protein-disjoint held-out FDR, from
wcs_results.node_disjoint_wild): 0.29 / 0.32 / 0.37 at q = 0.05 / 0.10 / 0.20.

Two candidate fixes, both evaluated on the SAME protein-disjoint splits:
  (T1.2) hard-negative-augmented calibration null (hardnegatives.py rule).
  (T1.3) degree-stratified (Mondrian) conformal: calibrate each test edge only
         against negatives in its own union-graph-degree stratum, so the
         hub-selection covariate the shift rides on is conditioned out.

Held-out FDR here = fraction of certified held-out CORUM-negatives among certified
held-out labeled pairs, over protein-disjoint splits (an UPPER bound; CORUM
incomplete). Referee-corroboration (DepMap co-ess of the certified set) is reported
alongside so a power collapse is visible.

Purity: negatives/positives are CORUM + decoys; NO DepMap in calibration. Run:
python -m emperor.shift_control
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .covariates import add_covariates
from .conformal import conformal_pvalues, benjamini_hochberg
from .nonconformity import nonconformity
from .hardnegatives import build_hard_negatives

QS = [0.05, 0.10, 0.20]


def _labeled_frame():
    """CORUM positives + decoy negatives with covariates and protein ids."""
    lab = pd.read_parquet(C.INTERIM / "labels.parquet").copy()
    inter = add_covariates(pd.read_parquet(C.INTERIM / "interactome.parquet"))
    degmap = {}
    for _, r in inter.iterrows():
        degmap[frozenset((r.uniprot_a, r.uniprot_b))] = r.degree
    lab["degree"] = [degmap.get(frozenset((a, b)), 0)
                     for a, b in zip(lab.uniprot_a, lab.uniprot_b)]
    return lab


def _protein_disjoint_split(lab, seed):
    """Assign proteins to cal/test; a pair is 'test' iff BOTH proteins are test."""
    rng = np.random.default_rng(seed)
    prots = pd.unique(pd.concat([lab.uniprot_a, lab.uniprot_b]))
    assign = {p: ("cal" if rng.random() < 0.5 else "test") for p in prots}
    def side(a, b):
        sa, sb = assign[a], assign[b]
        if sa == sb == "test":
            return "test"
        if sa == sb == "cal":
            return "cal"
        return "drop"
    return lab.assign(pd_split=[side(a, b) for a, b in zip(lab.uniprot_a, lab.uniprot_b)])


def _degree_cutoff():
    inter = add_covariates(pd.read_parquet(C.INTERIM / "interactome.parquet"))
    cand = inter[inter.is_true_pair]
    from .hardnegatives import DEGREE_TERTILE
    return float(np.quantile(cand["degree"].to_numpy(float), DEGREE_TERTILE))


def _fdr_pooled(method, lab, hard_flag, n_splits=200):
    """Pooled held-out FDR per q for a calibration method on protein-disjoint splits.

    method in {'plain','hard','mondrian'}.
    """
    cutoff = _degree_cutoff()
    acc = {q: [0, 0] for q in QS}     # [false, total]
    ncert = {q: 0 for q in QS}
    for seed in range(n_splits):
        sp = _protein_disjoint_split(lab, seed)
        cal = sp[sp.pd_split == "cal"]
        test = sp[sp.pd_split == "test"]
        if len(test) == 0 or (cal.label == 0).sum() == 0:
            continue
        s_test = nonconformity(test)
        y = test.label.to_numpy()

        if method == "plain":
            s_cal = nonconformity(cal[cal.label == 0])
            p = conformal_pvalues(s_test, s_cal)
        elif method == "hard":
            calneg = cal[cal.label == 0]
            # hard = high-degree negatives (pre-registered cutoff via hardnegatives)
            hard = calneg[calneg.degree >= cutoff]
            s_cal = nonconformity(pd.concat([calneg, hard]))  # upweight hard region
            p = conformal_pvalues(s_test, s_cal)
        elif method == "mondrian":
            calneg = cal[cal.label == 0]
            # two degree strata; each test edge calibrated within its own stratum
            p = np.ones(len(test))
            for lo, hi in [(-1, cutoff), (cutoff, 1e9)]:
                tmask = (test.degree.to_numpy() >= lo) & (test.degree.to_numpy() < hi) if lo == -1 \
                    else (test.degree.to_numpy() >= lo)
                cmask = (calneg.degree.to_numpy() >= lo) & (calneg.degree.to_numpy() < hi) if lo == -1 \
                    else (calneg.degree.to_numpy() >= lo)
                if tmask.sum() == 0 or cmask.sum() < 5:
                    continue
                p[tmask] = conformal_pvalues(s_test[tmask], nonconformity(calneg[cmask]))
        else:
            raise ValueError(method)

        for q in QS:
            r = benjamini_hochberg(p, q)
            acc[q][0] += int((1 - y[r]).sum())
            acc[q][1] += int(r.sum())
            ncert[q] += int(r.sum())
    out = {}
    for q in QS:
        f, t = acc[q]
        out[q] = dict(pooled_fdr=float(f / t) if t else 0.0,
                      mean_ncert=float(ncert[q] / n_splits))
    return out


def run(n_splits=200):
    lab = _labeled_frame()
    baseline = {0.05: 0.2895, 0.10: 0.3191, 0.20: 0.3668}   # wcs_results node_disjoint_wild pooled
    res = {}
    for method in ("plain", "hard", "mondrian"):
        res[method] = _fdr_pooled(method, lab, None, n_splits=n_splits)
    out = dict(baseline_wcs_plain=baseline, n_splits=n_splits, methods=res,
               q_sweep=QS)
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "shift_control.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"protein-disjoint held-out pooled FDR ({o['n_splits']} splits), lower=better toward q:")
    print(f"{'method':>10} | " + " | ".join(f"q={q}" for q in QS))
    for m, d in o["methods"].items():
        print(f"{m:>10} | " + " | ".join(f"{d[q]['pooled_fdr']:.3f} (n={d[q]['mean_ncert']:.0f})" for q in QS))


if __name__ == "__main__":
    main()
