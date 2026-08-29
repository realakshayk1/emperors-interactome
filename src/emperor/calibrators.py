"""calibrators.py — the feasibility bound, and how much of it the calibrator controls.

dependence.py recomputes the audit under BY and under e-BH with a *harmonic*-calibrated
conformal e-value, and reports zero certifications for both. That zero is not a property of
e-BH. It is a property of the calibrator.

Proposition (feasibility). Fix m hypotheses and level q, and write H_k = sum_{i<=k} 1/i.
  (i)  If every p_j >= F, the BY step-up condition p_(i) <= i q / (m H_m) fails for every
       i < ceil(m H_m F / q).
  (ii) If every e_j <= e_max, the e-BH condition e_(k) >= m / (k q) fails for every
       k < ceil(m / (q e_max)).
In either case the procedure rejects at least k* hypotheses or none, and rejects none
whatever the data once k* > m.

For marginal conformal p-values F = 1/(n_cal+1) and H_m is fixed by the number of tests, so
BY's bound admits no remedy. For conformal e-values e_max belongs to the calibrator:

  harmonic      e_j = (n+1) / (R_j H_{n+1})              e_max = (n+1)/H_{n+1}
  threshold(t)  e_j = ((n+1)/t) * 1{R_j <= t}            e_max = (n+1)/t

Both are valid e-values under exchangeability, since R_j is stochastically at least uniform
on {1..n+1}, giving P(R_j <= t) <= t/(n+1) and hence E[e_j] <= 1. The harmonic calibrator
pays the arbitrary-dependence logarithm a second time, by hand; the threshold calibrator does
not, and on this map it recovers BH's discovery count without BH's PRDS assumption.

The threshold index t is a design parameter and must be fixed before the data are seen.
This module reports t in {1, 2} so the dependence of the result on that choice is visible.

Run: python -m emperor.calibrators
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from . import config as C
from .conformal import benjamini_hochberg
from .dependence import _harmonic, benjamini_yekutieli, ebh

QS = (0.05, 0.10, 0.20)
THRESHOLDS = (1, 2)


def threshold_evalues(rank: np.ndarray, n_cal: int, t: int) -> np.ndarray:
    """Threshold-calibrated conformal e-values (Bashari et al. 2023)."""
    return np.where(rank <= t, (n_cal + 1.0) / t, 0.0)


def harmonic_evalues(rank: np.ndarray, n_cal: int) -> np.ndarray:
    """Harmonic-calibrated conformal e-values, as used in dependence.py."""
    return (n_cal + 1.0) / (rank * _harmonic(n_cal + 1))


def run() -> dict:
    cert = pd.read_parquet(C.PROCESSED / "certified.parquet")
    dep = json.loads((C.PROCESSED / "dependence_robustness.json").read_text())

    p = cert["conf_pvalue"].to_numpy(dtype=float)
    hi = cert["high_conf"].to_numpy().astype(bool)
    m = p.size
    n_cal = int(dep["n_cal_neg"])
    n1 = n_cal + 1

    # conformal p-values live on the grid {1/(n+1), ..., 1}, so the calibration rank is exact
    rank = np.rint(p * n1).astype(int)
    H_m, H_n1 = _harmonic(m), _harmonic(n1)

    def counts(reject: np.ndarray) -> tuple[int, int]:
        return int(reject.sum()), int((hi & ~reject).sum())

    out: dict = {
        "m": m, "n_cal_plus_1": n1, "n_high_conf": int(hi.sum()),
        "H_m": H_m, "H_n_cal_plus_1": H_n1,
        "n_at_rank": {str(t): int((rank <= t).sum()) for t in (1, 2, 3)},
        "k_star_BY": {str(q): math.ceil(m * H_m / (q * n1)) for q in QS},
        "procedures": {},
    }

    for name, sel, dep_kind in (("BH", benjamini_hochberg, "PRDS"),
                                ("BY", benjamini_yekutieli, "arbitrary")):
        rec = {"dependence": dep_kind, "certified": {}, "dropped": {}}
        for q in QS:
            c, d = counts(sel(p, q))
            rec["certified"][str(q)], rec["dropped"][str(q)] = c, d
        out["procedures"][name] = rec

    cals = [("harmonic", harmonic_evalues(rank, n_cal), n1 / H_n1)]
    cals += [(f"threshold_t{t}", threshold_evalues(rank, n_cal, t), n1 / t) for t in THRESHOLDS]

    for name, e, e_max in cals:
        rec = {"dependence": "arbitrary", "e_max": float(e_max),
               "k_star": {str(q): math.ceil(m / (q * e_max)) for q in QS},
               "certified": {}, "dropped": {}}
        for q in QS:
            c, d = counts(ebh(e, q))
            rec["certified"][str(q)], rec["dropped"][str(q)] = c, d
        out["procedures"]["e-BH " + name] = rec

    (C.PROCESSED / "calibrator_comparison.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    o = run()
    print(f"m={o['m']}  n_cal+1={o['n_cal_plus_1']}  "
          f"at floor={o['n_at_rank']['1']}  rank<=2={o['n_at_rank']['2']}")
    print(f"{'procedure':24s}{'e_max':>9s}   certified 0.05/0.10/0.20      k*")
    for name, rec in o["procedures"].items():
        c = "/".join(str(rec["certified"][str(q)]) for q in QS)
        em = f"{rec['e_max']:.1f}" if "e_max" in rec else "-"
        ks = ("/".join(str(rec["k_star"][str(q)]) for q in QS) if "k_star" in rec
              else "/".join(str(o["k_star_BY"][str(q)]) for q in QS) if name == "BY" else "-")
        print(f"  {name:22s}{em:>9s}   {c:>22s}      {ks}")
    print("-> data/processed/calibrator_comparison.json")
