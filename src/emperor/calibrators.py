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
T_SWEEP = range(1, 26)   # report the whole grid, not whichever cell flatters the result


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

    # --- the whole threshold grid, and where BH itself stopped --------------------
    # Threshold-t e-BH rejects exactly {R_j <= t} when it rejects at all, and BH rejects a
    # rank prefix, so setting t to BH's realized cutoff makes the two sets identical by
    # construction. Recording the cutoff makes that coincidence visible instead of
    # letting a tuned t look like a discovery.
    out["bh_cutoff"] = {}
    for q in QS:
        r = benjamini_hochberg(p, q)
        cut = int(rank[r].max()) if r.any() else None
        out["bh_cutoff"][str(q)] = {"largest_rejected_rank": cut,
                                    "n_at_or_below": int((rank <= cut).sum()) if cut else 0}

    # --- the whole admissible calibrator class -----------------------------------
    # A calibrator maps rank R to an e-value e(R). Exchangeability constrains the rank only
    # through P(R <= t) <= t/(n_cal+1), under which
    #     sup E[e(R)] = (1/(n_cal+1)) * sum_r max_{s>=r} e(s).
    # For NON-INCREASING e that supremum is the budget sum_r e(r) <= n_cal+1, so the budget
    # is validity. Off the monotone class it is not: a calibrator can spend exactly the
    # budget and still reach E[e] > 1, so the budget alone must never be used as an
    # admissibility test (scripts/check_calibrator_validity.py exhibits one at 15/7).
    # Restricting to non-increasing e costs nothing, because the right-running maximum is
    # valid whenever e is, dominates it pointwise, and e-BH only grows in its e-values.
    # Monotonicity then makes e-BH reject a down-set {R <= t}, with t*e(t) <= n_cal+1 and so
    # e_max <= (n_cal+1)/t, and sweeping t over 1..n_cal+1 exhausts the class.
    best, feasible = {}, {}
    for q in QS:
        b, feas = 0, []
        for t in range(1, n1 + 1):
            e_max_t = n1 / t
            if int((rank <= t).sum()) >= math.ceil(m / (q * e_max_t)):
                feas.append(t)
                b = max(b, int(ebh(threshold_evalues(rank, n_cal, t), q).sum()))
        best[str(q)] = b
        feasible[str(q)] = feas
    out["class_optimum"] = {
        "max_certified_over_monotone_calibrators": best,
        "feasible_t": {q: v if len(v) <= 20 else v[:20] for q, v in feasible.items()},
        "n_feasible_t": {q: len(v) for q, v in feasible.items()},
        "equals_BH": {q: best[q] == out["procedures"]["BH"]["certified"][q] for q in best},
    }

    out["t_sweep"] = []
    for t in T_SWEEP:
        e = threshold_evalues(rank, n_cal, t)
        row = {"t": t, "n_at_or_below": int((rank <= t).sum()),
               "e_max": float(n1 / t),
               "k_star": {str(q): math.ceil(m / (q * (n1 / t))) for q in QS},
               "certified": {str(q): int(ebh(e, q).sum()) for q in QS}}
        out["t_sweep"].append(row)
    out["n_t_certifying_at_0.05"] = sum(1 for r in out["t_sweep"]
                                        if r["certified"]["0.05"] > 0)

    # --- the other lever: the size of the hypothesis set --------------------------
    # k*_BY = ceil(m H_m F / q) grows with m, so restricting the audit to the published
    # tier is a design choice that moves BY's bound as surely as the calibrator moves
    # e-BH's. The tier is defined by the map's authors, independently of any p-value.
    p_hi = p[hi]
    out["restricted_to_tier"] = {
        "m": int(p_hi.size), "H_m": _harmonic(p_hi.size),
        "contains_all_rank1": int((rank[hi] <= 1).sum()) == int((rank <= 1).sum()),
        "contains_all_rank2": int((rank[hi] <= 2).sum()) == int((rank <= 2).sum()),
        "k_star_BY": {str(q): math.ceil(p_hi.size * _harmonic(p_hi.size) / (q * n1))
                      for q in QS},
        "BH": {str(q): int(benjamini_hochberg(p_hi, q).sum()) for q in QS},
        "BY": {str(q): int(benjamini_yekutieli(p_hi, q).sum()) for q in QS},
    }

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
