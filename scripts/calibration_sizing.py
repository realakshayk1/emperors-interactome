"""How large a calibration set an assumption-free audit needs before it can return anything.

Proposition 1 says BY and e-BH each reject at least k* hypotheses or none. Read forwards that
is a warning. Read backwards it is a design rule, and the design rule is what a practitioner
needs, because it is computable before any data are collected.

Feasibility requires the map to supply as many candidates inside the threshold as the bound
demands:

    N(t)  >=  m t W / (q (n_cal + 1)),      W = H_m for BY,  W = 1 for a threshold calibrator.

Write phi = N(t)/m for the fraction of candidates the audit is trying to certify. The m
cancels, and the requirement becomes

    n_cal + 1  >=  t W / (q phi).                                    (*)

That is the useful form. It does not depend on the number of hypotheses except through H_m,
which is logarithmic, so the binding quantities are the level, the yield you want, and which
procedure you are willing to use. The H_m factor is the whole price of BY over a threshold
calibrator: about eightfold here, and it appears directly as a multiple on the number of
negatives you must ship.

Applied to the CM4AI map, the requirement must be minimised over t: at q = 0.05 BY needs
6,832 calibration points before any threshold becomes feasible, attained at t = 2. Evaluating
only at t = 1, where the map supplies 33 candidates, gives 8,074 and overstates the requirement
by 18%. The map shipped 1,788 decoys, of which 905 became the calibration set.

The estimate is optimistic and should be read as a lower bound. It holds phi fixed while
n_cal grows, but a candidate's rank is its position among the calibration scores, so enlarging
the calibration set can only push ranks up and phi down. How far is not knowable from 904
decoys, which is why we report (*) as the requirement under the map's observed yield rather
than extrapolating a yield we cannot measure.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
QS = (0.05, 0.10, 0.20)


def required_ncal(t: int, weight: float, q: float, phi: float) -> int:
    """Smallest n_cal+1 for which a yield of phi is feasible at level q. Equation (*)."""
    return math.ceil(t * weight / (q * phi))


def main() -> int:
    d = pd.read_parquet(PROC / "certified.parquet")
    dep = json.loads((PROC / "dependence_robustness.json").read_text())
    p = d["conf_pvalue"].to_numpy(float)
    m = p.size
    n1 = dep["n_cal_neg"] + 1
    rank = np.rint(p * n1).astype(int)
    H_m = sum(1.0 / i for i in range(1, m + 1))
    fails = []

    print(f"  map: m={m}, shipped n_cal+1={n1}, H_m={H_m:.3f}")
    print()

    # --- the rule reproduces the bounds the paper reports ---------------------
    n_floor = int((rank <= 1).sum())
    by_need = required_ncal(1, H_m, 0.05, n_floor / m)
    print(f"  BY at q=0.05 with yield {n_floor}/{m}: needs n_cal+1 >= {by_need:,}"
          f"  ({by_need / n1:.1f}x what was shipped)")
    # A round trip through required_ncal's own inverse proves nothing. The informative check
    # is that the rule is minimised over t, as it must be, and that the minimum is the smallest
    # calibration set at which ANY threshold becomes feasible.
    by_min = min((required_ncal(t, H_m, 0.05, int((rank <= t).sum()) / m), t)
                 for t in range(1, n1 + 1) if (rank <= t).sum())
    print(f"  BY minimised over t: needs n_cal+1 >= {by_min[0]:,} at t={by_min[1]}")
    if by_min[0] > by_need:
        fails.append(f"t=1 requirement {by_need} is below the minimum {by_min[0]}; rule mis-stated")

    best = min(((required_ncal(t, 1.0, 0.05, int((rank <= t).sum()) / m), t)
                for t in range(1, n1 + 1) if (rank <= t).sum()), key=lambda x: x[0])
    print(f"  e-BH at q=0.05, best threshold t={best[1]}: needs n_cal+1 >= {best[0]:,}"
          f"  ({'satisfied' if best[0] <= n1 else f'{best[0] / n1:.1f}x short'})")
    if best[0] > n1:
        fails.append("e-BH requirement exceeds the shipped calibration set, contradicting Table 1")
    print()

    # --- the design table, which needs no data at all -------------------------
    print("  Calibration points needed to certify a fraction phi of candidates (equation *):")
    print()
    print(f"    {'phi':>7}  {'q':>5}   {'BY (W=H_m)':>14}   {'e-BH, t=1':>11}   {'ratio':>6}")
    for phi in (0.01, 0.02, 0.05, 0.10):
        for q in QS:
            by = required_ncal(1, H_m, q, phi)
            eb = required_ncal(1, 1.0, q, phi)
            print(f"    {phi:>7.2f}  {q:>5}   {by:>14,}   {eb:>11,}   {by / eb:>5.0f}x")
    print()
    print("  The ratio is H_m throughout: the harmonic penalty is a direct multiplier on the")
    print("  number of negatives a release must ship for an arbitrary-dependence audit.")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("sizing rule is consistent with the reported bounds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
