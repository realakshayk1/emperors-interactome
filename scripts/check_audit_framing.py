"""Check the framing claims the papers make about the audit, which are not single artifact
lookups and so are invisible to verify_claims.py.

Each of these was a defect found by review rather than by the checker, and each is the kind of
statement that reproduces perfectly while describing the wrong population or the wrong family:

  1. BH on the conformal p-values is a score threshold. The p-value is monotone in the released
     score, so the certified set is exactly the top-R candidates by score and the audit
     contributes the cut point, not the ranking.
  2. The "22% of the tier fails" headline is a property of the testing family. No tier edge has
     a p-value above 0.041, so the 161 tested as their own family certify entirely.
  3. The 0.60 sigma calibration-to-candidate shift is dominated by the tier under audit; the
     residual over the remaining candidates is below the measured breaking point.
  4. Randomised calibrators break the class identity and beat BH in expectation, which is
     why Proposition 3 is stated for deterministic ones. A sweep over a few fixed firing
     probabilities suggested the opposite; optimising the probability does not.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
QS = (0.05, 0.10, 0.20)


def _load():
    d = pd.read_parquet(PROC / "certified.parquet")
    dep = json.loads((PROC / "dependence_robustness.json").read_text())
    return d, dep["n_cal_neg"] + 1


def _bh(pv: np.ndarray, q: float) -> int:
    m = pv.size
    s = np.sort(pv)
    ok = np.nonzero(s <= np.arange(1, m + 1) * q / m)[0]
    return int(ok[-1] + 1) if ok.size else 0


def _ebh(e: np.ndarray, q: float) -> int:
    m = e.size
    s = np.sort(e)[::-1]
    ok = np.nonzero(s >= m / (np.arange(1, m + 1) * q))[0]
    if not ok.size:
        return 0
    k = int(ok[-1] + 1)
    return int((e >= m / (k * q)).sum())


def main() -> int:
    d, n1 = _load()
    p = d["conf_pvalue"].to_numpy(float)
    score = d["score"].to_numpy(float)
    hi = d["high_conf"].to_numpy(bool)
    m = p.size
    rank = np.rint(p * n1).astype(int)
    fails = []

    # 1 -- the audit is a threshold on the released score
    order = np.argsort(score)
    if not np.all(np.diff(p[order]) <= 1e-12):
        fails.append("conformal p-value is not monotone in score")
    for q in QS:
        r = _bh(p, q)
        top = np.zeros(m, bool)
        top[np.argsort(-score)[:r]] = True
        cert = d[f"certified@{q if q != 0.10 else 0.1}"].to_numpy(bool)
        if not (cert == top).all():
            fails.append(f"certified@{q} is not the top-{r} by score")
    print(f"  certified set equals the top-R by score at every q  ({', '.join(str(_bh(p, q)) for q in QS)})")

    # 2 -- the headline depends on the testing family
    print(f"  max conformal p over the {int(hi.sum())} tier edges = {p[hi].max():.4f}")
    for q in QS:
        tier_alone = _bh(p[hi], q)
        if tier_alone != int(hi.sum()):
            fails.append(f"tier-as-own-family at q={q} certifies {tier_alone}, expected all")
        full = _bh(p, q)
        thr = np.sort(p)[full - 1]
        dropped = int(hi.sum()) - int((p[hi] <= thr).sum())
        print(f"    q={q}: family m={m} drops {dropped} of {int(hi.sum())}; family m={int(hi.sum())} drops 0")

    # 3 -- the shift is dominated by the tier
    w = json.loads((PROC / "wcs_results.json").read_text())["shift"]
    mu_cal = w["score_mean_calneg"]
    # Read the sd rather than back-deriving it from the rounded 0.600 we are checking;
    # dividing by the answer makes the assertion below unfalsifiable.
    sd = json.loads((PROC / "gamma_seed.json").read_text())["nonconf_sd"]
    gap_all = score.mean() - mu_cal
    gap_rest = score[~hi].mean() - mu_cal
    # The gap decomposes by pool share, which is what "92% comes from the tier" means; the
    # unweighted ratio of the residual to the total decomposes nothing.
    n_hi = int(hi.sum())
    share = (n_hi / m) * (score[hi].mean() - mu_cal) / gap_all
    rest_share = ((m - n_hi) / m) * (gap_rest) / gap_all
    print(f"  shift {gap_all / sd:.3f} sd total; tier contributes {100 * share:.1f}%, "
          f"the other {m - n_hi} contribute {100 * rest_share:.1f}%; "
          f"residual gap {gap_rest / sd:.3f} sd")
    # share + rest_share == 1 is algebra, not a test. The claim with content is that the tier
    # contributes the reported fraction while holding a tenth of the candidates.
    if not (0.09 <= n_hi / m <= 0.10):
        fails.append(f"tier is {n_hi}/{m} of candidates, not the ~10% the decomposition assumes")
    if not (0.915 <= share <= 0.925):
        fails.append(f"tier share of the shift is {share:.4f}, expected ~0.92")
    if abs(gap_all / sd - 0.600) > 0.005:
        fails.append(f"total shift is {gap_all / sd:.3f} sd, expected 0.600")

    # the identifying curve evaluated at the residual rather than the full gap
    rows = json.loads((PROC / "identifying_experiment.json").read_text())["rows"]
    lo, hi_f = rows[0]["fdr@0.1"], rows[1]["fdr@0.1"]
    at_res = lo + (hi_f - lo) * ((gap_rest / sd) / rows[1]["delta"])
    print(f"  identifying curve at the residual shift: FDR {at_res:.3f} (nominal 0.10)")
    if at_res >= 0.10:
        fails.append(f"residual-shift FDR is {at_res:.3f}, not below nominal")

    # 4 -- randomised calibrators escape the deterministic ceiling
    # The calibrator is valid for EVERY firing probability, since pi cancels in the
    # expectation. pi <= q*N(t)*n1/(m*t) is the separate condition for e-BH to take all of
    # N(t) when it fires; expected yield is then pi*N(t). Optimise over t rather than fixing
    # a few c, since a narrow sweep is how the opposite (and false) claim survived here.
    print("  randomised calibrator, firing probability optimised:")
    for q in QS:
        bh_q = _bh(p, q)
        best, arg = 0.0, None
        for t in range(1, n1 + 1):
            n_t = int((rank <= t).sum())
            if not n_t:
                continue
            pi = min(1.0, q * n_t * n1 / (m * t))
            if pi * n_t > best:
                best, arg = pi * n_t, (t, pi)
        print(f"    q={q}: expected {best:.1f} at t={arg[0]} (pi={arg[1]:.4f}) vs BH {bh_q}")
        if best <= bh_q:
            fails.append(f"randomised optimum {best:.1f} at q={q} does not exceed BH {bh_q}")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("all audit-framing claims hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
