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
  4. Randomised calibrators break the class identity but lose in expectation, which is why
     Proposition 3 is stated for deterministic ones.
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
    sd = (w["score_mean_wild"] - mu_cal) / 0.600
    gap_all = score.mean() - mu_cal
    gap_rest = score[~hi].mean() - mu_cal
    share = 1 - gap_rest / gap_all
    print(f"  shift {gap_all / sd:.3f} sd total; tier contributes {100 * share:.1f}%; "
          f"residual {gap_rest / sd:.3f} sd")
    if not (0.90 <= share <= 0.92):
        fails.append(f"tier share of the shift is {share:.3f}, expected ~0.91")

    # the identifying curve evaluated at the residual rather than the full gap
    rows = json.loads((PROC / "identifying_experiment.json").read_text())["rows"]
    lo, hi_f = rows[0]["fdr@0.1"], rows[1]["fdr@0.1"]
    at_res = lo + (hi_f - lo) * ((gap_rest / sd) / rows[1]["delta"])
    print(f"  identifying curve at the residual shift: FDR {at_res:.3f} (nominal 0.10)")
    if at_res >= 0.10:
        fails.append(f"residual-shift FDR is {at_res:.3f}, not below nominal")

    # 4 -- randomised calibrators win conditionally and lose in expectation
    print("  randomised calibrator, fires with probability 1/c:")
    for q in (0.05, 0.10):
        bh_q = _bh(p, q)
        for c in (2, 3, 4):
            best = max(_ebh(np.where(rank <= t, c * n1 / t, 0.0), q) for t in range(1, n1 + 1))
            exp = best / c
            flag = "ok" if exp < bh_q else "BEATS BH IN EXPECTATION"
            print(f"    q={q} c={c}: {best} when it fires, {exp:.1f} expected vs BH {bh_q}  {flag}")
            if exp >= bh_q:
                fails.append(f"randomised c={c} at q={q} beats BH in expectation")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("all audit-framing claims hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
