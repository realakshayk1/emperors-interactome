"""gamma_seed.py — map the measured real-data shift onto the §19 identifying curve
(PLAN_V3 markup #5 seed for the Path-B conditional guarantee).

Reads the extended-program inputs (data/external/{identifying_experiment,wcs_results}.json)
and the calibration labels, expresses the measured calibration->wild mean-shift in
units of the nonconformity SD, interpolates the identifying-experiment δ->realized-FDR
curve at that shift, and finds δ* (the shift at which control breaks). Emits
data/processed/gamma_seed.json, consumed by sensitivity.py.

Purity: uses AF scores + calibration negatives only — no DepMap. Run: python -m emperor.gamma_seed
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C

EXT = C.ROOT / "data" / "external"


def run():
    ident = json.loads((EXT / "identifying_experiment.json").read_text())
    wcs = json.loads((EXT / "wcs_results.json").read_text())["shift"]
    rows = ident["rows"]
    deltas = np.array([r["delta"] for r in rows])

    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    sd_cal = (1.0 - lab[(lab.label == 0) & (lab.split == "cal")]["score"]).std()
    delta_real = abs(wcs["score_mean_wild"] - wcs["score_mean_calneg"]) / sd_cal

    est, dstar = {}, {}
    for q in ["0.05", "0.1", "0.2"]:
        fq = np.array([r[f"fdr@{q}"] for r in rows])
        est[q] = float(np.interp(delta_real, deltas, fq))
        qv = float(q)
        above = np.where(fq > qv)[0]
        if len(above) and above[0] > 0:
            i = above[0]
            dstar[q] = float(deltas[i - 1] + (qv - fq[i - 1]) * (deltas[i] - deltas[i - 1]) /
                             (fq[i] - fq[i - 1]))
        else:
            dstar[q] = 0.0 if (len(above) and above[0] == 0) else float(deltas[-1])

    out = dict(delta_real_sigma=float(delta_real),
               mean_gap=float(abs(wcs["score_mean_wild"] - wcs["score_mean_calneg"])),
               nonconf_sd=float(sd_cal),
               estimated_realized_fdr=est, delta_star=dstar,
               basis="interp on identifying-experiment curve (n_splits=400)",
               shift_source="wcs_results.shift: score_mean_calneg vs score_mean_wild / nonconf sd")
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "gamma_seed.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"measured shift = {o['delta_real_sigma']:.3f} sigma (gap {o['mean_gap']:.4f} / sd {o['nonconf_sd']:.4f})")
    print("estimated realized FDR:", {k: round(v, 3) for k, v in o["estimated_realized_fdr"].items()})
    print("delta*:", {k: round(v, 3) for k, v in o["delta_star"].items()})


if __name__ == "__main__":
    main()
