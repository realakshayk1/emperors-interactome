"""calibrate.py — Day-1 calibration pre-check (METHODS §7, SPEC G1/A1).

Isotonic-fit raw AF-M `score` -> CORUM label on the calibration split; evaluate
Expected Calibration Error (ECE) and reliability on the held-out test split.
A raw metric with ECE_raw >> ECE_isotonic is miscalibrated -> artifacts plausibly
exist (Branch A). Writes reliability data to data/processed/calibration.json.
Run: python -m emperor.calibrate
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score

from . import config as C


def ece(y: np.ndarray, p: np.ndarray, nbins: int = 10) -> float:
    bins = np.linspace(0, 1, nbins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, nbins - 1)
    e, n = 0.0, len(y)
    for b in range(nbins):
        m = idx == b
        if m.sum():
            e += m.sum() / n * abs(y[m].mean() - p[m].mean())
    return float(e)


def reliability_bins(y, p, nbins=10):
    bins = np.linspace(0, 1, nbins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, nbins - 1)
    rows = []
    for b in range(nbins):
        m = idx == b
        if m.sum():
            rows.append(dict(bin_mid=(bins[b] + bins[b + 1]) / 2,
                             conf=float(p[m].mean()), acc=float(y[m].mean()),
                             n=int(m.sum())))
    return rows


def run():
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    cal, test = lab[lab.split == "cal"], lab[lab.split == "test"]
    yte = test.label.values.astype(float)
    ste = test.score.values.astype(float)

    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(cal.score.values, cal.label.values)
    pte = iso.predict(ste)

    res = dict(
        n_test=int(len(test)), n_pos=int(yte.sum()),
        auroc=float(roc_auc_score(yte, ste)),
        ece_raw=ece(yte, ste), ece_isotonic=ece(yte, pte),
        reliability_raw=reliability_bins(yte, ste),
        reliability_isotonic=reliability_bins(yte, pte),
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "calibration.json").write_text(json.dumps(res, indent=2))
    return res


def main():
    r = run()
    print(f"test n={r['n_test']} pos={r['n_pos']}  AUROC={r['auroc']:.3f}")
    print(f"ECE raw={r['ece_raw']:.3f}  ECE isotonic={r['ece_isotonic']:.3f}")
    branch = "A (artifacts likely)" if r["ece_raw"] > 2 * r["ece_isotonic"] else "B (well-calibrated)"
    print(f"pre-check -> Branch {branch}")


if __name__ == "__main__":
    main()
