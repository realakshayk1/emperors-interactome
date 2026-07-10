"""
§A — The identifying experiment: exchangeability is the CAUSE of distribution-free control.

The rest of this project shows the real CORUM/DepMap calibration null is NOT
exchangeable with the wild candidate pool, and that under complex/node-disjoint
splits distribution-free FDR control does not hold (the honest negative). That is
an observation. This experiment turns it into a CAUSAL claim by isolating the one
variable responsible.

Design (single-knob toggle). In the synthetic setting where the guarantee is
PROVABLE (the unit-test setting), hold everything fixed — number of true
positives, number of nulls, calibration size, the BH procedure, q, and all
random draws' shapes — and vary ONLY one scalar: delta, the distributional shift
between the calibration-negative law and the test-null law.

    cal_neg  ~ N(0, 1)          (what the auditor calibrates on)
    test_null~ N(-delta, 1)     <- the ONLY thing delta changes: the true null,
                                   shifted toward the interaction-like signal
    test_pos ~ N(-2.5, 1)       (fixed; the interaction-like signal)

delta = 0  => cal and test nulls are exchangeable  => conformal p-values of nulls
             are (super-)uniform => BH controls FDR <= q. CONTROL ON.
delta > 0  => the true null looks MORE interaction-like than the null the auditor
             calibrated on => conformal p-values anti-conservative => realized FDR
             climbs above q. CONTROL OFF. (nonconformity is oriented so that
             lower = more interaction-like, matching the project convention.)

Because delta is the only thing that moves, any change in realized FDR is
attributable to exchangeability alone. This is the identifying experiment:
distribution-free control is present if and only if the null is exchangeable.
"""
import json
import numpy as np
from . import config as C
from .conformal import conformal_pvalues, benjamini_hochberg


def run(n_splits: int = 400, deltas=None, n_pos: int = 200, n_null: int = 800,
        n_cal: int = 1500):
    if deltas is None:
        deltas = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    qs = sorted(set(C.Q_SWEEP) | {C.Q})
    rows = []
    for delta in deltas:
        per_q = {q: [] for q in qs}
        for seed in range(n_splits):
            rng = np.random.default_rng(20_000 + seed)
            cal_neg   = rng.normal(0.0,    1.0, n_cal)    # what the auditor calibrates on
            test_null = rng.normal(-delta, 1.0, n_null)   # true null, shifted TOWARD the signal
            test_pos  = rng.normal(-2.5,  1.0, n_pos)     # fixed signal
            s_test = np.concatenate([test_pos, test_null])
            is_null = np.concatenate([np.zeros(n_pos, bool), np.ones(n_null, bool)])
            p = conformal_pvalues(s_test, cal_neg)
            for q in qs:
                rej = benjamini_hochberg(p, q)
                per_q[q].append(float(is_null[rej].mean()) if rej.sum() else 0.0)
        row = dict(delta=float(delta))
        for q in qs:
            arr = np.array(per_q[q])
            row[f"fdr@{q}"] = float(arr.mean())
            row[f"se@{q}"]  = float(arr.std() / np.sqrt(len(arr)))
        row["control_holds@primary"] = bool(row[f"fdr@{C.Q}"] <= C.Q + 0.01)
        rows.append(row)
    out = dict(
        design="synthetic identifying experiment: toggle only cal/test null exchangeability",
        toggle="delta = mean shift of test-null law vs calibration-null law (in sigma)",
        held_fixed=["n_pos", "n_null", "n_cal", "test_pos law", "cal_neg law",
                    "BH procedure", "q", "conformal p-value construction"],
        q_primary=C.Q, q_sweep=qs, n_splits=n_splits, rows=rows,
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "identifying_experiment.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    o = run()
    qp = o["q_primary"]
    for r in o["rows"]:
        state = "ON" if r["control_holds@primary"] else "OFF"
        print(f"delta={r['delta']:.2f}  FDR@{qp}={r[f'fdr@{qp}']:.3f}  control={state}")
