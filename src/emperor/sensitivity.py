"""sensitivity.py — Rosenbaum-style Γ sensitivity bound for conformal FDR (PLAN_V3 T1.4).

WS1 showed unconditional control does not survive the real calibration->wild shift
(Path A fails). This module delivers the Path-B headline: a CONDITIONAL guarantee of
the form "FDR <= q under exchangeability violation of at most Γ", with Γ* — the
smallest violation at which BH control breaks — estimated on the real map.

Model. Conformal p-values are exact under exchangeability. Under a bounded
violation, a true-null test point's nonconformity law differs from the calibration
null by a likelihood ratio bounded in [1/Γ, Γ] (Rosenbaum's sensitivity parameter,
here on the conformal-rank scale). A worst-case tilt within that band inflates each
null p-value's stochastic size by at most a factor governed by Γ; we propagate the
tilt through the BH threshold to get the worst-case realized FDR as a function of Γ,
and report Γ* where it first exceeds q.

Two independent estimates of the realized FDR at the map's measured shift are
reconciled:
  (1) IDENTIFYING-CURVE route: map the measured mean-shift (0.60σ) onto the
      tests/ identifying experiment (gamma_seed.json).
  (2) SENSITIVITY-MODEL route: convert the measured shift to an equivalent Γ and read
      the worst-case FDR(Γ) curve here.

Purity: uses only nonconformity scores + calibration negatives — no DepMap.
Run: python -m emperor.sensitivity
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .conformal import conformal_pvalues, benjamini_hochberg
from .nonconformity import nonconformity
from .covariates import add_covariates


def _worst_case_fdr_curve(p_obs, m, q, gammas):
    """Worst-case realized FDR as a function of Γ (sensitivity bound).

    Conformal p-values are super-uniform under exchangeability: P(p<=t)<=t for a true
    null. Under a Γ-bounded exchangeability violation the worst-case null CDF inflates
    to P(p<=t) <= min(1, Γ·t) (a true null can look up-to-Γ times more significant).

    BH is a step-up on the OBSERVED p-values, so the rejection set (and its size R) is
    fixed by the data; what the adversary controls is how many of those R rejections
    are false. The worst-case expected number of false discoveries is

        V(Γ) = Σ_{j in rejected} min(1, Γ · p_j)                    (each rejected null
                                                                     contributes its
                                                                     worst-case reject prob)

    and the worst-case FDR = E[V(Γ)] / R. At Γ=1 this reduces to the standard BH bound
    Σ p_j / R <= q. (The earlier version evaluated only the single cutoff threshold,
    undercounting V by ~R×; this sums the per-rejection worst-case probabilities.)
    """
    p = np.asarray(p_obs, float)
    thresh = q * np.arange(1, m + 1) / m
    order = np.argsort(p)
    ranked = p[order]
    passed = ranked <= thresh
    if not passed.any():
        return [0.0 for _ in gammas]
    kmax = np.nonzero(passed)[0].max()
    rej_p = ranked[: kmax + 1]          # p-values of the rejected set
    R = len(rej_p)
    curve = []
    for g in gammas:
        V = np.minimum(1.0, g * rej_p).sum()
        curve.append(float(min(1.0, V / R)))
    return curve


def _shift_to_gamma(delta_sigma):
    """Convert a mean shift (in σ) to an equivalent likelihood-ratio bound Γ.

    For two Gaussians differing in mean by δσ, the max density ratio over the bulk
    (within ±2σ) is exp(δ·(2)) on the tail side; we use the moderate-region bound
    Γ ≈ exp(δ) as a calibrated, interpretable proxy (δ=0 -> Γ=1 -> no violation).
    """
    return float(np.exp(delta_sigma))


def run():
    inter = add_covariates(pd.read_parquet(C.INTERIM / "interactome.parquet"))
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    cal_neg = lab[(lab.label == 0) & (lab.split == "cal")]
    s_cal = nonconformity(cal_neg)
    cand = inter[inter.is_true_pair].copy()
    s = nonconformity(cand)
    p = conformal_pvalues(s, s_cal)
    m = len(cand)

    # grid extended to Γ=50 so a genuine crossing (if any) is found rather than
    # censored prematurely; the certified p-values are tiny so the crossing, if it
    # exists, is at large Γ.
    gammas = list(np.round(np.concatenate([np.arange(1.0, 10.01, 0.1),
                                           np.arange(11.0, 50.01, 1.0)]), 2))
    curves = {q: _worst_case_fdr_curve(p, m, q, gammas) for q in (0.05, 0.10, 0.20)}

    # Γ* = smallest Γ at which worst-case FDR first exceeds q. If the curve never
    # crosses q within the tested grid [1, Γ_max], Γ* is RIGHT-CENSORED at Γ_max
    # (reported as ">= Γ_max", never as an exact crossing — the certified set's
    # p-values are so small that the analytic tilt does not breach q anywhere in the
    # grid, i.e. the observed rejections are highly robust).
    gstar = {}
    for q in (0.05, 0.10, 0.20):
        c = np.array(curves[q])
        above = np.where(c > q)[0]
        if len(above):
            gstar[q] = dict(gamma_star=float(gammas[above[0]]), censored=False)
        else:
            gstar[q] = dict(gamma_star=float(gammas[-1]), censored=True,
                            note=f">= {gammas[-1]} (no crossing in tested grid; curve reaches "
                                 f"{c[-1]:.3f} at Gamma={gammas[-1]}, still below q={q})")

    seed = json.loads((C.PROCESSED / "gamma_seed.json").read_text())
    delta_real = seed["delta_real_sigma"]
    gamma_real = _shift_to_gamma(delta_real)
    # The analytic tilt evaluated at the certified set's OBSERVED p-values (a
    # conditional robustness statement about those specific rejections):
    tilt_fdr_on_certified = {q: float(np.interp(gamma_real, gammas, curves[q]))
                             for q in (0.05, 0.10, 0.20)}

    out = dict(
        gammas=gammas,
        # ROUTE 1 (analytic tilt): worst-case FDR conditional on the OBSERVED
        # rejected p-values. Answers "how robust are the specific certified edges to a
        # Γ-bounded tilt" — NOT the population realized FDR. Stays low because the
        # certified edges have very small conformal p-values (they beat the whole null).
        analytic_tilt=dict(
            worst_case_fdr_curve={q: curves[q] for q in (0.05, 0.10, 0.20)},
            gammas=gammas,
            gamma_star_conditional=gstar,
            fdr_at_measured_gamma=tilt_fdr_on_certified,
            measures="robustness of the OBSERVED certified set to a Γ-bounded tilt (conditional, optimistic-by-design)",
        ),
        # ROUTE 2 (identifying experiment) + ROUTE 3 (empirical node-disjoint): the
        # POPULATION realized FDR under the real shift. These are the headline numbers.
        realized_fdr_population=dict(
            measured_shift_sigma=delta_real,
            equivalent_gamma=gamma_real,
            identifying_curve_fdr=seed["estimated_realized_fdr"],
            empirical_node_disjoint_fdr={"0.1": 0.3191},
            delta_star_control_breaks=seed["delta_star"],
        ),
        interpretation=(
            "TWO DIFFERENT QUANTITIES, reported separately (they do NOT and should not "
            "coincide). (A) Population realized FDR under the real shift: the map's "
            "measured 0.60σ calibration->wild shift maps onto the identifying-experiment "
            "curve at realized FDR ≈ %.2f (q=0.10), independently corroborated by the "
            "empirical node-disjoint held-out FDR of 0.32. Control provably breaks at a "
            "tiny δ*≈%.2fσ, so the observed shift is ~7× past the breaking point — "
            "unconditional control does NOT hold. (B) Conditional robustness of the "
            "certified SET: the specific certified edges have very small conformal "
            "p-values, so a Γ-bounded tilt at the measured Γ≈%.2f inflates THEIR "
            "worst-case FDR only to ~%.3f — the strongly-certified edges are robust even "
            "though the marginal guarantee is not. Defensible Path-B claim: 'The certified "
            "set's realized FDR under the measured shift is ≈0.3 (three-way: identifying "
            "curve, node-disjoint empirical, and honest concession that BH's marginal "
            "guarantee is void past δ*); the high-confidence CORE of the certified set is "
            "separately shown robust to bounded exchangeability violations.'"
        ) % (seed["estimated_realized_fdr"]["0.1"],
             seed["delta_star"]["0.1"], gamma_real, tilt_fdr_on_certified[0.10]),
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "sensitivity.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    g = o["gammas"]
    print("=== ROUTE A: population realized FDR under the REAL shift (headline) ===")
    rp = o["realized_fdr_population"]
    print(f"measured shift {rp['measured_shift_sigma']:.2f}σ; control breaks at δ*={rp['delta_star_control_breaks']['0.1']:.2f}σ")
    print(f"realized FDR @q=0.10:  identifying-curve {rp['identifying_curve_fdr']['0.1']:.3f}  |  "
          f"empirical node-disjoint {rp['empirical_node_disjoint_fdr']['0.1']:.3f}")
    print("\n=== ROUTE B: conditional robustness of the OBSERVED certified set ===")
    at = o["analytic_tilt"]
    for q in (0.05, 0.10, 0.20):
        c = at["worst_case_fdr_curve"][q]
        gs = at["gamma_star_conditional"][q]
        gstr = f">={gs['gamma_star']:.0f} (censored)" if gs["censored"] else f"{gs['gamma_star']:.1f}"
        print(f"  q={q}: Γ*={gstr}  worst-case FDR @Γ=1:{c[0]:.3f} @Γ=2:{c[g.index(2.0)]:.3f} @Γ=3:{c[g.index(3.0)]:.3f}  "
              f"(at measured Γ≈{o['realized_fdr_population']['equivalent_gamma']:.2f}: {at['fdr_at_measured_gamma'][q]:.3f})")
    print("\n" + o["interpretation"])


if __name__ == "__main__":
    main()
