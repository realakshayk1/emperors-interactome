#!/usr/bin/env python
"""Verify that every numeric claim in the manuscript traces to a committed result artifact.

The paper's central argument is that stated confidence should be checkable. This script
applies that standard to the paper itself: each claim below names the artifact field it
comes from, and the script re-reads that artifact and compares.

This exists because the earlier hand-written fact-check (reports/FACTCHECK.md) compared the
manuscript against *transcribed* artifact values and passed a rounding error as a result.
Re-reading the artifact at check time removes that failure mode.

Usage:
    python scripts/verify_claims.py          # check all claims
    python scripts/verify_claims.py --list   # print the claim -> artifact table
    python scripts/verify_claims.py -v       # show passing claims too

Exit code is 0 when every claim matches, 1 otherwise, so this works as a CI gate or a
pre-release check.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
_CACHE: dict = {}


def load(name: str):
    if name not in _CACHE:
        path = PROCESSED / name
        if not path.exists():
            raise FileNotFoundError(path)
        _CACHE[name] = json.loads(path.read_text(encoding="utf-8"))
    return _CACHE[name]


def dig(obj, path: str):
    """Walk a '/'-separated path. '/' rather than '.' because several artifacts key their
    results by target FDR ('0.05', '0.1', '0.2'), which a dotted separator would split."""
    cur = obj
    for part in path.split("/"):
        if not part:
            continue
        cur = cur[int(part)] if isinstance(cur, list) else cur[part]
    return cur


class Claim:
    def __init__(self, text, file, path, expect, tol=None, transform=None, note=""):
        self.text, self.file, self.path, self.expect = text, file, path, expect
        self.tol = tol if tol is not None else self._auto_tol(expect)
        self.transform, self.note = transform, note

    @staticmethod
    def _auto_tol(expect):
        """Half a unit in the last decimal place the paper actually quotes."""
        if isinstance(expect, int):
            return 0
        s = repr(float(expect))
        return 0.5 * 10 ** (-len(s.split(".")[1])) if "." in s else 0.5

    def check(self):
        actual = dig(load(self.file), self.path)
        if self.transform is not None:
            actual = self.transform(actual)
        if isinstance(self.expect, int) and isinstance(actual, int):
            return actual == self.expect, actual
        return math.isclose(float(actual), float(self.expect), abs_tol=self.tol), actual


CLAIMS = [
    # -- the audit ---------------------------------------------------------
    Claim("candidate edges = 1,666", "audit_summary.json", "n_candidates", 1666),
    Claim("published high-confidence tier = 161", "audit_summary.json", "n_high_conf", 161),
    Claim("certified at q=0.10 = 132", "audit_summary.json", "certified_by_q/0.1", 132),
    Claim("certified at q=0.05 = 78", "audit_summary.json", "certified_by_q/0.05", 78),
    Claim("certified at q=0.20 = 177", "audit_summary.json", "certified_by_q/0.2", 177),
    Claim("high-conf dropped at q=0.10 = 35", "audit_summary.json", "high_conf_dropped_by_q/0.1", 35),
    Claim("high-conf dropped at q=0.05 = 83", "audit_summary.json", "high_conf_dropped_by_q/0.05", 83),
    Claim("high-conf dropped at q=0.20 = 12", "audit_summary.json", "high_conf_dropped_by_q/0.2", 12),
    Claim("dropped fraction = 22%", "audit_summary.json", "high_conf_dropped_by_q/0.1", 0.22,
          tol=0.005, transform=lambda v: v / 161, note="35/161 = 21.7%, quoted as 22%"),

    # -- held-out co-essentiality referee ----------------------------------
    Claim("certified co-essential = 41%", "validation.json", "sets/certified/frac_coess", 0.41),
    Claim("raw high-conf co-essential = 37%", "validation.json", "sets/raw_high_conf/frac_coess", 0.37),
    Claim("dropped co-essential = 17%", "validation.json", "sets/dropped/frac_coess", 0.17),
    Claim("certified vs dropped p = 0.016", "validation.json", "permutation_vs_dropped/p_frac", 0.016),
    Claim("certified vs raw high-conf p = 0.51 (n.s.)", "validation.json", "permutation/p_frac", 0.51),
    Claim("certified edges with referee coverage = 111", "validation.json", "sets/certified/n_scored", 111),
    Claim("dropped edges with referee coverage = 30", "validation.json", "sets/dropped/n_scored", 30),

    # -- the exchangeability violation -------------------------------------
    Claim("measured shift = 0.60 sigma", "sensitivity.json",
          "realized_fdr_population/measured_shift_sigma", 0.60),
    Claim("delta* (control breaks) at q=0.10 = 0.08 sigma", "sensitivity.json",
          "realized_fdr_population/delta_star_control_breaks/0.1", 0.08),
    Claim("realized FDR, identifying curve, q=0.10 = 0.295", "sensitivity.json",
          "realized_fdr_population/identifying_curve_fdr/0.1", 0.295),
    Claim("realized FDR, node-disjoint empirical, q=0.10 = 0.319", "sensitivity.json",
          "realized_fdr_population/empirical_node_disjoint_fdr/0.1", 0.319),

    # -- conditional robustness --------------------------------------------
    Claim("Gamma* at q=0.05 = 29.0", "sensitivity.json",
          "analytic_tilt/gamma_star_conditional/0.05/gamma_star", 29.0),
    Claim("Gamma* at q=0.10 = 31.0", "sensitivity.json",
          "analytic_tilt/gamma_star_conditional/0.1/gamma_star", 31.0),
    Claim("Gamma* at q=0.20 = 31.0", "sensitivity.json",
          "analytic_tilt/gamma_star_conditional/0.2/gamma_star", 31.0),
    Claim("measured-shift equivalent Gamma = 1.82", "sensitivity.json",
          "realized_fdr_population/equivalent_gamma", 1.82),
    Claim("worst-case FDR at measured Gamma, q=0.10 = 0.006", "sensitivity.json",
          "analytic_tilt/fdr_at_measured_gamma/0.1", 0.006),

    # -- attempted repairs -------------------------------------------------
    Claim("hard-negative null KS = 0.156", "hard_negatives.json", "ks_hard_vs_wild/ks", 0.156),
    Claim("plain null KS = 0.129", "hard_negatives.json", "ks_all_vs_wild/ks", 0.129),
    Claim("native decoys shipped = 1,788", "hard_negatives.json", "n_decoy_total", 1788),
    Claim("hard decoys = 1,058", "hard_negatives.json", "n_hard", 1058),
    Claim("soft decoys = 730", "hard_negatives.json", "n_soft", 730),
    Claim("Mondrian pooled held-out FDR = 0.198", "mondrian_summary.json",
          "pooled_overall_heldout_fdr", 0.198),

    # -- dependence --------------------------------------------------------
    Claim("calibration negatives n_cal = 904", "dependence_robustness.json", "n_cal_neg", 904),
    Claim("harmonic penalty H_m = 8.0", "dependence_robustness.json",
          "diagnostic/harmonic_penalty_H_m", 8.0, tol=0.05),
    Claim("conformal p-value floor = 1/905", "dependence_robustness.json",
          "diagnostic/min_possible_p", 1 / 905, tol=1e-9),
    Claim("candidates attaining the floor = 33", "dependence_robustness.json",
          "diagnostic/n_tied_at_min_p", 33),

    # -- orthogonal physical channel ---------------------------------------
    Claim("IntAct certified with evidence = 101", "experimental_ppi_referee.json",
          "results/certified/n_phys", 101),
    Claim("IntAct certified rate = 77%", "experimental_ppi_referee.json",
          "results/certified/obs_rate", 0.77),
    Claim("IntAct degree-matched null = 34%", "experimental_ppi_referee.json",
          "results/certified/matched_null_mean", 0.34),
    Claim("IntAct matched-null OR = 6.2", "experimental_ppi_referee.json",
          "results/certified/odds_ratio_vs_matched_null", 6.2),
    Claim("IntAct naive OR vs background = 21.9", "experimental_ppi_referee.json",
          "results/certified/odds_ratio_vs_background", 21.9),
    Claim("IntAct background rate = 13%", "experimental_ppi_referee.json", "background_rate", 0.13),
    Claim("IntAct dropped rate = 40%", "experimental_ppi_referee.json", "results/dropped/obs_rate", 0.40),
    Claim("IntAct dropped matched null = 33%", "experimental_ppi_referee.json",
          "results/dropped/matched_null_mean", 0.33),
    Claim("IntAct dropped p = 0.85 (n.s.)", "experimental_ppi_referee.json",
          "results/dropped/perm_p", 0.85),

    # -- held-out member recovery ------------------------------------------
    Claim("LOCO complexes = 73", "loco_validation.json", "n_complexes", 73),
    Claim("LOCO member trials = 216", "loco_validation.json", "n_member_trials", 216),
    Claim("LOCO impostor trials = 276", "loco_validation.json", "n_impostor_trials", 276),
    Claim("LOCO member recovery = 49.5%", "loco_validation.json", "observed_member_recovery", 0.495),
    Claim("LOCO impostor recovery = 23.2%", "loco_validation.json", "observed_impostor_recovery", 0.232),
    Claim("LOCO odds ratio = 3.3", "loco_validation.json", "odds_ratio_member_vs_impostor", 3.3),

    # -- multiplicity: the feasibility bounds and Table 1 of the e-values paper ---
    # These are the claims a reviewer is most likely to dispute, and until now none of
    # them were checked. calibrator_comparison.json is regenerated by the same code path
    # that produces the table.
    Claim("tests m = 1,666", "calibrator_comparison.json", "m", 1666),
    Claim("calibration set + 1 = 905", "calibrator_comparison.json", "n_cal_plus_1", 905),
    Claim("H_m = 7.996", "calibrator_comparison.json", "H_m", 7.996),
    Claim("H_(n_cal+1) = 7.386", "calibrator_comparison.json", "H_n_cal_plus_1", 7.386),
    Claim("candidates attaining the floor = 33", "calibrator_comparison.json",
          "n_at_rank/1", 33),
    Claim("candidates with calibration rank <= 2 = 78", "calibrator_comparison.json",
          "n_at_rank/2", 78),

    Claim("k*_BY at q=0.05 = 295", "calibrator_comparison.json", "k_star_BY/0.05", 295),
    Claim("k*_BY at q=0.10 = 148", "calibrator_comparison.json", "k_star_BY/0.1", 148),
    Claim("k*_BY at q=0.20 = 74", "calibrator_comparison.json", "k_star_BY/0.2", 74),

    Claim("BH certifies 78 at q=0.05", "calibrator_comparison.json",
          "procedures/BH/certified/0.05", 78),
    Claim("BH certifies 132 at q=0.10", "calibrator_comparison.json",
          "procedures/BH/certified/0.1", 132),
    Claim("BY certifies 0 at q=0.05", "calibrator_comparison.json",
          "procedures/BY/certified/0.05", 0),
    Claim("BY certifies 0 at q=0.10", "calibrator_comparison.json",
          "procedures/BY/certified/0.1", 0),

    Claim("harmonic calibrator e_max = 122.5", "calibrator_comparison.json",
          "procedures/e-BH harmonic/e_max", 122.5, tol=0.05),
    Claim("harmonic e-BH certifies 0 at q=0.10", "calibrator_comparison.json",
          "procedures/e-BH harmonic/certified/0.1", 0),
    Claim("k*_eBH harmonic at q=0.10 = 136", "calibrator_comparison.json",
          "procedures/e-BH harmonic/k_star/0.1", 136),

    Claim("threshold t=1 e_max = 905", "calibrator_comparison.json",
          "procedures/e-BH threshold_t1/e_max", 905.0, tol=0.05),
    Claim("threshold t=1 e-BH certifies 33 at q=0.10", "calibrator_comparison.json",
          "procedures/e-BH threshold_t1/certified/0.1", 33,
          note="the claim that the harmonic zero is a calibrator artifact rests on this"),

    Claim("threshold t=2 e_max = 452.5", "calibrator_comparison.json",
          "procedures/e-BH threshold_t2/e_max", 452.5, tol=0.05),
    Claim("threshold t=2 e-BH certifies 78 at q=0.05", "calibrator_comparison.json",
          "procedures/e-BH threshold_t2/certified/0.05", 78,
          note="equals BH's count at the same level, without PRDS"),
    Claim("k*_eBH threshold t=2 at q=0.05 = 74", "calibrator_comparison.json",
          "procedures/e-BH threshold_t2/k_star/0.05", 74),

    # -- class optimality: the paper's headline -----------------------------------
    # Swept over every down-set t <= n_cal+1, which by Proposition 2 exhausts the
    # monotone deterministic calibrators.
    Claim("class maximum at q=0.05 = 78", "calibrator_comparison.json",
          "class_optimum/max_certified_over_monotone_calibrators/0.05", 78,
          note="equals BH's count; e-BH cannot beat BH on this map"),
    Claim("class maximum at q=0.10 = 132", "calibrator_comparison.json",
          "class_optimum/max_certified_over_monotone_calibrators/0.1", 132),
    Claim("class maximum at q=0.20 = 177", "calibrator_comparison.json",
          "class_optimum/max_certified_over_monotone_calibrators/0.2", 177),
    Claim("class maximum equals BH at q=0.05", "calibrator_comparison.json",
          "class_optimum/equals_BH/0.05", 1, transform=lambda b: int(bool(b))),
    Claim("class maximum equals BH at q=0.10", "calibrator_comparison.json",
          "class_optimum/equals_BH/0.1", 1, transform=lambda b: int(bool(b))),
    Claim("class maximum equals BH at q=0.20", "calibrator_comparison.json",
          "class_optimum/equals_BH/0.2", 1, transform=lambda b: int(bool(b))),
    Claim("exactly one feasible calibrator at q=0.05", "calibrator_comparison.json",
          "class_optimum/n_feasible_t/0.05", 1),
    Claim("the unique feasible t at q=0.05 is 2", "calibrator_comparison.json",
          "class_optimum/feasible_t/0.05", 2, transform=lambda v: v[0]),
    Claim("seven feasible calibrators at q=0.10", "calibrator_comparison.json",
          "class_optimum/n_feasible_t/0.1", 7),

    # -- the two design levers, and the coincidence we flag as a caution ---------
    Claim("BH's largest accepted rank at q=0.05 = 2", "calibrator_comparison.json",
          "bh_cutoff/0.05/largest_rejected_rank", 2,
          note="equals the threshold t that reproduces BH's count; the paper flags this"),
    Claim("BH's largest accepted rank at q=0.10 = 7", "calibrator_comparison.json",
          "bh_cutoff/0.1/largest_rejected_rank", 7),
    Claim("BH's largest accepted rank at q=0.20 = 19", "calibrator_comparison.json",
          "bh_cutoff/0.2/largest_rejected_rank", 19),
    Claim("exactly one t in 1..25 certifies anything at q=0.05", "calibrator_comparison.json",
          "n_t_certifying_at_0.05", 1),
    Claim("t=1 certifies 0 at q=0.05", "calibrator_comparison.json",
          "t_sweep/0/certified/0.05", 0),
    Claim("t=1 certifies 33 at q=0.10", "calibrator_comparison.json",
          "t_sweep/0/certified/0.1", 33),
    Claim("t=2 certifies 78 at q=0.05", "calibrator_comparison.json",
          "t_sweep/1/certified/0.05", 78),
    Claim("t=2 k* at q=0.05 = 74", "calibrator_comparison.json",
          "t_sweep/1/k_star/0.05", 74),

    # -- the other lever: restricting the hypothesis set to the published tier ----
    Claim("restricted hypothesis set m = 161", "calibrator_comparison.json",
          "restricted_to_tier/m", 161),
    Claim("restricted k*_BY at q=0.05 = 21", "calibrator_comparison.json",
          "restricted_to_tier/k_star_BY/0.05", 21),
    Claim("BY on the published tier certifies 124 at q=0.05", "calibrator_comparison.json",
          "restricted_to_tier/BY/0.05", 124,
          note="refutes any claim that BY admits no remedy on this map"),
    Claim("BY on the published tier certifies 133 at q=0.10", "calibrator_comparison.json",
          "restricted_to_tier/BY/0.1", 133),

    # -- Section 3 figures the paper now quotes precisely -------------------------
    Claim("delta* at q=0.10 = 0.077 sigma", "sensitivity.json",
          "realized_fdr_population/delta_star_control_breaks/0.1", 0.077),
    Claim("measured shift is 7.8x delta*", "sensitivity.json",
          "realized_fdr_population", 7.8, tol=0.05,
          transform=lambda d: d["measured_shift_sigma"] / d["delta_star_control_breaks"]["0.1"]),
    Claim("realized FDR at nominal q=0.05 = 0.181", "sensitivity.json",
          "realized_fdr_population/identifying_curve_fdr/0.05", 0.181),

    # -- shift attribution: Table 1 and the density-ratio diagnosis --------
    # These are the paper's central diagnosis (which axis the shift lives on, and how
    # much of it a 1-D score reweighting can see). They were previously unverified.
    Claim("KS endpoint degree = 0.256", "shift_attribution.json", "per_covariate_ks/degree/ks", 0.256),
    Claim("KS ipTM = 0.132", "shift_attribution.json", "per_covariate_ks/iptm_mean/ks", 0.132),
    Claim("KS confidence score = 0.129", "shift_attribution.json", "per_covariate_ks/score/ks", 0.129),
    Claim("KS ipTM-pTM gap = 0.120", "shift_attribution.json", "per_covariate_ks/iptm_ptm_gap/ks", 0.120,
          tol=0.0005),
    Claim("KS pTM = 0.061", "shift_attribution.json", "per_covariate_ks/ptm_mean/ks", 0.061),
    Claim("density-ratio AUC, full covariate = 0.64", "shift_attribution.json",
          "density_ratio/full_covariate/auc", 0.64),
    Claim("density-ratio AUC, score only = 0.55", "shift_attribution.json",
          "density_ratio/score_only/auc", 0.55,
          note="0.5548 rounds to 0.55; the manuscript once said 0.56 via a double rounding"),
    Claim("divergence invisible to score reweighting = 62%", "shift_attribution.json",
          "density_ratio/divergence_invisible_to_score_frac", 0.62, tol=0.005),

    # -- semi-synthetic benchmark: the prevalence wedge ---------------------
    Claim("benchmark seeds = 100", "benchmark_synth.json", "n_seeds", 100),
    Claim("benchmark-cutoff FDR at pi=0.30, q=0.10 = 0.20", "benchmark_synth.json",
          "rows/1/benchmark_fdr", 0.20, tol=0.005),
    Claim("benchmark-cutoff FDR at pi=0.02, q=0.10 = 0.84", "benchmark_synth.json",
          "rows/10/benchmark_fdr", 0.84),
    Claim("benchmark-cutoff FDR at pi=0.02, q=0.20 = 0.92", "benchmark_synth.json",
          "rows/11/benchmark_fdr", 0.92),
    Claim("conformal holds FDR<=q in all twelve prevalence x q cells", "benchmark_synth.json",
          "rows", 12, transform=lambda rows: sum(r["conformal_fdr"] <= r["q"] for r in rows),
          note="the paper's claim that the same code controls FDR where exchangeability holds"),

    # -- identifying experiment: the exchangeable-null control --------------
    # Appendix A.1 rests the "our implementation is correct" defense on these two.
    Claim("identifying experiment splits = 400", "identifying_experiment.json", "n_splits", 400),
    Claim("realized FDR at delta=0, q=0.10 = 0.0798", "identifying_experiment.json",
          "rows/0/fdr@0.1", 0.0798),

    # -- second real map ---------------------------------------------------
    Claim("Predictomes proteins = 20,196", "secondmap_audit.json", "n_proteins", 20196),
    Claim("Predictomes pairs = 1,614,047", "secondmap_audit.json", "n_pairs", 1614047),
    Claim("Predictomes SPOC>=0.9 tier = 12,767", "secondmap_audit.json", "n_high_conf", 12767),
    Claim("Predictomes certified at q=0.10 = 12,420", "secondmap_audit.json", "certified_by_q/0.1", 12420),
    Claim("Predictomes dropped at q=0.10 = 347", "secondmap_audit.json", "dropped_by_q/0.1", 347),
    Claim("Predictomes certified co-essential = 0.43", "secondmap_audit.json",
          "referee/certified_coess", 0.43),
    Claim("Predictomes dropped co-essential = 0.50", "secondmap_audit.json",
          "referee/dropped_coess", 0.50),
    Claim("Predictomes referee coverage = 320 edges", "secondmap_audit.json", "referee/n_covered", 320),
    Claim("Predictomes dropped edges with referee coverage = 12", "secondmap_audit.json",
          "referee/n_dropped", 12, note="the 0.50 figure rests on n=12; the paper must state this"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="print the claim table and exit")
    ap.add_argument("-v", "--verbose", action="store_true", help="show passing claims too")
    args = ap.parse_args()

    if args.list:
        for c in CLAIMS:
            print(f"{c.text:58s} <- {c.file}:{c.path}")
        return 0

    failures, missing = [], []
    for c in CLAIMS:
        try:
            ok, actual = c.check()
        except FileNotFoundError as e:
            missing.append((c, f"artifact not found: {Path(str(e)).name}"))
            continue
        except (KeyError, IndexError, TypeError) as e:
            missing.append((c, f"field not found ({c.path}): {e!r}"))
            continue
        if ok:
            if args.verbose:
                print(f"  ok   {c.text:58s} = {actual}")
        else:
            failures.append((c, actual))

    for c, actual in failures:
        print(f"  FAIL {c.text}")
        print(f"       claimed {c.expect}  but  {c.file}:{c.path} = {actual}")
        if c.note:
            print(f"       note: {c.note}")
    for c, why in missing:
        print(f"  MISS {c.text}\n       {why}")

    checked = len(CLAIMS) - len(missing)
    print(f"\n{checked - len(failures)}/{checked} claims verified "
          f"({len(failures)} mismatched, {len(missing)} unresolvable) "
          f"across {len({c.file for c in CLAIMS})} artifacts.")
    return 1 if (failures or missing) else 0


if __name__ == "__main__":
    sys.exit(main())
