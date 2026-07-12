"""hardnegatives.py — hard-negative calibration set (PLAN_V3 T1.2, markup #1).

The T1.1 diagnosis: ~62% of the calibration->wild shift is invisible to a 1-D
`score` reweighting, and the dominant unmodeled axis is union-graph DEGREE — the
wild candidates are hub-enriched by AP-MS+AF pre-screening, while random decoys
are not. So the random-decoy null does not span the "selected-looking but
non-interacting" region the candidates occupy; no reweighting of that null can fix
control (need a better null).

Fix: reweight the calibration null toward its HARD region — non-CORUM-co-complex
decoy pairs whose union-graph degree is in the top tertile of the CANDIDATE degree
distribution (high-connectivity pairs that look selected like the candidates but
are not annotated interactors). NOTE (DECISIONS 2026-07-10 deviation): the
pre-registered rule said "augment with newly-drawn pairs," but new pairs are
unfolded (no AF score) and cannot enter a score-based null without GPU folding
(deferred WS2.1). Feasible substitute: upweight the hard region of the EXISTING
(scored) decoy set by an explicit replication factor. Intent — make the null span
the hub-selected region — is preserved using only scored pairs.

PRE-REGISTERED RULE (DECISIONS 2026-07-10, before KS-to-wild was measured):
  * degree tertile cutoff = candidate-degree 67th percentile (FIXED).
  * hard negatives are non-co-complex pairs (never CORUM same-complex).
  * KS(hard-null, wild) is a REPORTED DIAGNOSTIC, never minimized by tuning.
  * contamination sensitivity (epsilon in {1,5,10}%) is a required output.
Firewall: NO DepMap anywhere here.

Run: python -m emperor.hardnegatives
"""
from __future__ import annotations
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .covariates import add_covariates
from .conformal import conformal_pvalues, benjamini_hochberg
from .nonconformity import nonconformity

DEGREE_TERTILE = 0.67   # PRE-REGISTERED cutoff (candidate degree percentile)


def _corum_same_complex_set():
    """frozenset(pair) -> in-some-CORUM-complex, to exclude from hard negatives."""
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")
    pos = lab[lab.label == 1]
    return {frozenset((a, b)) for a, b in zip(pos.uniprot_a, pos.uniprot_b)}


def build_hard_negatives():
    """Return decoy rows augmented with a `hard` flag per the pre-registered rule."""
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    inter = add_covariates(inter)
    cand = inter[inter.is_true_pair]
    cutoff = np.quantile(cand["degree"].to_numpy(float), DEGREE_TERTILE)

    decoy = inter[~inter.is_true_pair].copy()
    same = _corum_same_complex_set()
    is_cocomplex = [frozenset((a, b)) in same for a, b in zip(decoy.uniprot_a, decoy.uniprot_b)]
    decoy["is_cocomplex"] = is_cocomplex
    # hard = high-degree AND not a known co-complex pair
    decoy["hard"] = (decoy["degree"] >= cutoff) & (~decoy["is_cocomplex"])
    return decoy, float(cutoff)


def run():
    decoy, cutoff = build_hard_negatives()
    inter = add_covariates(pd.read_parquet(C.INTERIM / "interactome.parquet"))
    cand = inter[inter.is_true_pair].copy()

    # Nonconformity distributions
    s_wild = nonconformity(cand)
    s_decoy_all = nonconformity(decoy)
    s_hard = nonconformity(decoy[decoy.hard])
    s_soft = nonconformity(decoy[~decoy.hard])

    # DIAGNOSTIC ONLY (never minimized): KS of each null vs wild candidates
    ks_all = stats.ks_2samp(s_decoy_all, s_wild)
    ks_hard = stats.ks_2samp(s_hard, s_wild)

    # Certify against (a) plain decoy null and (b) hard-negative-augmented null.
    # Augmented null = union of all decoys + hard negatives (hard upweighted by
    # inclusion; here we simply calibrate against the hard subset as the null that
    # matches the candidate region, plus report the union).
    def certify_counts(s_calneg):
        p = conformal_pvalues(s_wild, s_calneg)
        out = {}
        hc = cand["high_conf"].to_numpy()
        for q in sorted(set(C.Q_SWEEP) | {C.Q}):
            r = benjamini_hochberg(p, q)
            out[q] = dict(certified=int(r.sum()), high_conf_dropped=int((hc & ~r).sum()))
        return out

    plain = certify_counts(s_decoy_all)
    hard_only = certify_counts(s_hard)
    # Feasible operationalization of the pre-registered rule (see DECISIONS 2026-07-10
    # deviation note): we cannot add NEW random pairs to the null (unfolded -> no
    # score), so we UPWEIGHT the hard (high-degree) region of the existing decoy null
    # by an explicit integer replication factor. This makes the null's mass sit where
    # the wild candidates sit (the hub-selected region) without inventing scores.
    HARD_UPWEIGHT = 2
    s_augmented = np.concatenate([s_decoy_all] + [s_hard] * HARD_UPWEIGHT)
    augmented = certify_counts(s_augmented)

    # CONTAMINATION SENSITIVITY (required output, markup #1): suppose a fraction eps
    # of the "hard negatives" are truly interacting (real-but-unannotated). Removing
    # them from the null (they should not be nulls) LOWERS the null's nonconformity
    # floor -> larger conformal p-values -> FEWER certified. We bound the certified
    # count by dropping the eps-fraction most interaction-like (lowest-s) hard negs.
    contam = {}
    s_hard_sorted = np.sort(s_hard)  # ascending: lowest s = most interaction-like
    for eps in (0.01, 0.05, 0.10):
        k = int(np.floor(eps * len(s_hard)))
        s_clean = s_hard_sorted[k:] if k > 0 else s_hard_sorted  # drop k most signal-like
        s_aug_clean = np.concatenate([s_decoy_all, s_clean])
        p = conformal_pvalues(s_wild, s_aug_clean)
        hc = cand["high_conf"].to_numpy()
        r = benjamini_hochberg(p, C.Q)
        contam[eps] = dict(certified=int(r.sum()), high_conf_dropped=int((hc & ~r).sum()),
                           n_dropped_from_null=k)

    out = dict(
        degree_cutoff=cutoff,
        n_decoy_total=int(len(decoy)),
        n_hard=int(decoy.hard.sum()),
        n_soft=int((~decoy.hard).sum()),
        ks_all_vs_wild=dict(ks=float(ks_all.statistic), p=float(ks_all.pvalue)),
        ks_hard_vs_wild=dict(ks=float(ks_hard.statistic), p=float(ks_hard.pvalue)),
        ks_note="DIAGNOSTIC ONLY — the hard-negative rule was pre-registered before this was measured and is NOT tuned to reduce it",
        certified=dict(plain_decoy_null=plain, hard_region_null=hard_only,
                       augmented_null=augmented),
        contamination_sensitivity=dict(
            note="eps = fraction of hard negatives that are truly interacting; "
                 "removing them makes p-values LARGER -> FEWER certified (anti-conservative "
                 "bias direction if they are wrongly kept as nulls). Bounds the augmented-null count.",
            by_eps=contam,
        ),
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "hard_negatives.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"degree cutoff (cand 67th pct) = {o['degree_cutoff']:.1f}")
    print(f"decoys: {o['n_decoy_total']} total -> {o['n_hard']} HARD, {o['n_soft']} soft")
    print(f"KS(all-decoy null, wild) = {o['ks_all_vs_wild']['ks']:.3f}  [diagnostic]")
    print(f"KS(hard null,       wild) = {o['ks_hard_vs_wild']['ks']:.3f}  [diagnostic, not tuned]")
    print("\ncertified@0.10 (cert / high-conf dropped):")
    for name, d in o["certified"].items():
        print(f"  {name:20s}: {d[0.1]['certified']:>4} / {d[0.1]['high_conf_dropped']:>3}")
    print("\ncontamination sensitivity (augmented null, q=0.10):")
    for eps, d in o["contamination_sensitivity"]["by_eps"].items():
        print(f"  eps={eps}: certified={d['certified']}  dropped={d['high_conf_dropped']}  (removed {d['n_dropped_from_null']} from null)")


if __name__ == "__main__":
    main()
