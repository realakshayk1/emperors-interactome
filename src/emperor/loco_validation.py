"""loco_validation.py — leave-one-complex-member-out recovery (G-VALID, WS3.1 substitute).

Construct: does the conformal-certified edge set recover held-out true complex members
better than chance? For each CORUM complex with >=2 members in the candidate protein
universe, a member is "recovered" if it has at least one conformal-certified edge (q)
to ANOTHER member of the same complex. Recovery is genuinely held-out: certification
depends only on ipTM vs the calibration decoys (nonconformity), never on complex
membership, so hiding a member's label changes nothing about whether its edges certify.

Null: eligible impostors. For each complex, non-members that nonetheless have >=1
candidate edge INTO the complex are the fair comparison group (they are eligible to be
recovered, so the test is not the trivial "a random protein has no edge"). Significance
by within-complex member/impostor label permutation.

Purity: the referee never enters certification (ipTM vs decoys only). Reports member vs
impostor recovery, OR, permutation p.
Run: python -m emperor.loco_validation
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .conformal import certify


def _union_degree(inter: pd.DataFrame) -> dict:
    deg: dict = {}
    for a, b in zip(inter.gene_a, inter.gene_b):
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    return deg


def run(q: float = None, n_perm: int = 10000, seed: int = None, out_path=None):
    q = C.Q if q is None else q
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)

    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    lab = pd.read_parquet(C.INTERIM / "labels.parquet")

    # Candidate edge set = the true-pair candidates the audit operates on.
    cand = inter[inter.is_true_pair].copy()
    cal_neg = lab[(lab.label == 0) & (lab.split == "cal")]["score"].to_numpy()
    cal_neg = (1.0 - cal_neg)  # nonconformity of calibration negatives (w_phys=0 shipped run)

    # Certify candidate edges (membership-blind).
    cert = certify(cand, cal_neg, q=q, w_phys=0.0)
    cmask = cert[f"certified@{q}"].to_numpy()
    cert_edges = set()
    for keep, a, b in zip(cmask, cert.gene_a, cert.gene_b):
        if keep:
            cert_edges.add(frozenset((a, b)))

    # Complex membership (positives carry complex_id).
    pos = lab[(lab.label == 1) & lab.complex_id.notna()]
    cand_prots = set(cand.gene_a) | set(cand.gene_b)
    # map complex -> set of member genes that are in the candidate universe
    comp_members: dict = {}
    for cid, g_a, g_b in zip(pos.complex_id, pos.gene_a, pos.gene_b):
        s = comp_members.setdefault(cid, set())
        for g in (g_a, g_b):
            if g in cand_prots:
                s.add(g)
    comp_members = {c: m for c, m in comp_members.items() if len(m) >= 2}

    # Candidate-graph adjacency (union of true-pair candidate edges).
    from collections import defaultdict
    adj = defaultdict(set)
    for a, b in zip(cand.gene_a, cand.gene_b):
        adj[a].add(b); adj[b].add(a)

    def recovered(protein, members):
        """Does `protein` have a CERTIFIED edge to any OTHER member of `members`?"""
        for m in members:
            if m != protein and frozenset((protein, m)) in cert_edges:
                return True
        return False

    # Fair test with an ELIGIBLE-IMPOSTOR null. For each complex we form two groups:
    #   members   — true members (in candidate universe);
    #   impostors — NON-members that nonetheless have >=1 candidate edge INTO the complex
    #               (so they are eligible to be "recovered" — this removes the trivial
    #               confound that a random protein simply has no edge into the complex).
    # A trial is "recovered" if it has a certified within-complex edge. We test whether
    # true membership predicts recovery ABOVE eligible impostors, via a label-permutation
    # test that shuffles member/impostor labels WITHIN each complex (controls for each
    # complex's own certification density).
    per_complex = []  # (member_flags: np.bool_ array, recovered_flags: np.bool_ array)
    mem_rec = imp_rec = mem_n = imp_n = 0
    for cid, members in comp_members.items():
        impostors = set()
        for m in members:
            for nb in adj[m]:
                if nb not in members:
                    impostors.add(nb)
        trial_prots = list(members) + list(impostors)
        is_member = np.array([True] * len(members) + [False] * len(impostors))
        rec = np.array([recovered(p, members) for p in trial_prots], dtype=bool)
        if is_member.any():
            per_complex.append((is_member, rec))
            mem_rec += int(rec[is_member].sum()); mem_n += int(is_member.sum())
            imp_rec += int(rec[~is_member].sum()); imp_n += int((~is_member).sum())

    obs_member_rate = mem_rec / mem_n
    obs_impostor_rate = imp_rec / imp_n if imp_n else 0.0
    obs_diff = obs_member_rate - obs_impostor_rate

    # Within-complex label permutation: reshuffle which trials are "members" (keeping the
    # per-complex member count fixed), recompute the member-minus-impostor recovery gap.
    null_diff = np.empty(n_perm)
    for bx in range(n_perm):
        mr = mn = ir = ln = 0
        for is_member, rec in per_complex:
            k = int(is_member.sum()); tot = is_member.size
            idx = rng.permutation(tot)[:k]
            perm_mem = np.zeros(tot, dtype=bool); perm_mem[idx] = True
            mr += int(rec[perm_mem].sum()); mn += k
            ir += int(rec[~perm_mem].sum()); ln += tot - k
        null_diff[bx] = (mr / mn) - (ir / ln if ln else 0.0)
    perm_p = float((null_diff >= obs_diff).mean())

    def OR(p, qn):
        p = min(max(p, 1e-9), 1 - 1e-9); qn = min(max(qn, 1e-9), 1 - 1e-9)
        return (p / (1 - p)) / (qn / (1 - qn))

    orr = float(OR(obs_member_rate, obs_impostor_rate))
    out = dict(
        instrument="leave-one-complex-member-out recovery vs eligible-impostor null (G-VALID temporal substitute)",
        q=q, n_complexes=len(per_complex),
        n_member_trials=mem_n, n_impostor_trials=imp_n,
        n_certified_edges=int(cmask.sum()),
        observed_member_recovery=obs_member_rate,
        observed_impostor_recovery=obs_impostor_rate,
        recovery_gap=obs_diff,
        odds_ratio_member_vs_impostor=orr,
        permutation_p=perm_p, n_permutations=n_perm, seed=seed,
        null="within-complex member/impostor label permutation (impostor = non-member with >=1 candidate edge into the complex)",
        preregistered_thresholds=dict(min_or=2.0, max_p=0.05),
        gate_met=bool(orr >= 2.0 and perm_p < 0.05),
    )
    dest = out_path if out_path is not None else (C.PROCESSED / "loco_validation.json")
    from pathlib import Path
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"LOCO recovery vs eligible-impostor null ({o['n_certified_edges']} certified edges, "
          f"{o['n_complexes']} complexes)")
    print(f"  members   {o['n_member_trials']} trials -> recovery {o['observed_member_recovery']:.3f}")
    print(f"  impostors {o['n_impostor_trials']} trials -> recovery {o['observed_impostor_recovery']:.3f}")
    print(f"  gap={o['recovery_gap']:.3f}  OR={o['odds_ratio_member_vs_impostor']:.1f}  "
          f"perm_p={o['permutation_p']:.4f}  GATE {'MET' if o['gate_met'] else 'NOT MET'}")


if __name__ == "__main__":
    main()
