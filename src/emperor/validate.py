"""validate.py — held-out DepMap co-essentiality validation (METHODS §6, SPEC G3).

The certified set is convincing only if an INDEPENDENT signal agrees. DepMap
co-essentiality (GLS p-value/sign) is used ONLY here — never in the score or
labels (purity firewall). For three edge sets — certified@q, the paper's raw
high-confidence set, and dropped-artifact edges — we compute co-essentiality
enrichment and a permutation p-value.

Enrichment metrics per edge set:
  * frac_coess: fraction of edges that are a "strong co-essential pair"
    (GLS p <= COESS_SIG_THRESHOLD and positive sign).
  * mean_neglog_p: mean of -log10(GLS p) over edges with both genes in DepMap.
Permutation p: shuffle set membership across the union of scored edges.
Output: data/processed/validation.json. Run: python -m emperor.validate
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C


def _load_slice():
    z = np.load(C.INTERIM / "depmap_slice.npz", allow_pickle=True)
    genes = list(z["genes"])
    gidx = {g: i for i, g in enumerate(genes)}
    return z["gls_p"], z["gls_sign"], gidx


def _edge_coess(df, gls_p, gls_sign, gidx):
    """Return per-edge (neglog_p, is_strong) for edges with both genes in DepMap."""
    neglog, strong, mask = [], [], []
    for a, b in zip(df.gene_a, df.gene_b):
        if a in gidx and b in gidx and a != b:
            i, j = gidx[a], gidx[b]
            p = gls_p[i, j]
            s = gls_sign[i, j]
            nl = -np.log10(max(p, 1e-300))
            neglog.append(nl)
            strong.append(bool(p <= C.COESS_SIG_THRESHOLD and s > 0))
            mask.append(True)
        else:
            mask.append(False)
    return np.array(neglog), np.array(strong, dtype=bool), np.array(mask, dtype=bool)


def _stats(neglog, strong):
    return dict(n_scored=int(len(neglog)),
                frac_coess=float(strong.mean()) if len(strong) else 0.0,
                mean_neglog_p=float(neglog.mean()) if len(neglog) else 0.0)


def run():
    gls_p, gls_sign, gidx = _load_slice()
    cert = pd.read_parquet(C.PROCESSED / "certified.parquet")
    q = C.Q
    certcol = f"certified@{q}"

    sets = {
        "certified": cert[cert[certcol]],
        "raw_high_conf": cert[cert["high_conf"]],
        "dropped": cert[cert["high_conf"] & ~cert[certcol]],  # paper-HC but not certified
    }
    out = {"q": q, "sets": {}}
    edge_data = {}
    for name, df in sets.items():
        nl, strong, mask = _edge_coess(df, gls_p, gls_sign, gidx)
        out["sets"][name] = dict(n_edges=int(len(df)), **_stats(nl, strong))
        edge_data[name] = (nl, strong)

    # Permutation test: certified vs raw_high_conf on frac_coess and mean_neglog_p.
    # Pool the two sets' scored edges, reshuffle labels, recompute the difference.
    rng = np.random.default_rng(C.SEED)
    nl_c, st_c = edge_data["certified"]
    nl_r, st_r = edge_data["raw_high_conf"]
    obs_frac = st_c.mean() - st_r.mean()
    obs_nl = nl_c.mean() - nl_r.mean()
    pool_st = np.concatenate([st_c, st_r]).astype(float)
    pool_nl = np.concatenate([nl_c, nl_r])
    n_c = len(st_c)
    perm_frac = np.empty(C.N_PERMUTATIONS)
    perm_nl = np.empty(C.N_PERMUTATIONS)
    idx = np.arange(len(pool_st))
    for k in range(C.N_PERMUTATIONS):
        rng.shuffle(idx)
        a, b = idx[:n_c], idx[n_c:]
        perm_frac[k] = pool_st[a].mean() - pool_st[b].mean()
        perm_nl[k] = pool_nl[a].mean() - pool_nl[b].mean()
    out["permutation"] = dict(
        obs_frac_diff=float(obs_frac),
        p_frac=float((np.abs(perm_frac) >= abs(obs_frac)).mean()),
        obs_neglogp_diff=float(obs_nl),
        p_neglogp=float((np.abs(perm_nl) >= abs(obs_nl)).mean()),
        n_permutations=C.N_PERMUTATIONS,
    )
    # Also certified vs dropped (the sharper contrast)
    nl_d, st_d = edge_data["dropped"]
    if len(st_d):
        obs = st_c.mean() - st_d.mean()
        pool = np.concatenate([st_c, st_d]).astype(float)
        perm = np.empty(C.N_PERMUTATIONS)
        idx = np.arange(len(pool))
        for k in range(C.N_PERMUTATIONS):
            rng.shuffle(idx)
            perm[k] = pool[idx[:n_c]].mean() - pool[idx[n_c:]].mean()
        out["permutation_vs_dropped"] = dict(
            obs_frac_diff=float(obs), p_frac=float((np.abs(perm) >= abs(obs)).mean()))

    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "validation.json").write_text(json.dumps(out, indent=2))
    return out


def main():
    o = run()
    print(f"held-out DepMap co-essentiality enrichment (q={o['q']}):")
    for name, s in o["sets"].items():
        print(f"  {name:14s}: n={s['n_edges']:4d} scored={s['n_scored']:4d} "
              f"frac_coess={s['frac_coess']:.3f} mean_-log10p={s['mean_neglog_p']:.2f}")
    p = o["permutation"]
    print(f"\ncertified vs raw high-conf: Δfrac_coess={p['obs_frac_diff']:+.3f} "
          f"(perm p={p['p_frac']:.4f}); Δmean_-log10p={p['obs_neglogp_diff']:+.2f} "
          f"(perm p={p['p_neglogp']:.4f})")
    if "permutation_vs_dropped" in o:
        d = o["permutation_vs_dropped"]
        print(f"certified vs dropped: Δfrac_coess={d['obs_frac_diff']:+.3f} (perm p={d['p_frac']:.4f})")


if __name__ == "__main__":
    main()
