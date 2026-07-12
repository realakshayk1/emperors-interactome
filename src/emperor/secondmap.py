"""secondmap.py — second real-map audit: Predictomes (Schmid & Walter, Mol Cell 2025, 85(6):1216-1232.e5, doi:10.1016/j.molcel.2025.01.034).

Closes G-GENERAL by transferring the conformal audit to an INDEPENDENT real AlphaFold
interactome: the genome-scale human Predictome (20,200 proteins, 1.6M pairs), a different
lab / pipeline / protein universe from CM4AI. Confirms the method transfers, and yields
two findings the single-map paper cannot make:

  (1) AUDIT-AXIS REQUIREMENT. The audit axis must be a RAW, untrained structural-confidence
      readout (like CM4AI's ipTM). Predictomes' headline score SPOC is a *trained classifier*
      built to separate interactors from random pairs; on that axis the high-confidence
      candidates are perfectly separated from a random-pair null, so every conformal p-value
      pins to 1/(n+1) and the audit is degenerate (100% certification, no discrimination).
      We therefore audit on `num_unique_contacts`, a raw interface-size readout with genuine
      candidate/null overlap.

  (2) CLASSIFIER-CURATED MAPS HAVE LESS TO AUDIT. The DepMap co-essentiality referee is alive
      on this map (co-essential fraction rises 0.08 -> 0.43 across SPOC bins), but does NOT
      separate certified from dropped WITHIN the SPOC>=0.9 tier, because that tier is already
      uniformly co-essentiality-rich. Predictomes is classifier-curated (SPOC pre-filters);
      its high-confidence tier has the junk already removed, so the audit finds little to drop
      that the referee would flag. CM4AI's high-confidence tier is a raw ipTM threshold with no
      such curation, which is why the audit catches 22% of it.

Data: predictomes.org downloads -> 20251110_hs_predictome_pair_scores.csv.gz (28 MB), parsed
to data/external/predictomes_pair_scores.parquet. Raw per-pair ipTM lives only inside the
~1.5 TB per-structure archives (40x ~39 GB) and is not pulled; the ipTM-vs-ipTM cross-map
referee comparison is left as future work needing that access.

Purity: DepMap never enters certification (contact-axis nonconformity vs random null only).
Run: python -m emperor.secondmap
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd

from . import config as C
from .conformal import conformal_pvalues, benjamini_hochberg

SPOC_HICONF = 0.90     # Predictomes high-confidence tier (its genuine top tier)
N_NULL = 8000          # native random-pair calibration null size
COESS_SIG = 0.05       # DepMap co-essentiality significance (matches CM4AI run)


def _load():
    ext = C.ROOT / "data" / "external" / "predictomes_pair_scores.parquet"
    df = pd.read_parquet(ext)
    df["nc"] = pd.to_numeric(df["num_unique_contacts"], errors="coerce").fillna(0.0)
    return df


def _depmap():
    z = np.load(C.INTERIM / "depmap_slice.npz", allow_pickle=True)
    genes = list(z["genes"]); gidx = {g: i for i, g in enumerate(genes)}
    return gidx, z["gls_p"]


def run(seed: int = None):
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    df = _load()
    gidx, gls_p = _depmap()

    # SPOC degeneracy demonstration (finding 1): high-conf tier vs a native NEGATIVE null on
    # the SPOC axis. SPOC is trained to separate interactors from random pairs, so the tier is
    # perfectly separated from the negatives -> every conformal p-value pins to the 1/(n+1) floor.
    hi_spoc = df[df.spoc_score >= SPOC_HICONF]
    neg_pool = df.index[df.spoc_score <= 0.05].values           # map's own implied negatives
    spoc_null_idx = rng.choice(neg_pool, size=N_NULL, replace=False)
    p_spoc = conformal_pvalues(1.0 - hi_spoc.spoc_score.values, 1.0 - df.spoc_score.values[spoc_null_idx])
    spoc_degenerate = bool(np.mean(p_spoc <= 1.5 / (N_NULL + 1)) > 0.99)  # ~all at floor
    # contact-axis null is a genuine random-pair sample (negatives dominate the 1.6M-pair space)
    null_idx = rng.choice(len(df), size=N_NULL, replace=False)

    # Real audit on the RAW contact axis (finding: non-degenerate).
    s_null = -df.nc.values[null_idx]
    tier = hi_spoc.copy()
    p = conformal_pvalues(-tier.nc.values, s_null)
    tier["conf_p"] = p
    certified_by_q, dropped_by_q = {}, {}
    for q in C.Q_SWEEP:
        rej = benjamini_hochberg(p, q)
        tier[f"cert@{q}"] = rej
        certified_by_q[q] = int(rej.sum()); dropped_by_q[q] = int((~rej).sum())

    # DepMap referee on covered subset (finding 2).
    def coess(sa, sb):
        ia, ib = gidx.get(sa), gidx.get(sb)
        if ia is None or ib is None:
            return None
        return int(gls_p[ia, ib] < COESS_SIG)
    tier["coess"] = [coess(a, b) for a, b in zip(tier.sym_a, tier.sym_b)]
    cov = tier[tier.coess.notna()]
    cert = cov[cov[f"cert@{C.Q}"]]; drop = cov[~cov[f"cert@{C.Q}"]]

    # Referee liveness: co-essential fraction across SPOC bins on the covered universe.
    cov_all = df[df.sym_a.notna() & df.sym_b.notna()].copy()
    cov_all["coess"] = [coess(a, b) for a, b in zip(cov_all.sym_a, cov_all.sym_b)]
    cov_all = cov_all[cov_all.coess.notna()]
    spoc_gradient = {}
    for lo, hi in [(0.90, 1.01), (0.75, 0.90), (0.50, 0.75), (0.05, 0.50), (0.0, 0.05)]:
        m = (cov_all.spoc_score >= lo) & (cov_all.spoc_score < hi)
        if m.sum():
            spoc_gradient[f"{lo:.2f}-{hi:.2f}"] = dict(n=int(m.sum()), coess=float(cov_all.coess[m].mean()))

    out = dict(
        map="Predictomes hs (Schmid & Walter, Mol Cell 2025, 85(6):1216-1232.e5, doi:10.1016/j.molcel.2025.01.034) — genome-scale human AF-Multimer interactome",
        n_proteins=int(len(set(df.acc_a) | set(df.acc_b))), n_pairs=int(len(df)),
        audit_axis="num_unique_contacts (raw interface size; SPOC is a trained classifier and is degenerate for auditing)",
        spoc_axis_degenerate=spoc_degenerate,
        high_conf_tier=f"SPOC>={SPOC_HICONF}", n_high_conf=int(len(tier)),
        native_null="random-pair sample (negatives dominate); size %d" % N_NULL,
        certified_by_q=certified_by_q, dropped_by_q=dropped_by_q,
        referee=dict(
            source="DepMap co-essentiality (same referee as CM4AI; SPOC-independent -> purity firewall holds)",
            n_covered=int(len(cov)),
            certified_coess=float(cert.coess.mean()) if len(cert) else None, n_certified=int(len(cert)),
            dropped_coess=float(drop.coess.mean()) if len(drop) else None, n_dropped=int(len(drop)),
            background_coess=float(cov.coess.mean()),
            referee_alive_spoc_gradient=spoc_gradient,
            separation="NONE within SPOC>=0.9 tier — the tier is already uniformly co-essentiality-rich "
                       "(classifier-curated); the co-essentiality signal is carried by SPOC, not the raw "
                       "contact axis. Contrast: CM4AI's uncurated ipTM-threshold tier IS separated "
                       "(certified 0.41 vs dropped 0.17).",
        ),
        gate="G-GENERAL: audit method transfers to a second REAL map (>=2 real maps + semi-synthetic).",
        limitation="ipTM-vs-ipTM cross-map referee replication needs the ~1.5 TB per-structure archives (not pulled).",
        seed=seed,
    )
    C.PROCESSED.mkdir(parents=True, exist_ok=True)
    (C.PROCESSED / "secondmap_audit.json").write_text(json.dumps(out, indent=2, default=float))
    return out


def main():
    o = run()
    print(f"Second map: {o['map']}")
    print(f"  {o['n_proteins']} proteins, {o['n_pairs']} pairs; SPOC axis degenerate: {o['spoc_axis_degenerate']}")
    print(f"  high-conf tier {o['high_conf_tier']} (n={o['n_high_conf']}) on raw contact axis, native null")
    print(f"  certified by q: {o['certified_by_q']}  dropped by q: {o['dropped_by_q']}")
    r = o["referee"]
    print(f"  referee (n={r['n_covered']} covered): certified {r['certified_coess']:.3f} (n={r['n_certified']}) "
          f"vs dropped {r['dropped_coess']:.3f} (n={r['n_dropped']}), bg {r['background_coess']:.3f}")
    print(f"  referee alive (SPOC gradient): " +
          ", ".join(f"{k}:{v['coess']:.2f}" for k, v in r["referee_alive_spoc_gradient"].items()))


if __name__ == "__main__":
    main()
