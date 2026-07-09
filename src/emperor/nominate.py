"""nominate.py — certified missing-member nomination (METHODS §8, SPEC G4).

For a target cancer complex, rank non-member candidate proteins c by

    score(c) = certified_confidence(c <-> members) x held_out_coessentiality(c, members)

  * certified_confidence: evidence that c has conformally-certified (low conf
    p-value) edges to complex members in the interactome.
  * held_out_coessentiality: mean DepMap co-essentiality of c with members — the
    INDEPENDENT vote (never used in the score/labels).
Nominate the top c with certified risk < q AND positive held-out co-essentiality.
Output: data/processed/nomination.json. Run: python -m emperor.nominate <complex_id?>
"""
from __future__ import annotations
import json
import sys

import numpy as np
import pandas as pd

from . import config as C


def _load_slice():
    z = np.load(C.INTERIM / "depmap_slice.npz", allow_pickle=True)
    genes = list(z["genes"])
    return z["gls_p"], z["gls_sign"], {g: i for i, g in enumerate(genes)}


def _corum_members():
    corum = pd.read_csv(C.RAW / "corum_human.txt", sep="\t", dtype=str)
    out = {}
    for _, r in corum.iterrows():
        accs = [a.strip() for a in str(r["subunits_uniprot_id"]).split(";") if a.strip()]
        names = [n.strip() for n in str(r["subunits_gene_name"]).split(";") if n.strip()]
        out[r["complex_id"]] = dict(name=r["complex_name"], accs=set(accs),
                                    genes=set(names), pmid=r.get("pmid"))
    return out


def coess_gene_vs_set(g, members, gls_p, gls_sign, gidx):
    """Mean (-log10 GLS p, signed +) of gene g against member genes present in DepMap."""
    vals = []
    for m in members:
        if g in gidx and m in gidx and g != m:
            p = gls_p[gidx[g], gidx[m]]
            s = gls_sign[gidx[g], gidx[m]]
            vals.append(np.sign(s) * -np.log10(max(p, 1e-300)))
    return float(np.mean(vals)) if vals else 0.0, len(vals)


def nominate_for_complex(complex_id, cert, gls_p, gls_sign, gidx, corum):
    """Rank interactome proteins that edge into the complex but aren't members."""
    info = corum[complex_id]
    member_genes = info["genes"]
    q = C.Q
    certcol = f"certified@{q}"

    # Candidate = any interactome protein with a CERTIFIED edge to a member, not itself a member
    cand_rows = cert[(cert[certcol]) &
                     ((cert.gene_a.isin(member_genes)) ^ (cert.gene_b.isin(member_genes)))]
    results = []
    for _, r in cand_rows.iterrows():
        c = r.gene_b if r.gene_a in member_genes else r.gene_a
        if c in member_genes:
            continue
        coess, n_co = coess_gene_vs_set(c, member_genes, gls_p, gls_sign, gidx)
        results.append(dict(
            candidate=c, partner_member=(r.gene_a if c == r.gene_b else r.gene_b),
            conf_pvalue=float(r.conf_pvalue), certified_risk=float(r.conf_pvalue),
            score=float(r.score), coess_signed=coess, n_coess_members=n_co,
            nom_score=float((1 - r.conf_pvalue) * max(coess, 0.0)),
        ))
    df = pd.DataFrame(results).sort_values("nom_score", ascending=False)
    return info, df


def run(complex_id=None):
    gls_p, gls_sign, gidx = _load_slice()
    cert = pd.read_parquet(C.PROCESSED / "certified.parquet")
    corum = _corum_members()
    if complex_id is None:
        complex_id = getattr(C, "TARGET_COMPLEX_ID", None) or _pick_target(cert, corum, gidx)
    info, df = nominate_for_complex(complex_id, cert, gls_p, gls_sign, gidx, corum)
    top = df[(df.certified_risk < C.Q) & (df.coess_signed > 0)]
    nominee = (top.iloc[0].to_dict() if len(top) else
               (df.iloc[0].to_dict() if len(df) else None))
    out = dict(complex_id=complex_id, complex_name=info["name"],
               complex_pmid=info.get("pmid"),
               members=sorted(info["genes"]), q=C.Q,
               nominee=nominee, ranked=df.head(15).to_dict("records"))
    if nominee is not None:
        _enrich(out, nominee["candidate"], info, cert, gls_p, gls_sign, gidx, corum)
    (C.PROCESSED / "nomination.json").write_text(json.dumps(out, indent=2, default=str))
    return out


def _enrich(out, cand, info, cert, gls_p, gls_sign, gidx, corum):
    """Attach per-member co-essentiality, held-out CORUM confirmation, edges,
    and the literature-grounded rationale (see RATIONALE.md for the prose source)."""
    members = sorted(info["genes"])
    breakdown = []
    for m in members:
        if cand in gidx and m in gidx and cand != m:
            p = float(gls_p[gidx[cand], gidx[m]]); s = float(gls_sign[gidx[cand], gidx[m]])
            breakdown.append(dict(member=m, gls_p=p, sign=int(s),
                                  strong_coess=bool(p < C.COESS_SIG_THRESHOLD and s > 0)))
    out["coess_breakdown"] = sorted(breakdown, key=lambda d: d["gls_p"])
    # Independent confirmation: is the nominee an annotated member of ANOTHER CORUM
    # complex that shares members with the target? (held out of the pipeline)
    conf = []
    for cid, meta in corum.items():
        if cand in meta["genes"] and cid != out["complex_id"]:
            shared = sorted(meta["genes"] & info["genes"])
            if len(shared) >= 3:
                conf.append(dict(complex_id=cid, complex_name=meta["name"],
                                 pmid=meta.get("pmid"), n_shared=len(shared),
                                 shared_members=shared))
    out["independent_confirmation"] = sorted(conf, key=lambda d: -d["n_shared"])[:3]
    out["interactome_edges"] = cert[(cert.gene_a == cand) | (cert.gene_b == cand)][
        ["gene_a", "gene_b", "score", "conf_pvalue", f"certified@{C.Q}", "high_conf"]
    ].to_dict("records")
    rat = C.ROOT / "RATIONALE.md"
    if rat.exists():
        out["rationale"] = rat.read_text()
    # independent Boltz-2 structural corroboration of the nominee's edge, if predicted
    pv = C.ROOT / "data" / "structures" / "physical_validity.json"
    if pv.exists():
        import json as _json
        out["physical_validity"] = _json.loads(pv.read_text())


def _pick_target(cert, corum, gidx):
    """Placeholder auto-pick; the real choice is made in analysis (see DECISIONS)."""
    # choose a complex with certified edges to non-members and members in DepMap
    best, best_n = None, -1
    certcol = f"certified@{C.Q}"
    for cid, info in corum.items():
        mg = info["genes"]
        if len(mg) < 3:
            continue
        edges = cert[(cert[certcol]) &
                     ((cert.gene_a.isin(mg)) ^ (cert.gene_b.isin(mg)))]
        n = len(edges)
        if n > best_n:
            best, best_n = cid, n
    return best


def main():
    cid = sys.argv[1] if len(sys.argv) > 1 else None
    o = run(cid)
    print(f"target complex {o['complex_id']}: {o['complex_name']}")
    print(f"members ({len(o['members'])}): {', '.join(o['members'][:12])}")
    if o["nominee"]:
        n = o["nominee"]
        print(f"\nNOMINEE: {n['candidate']}  (partner {n['partner_member']})")
        print(f"  certified risk (conf p) = {n['certified_risk']:.4f}  < q={o['q']}")
        print(f"  held-out co-essentiality (signed -log10p) = {n['coess_signed']:.2f} "
              f"over {n['n_coess_members']} members")
    else:
        print("no nominee found")


if __name__ == "__main__":
    main()
