"""nominate_sets.py — nomination as FDR-controlled SELECTION.

The single flagship nomination (KANSL3 -> MLL1-WDR5) generalises: every conformally
CERTIFIED edge at level q that bridges a CORUM complex MEMBER to a NON-member is a
missing-member nomination, and — because the certified set is Benjamini-Hochberg
FDR-controlled at q — the resulting collection of nominations inherits the SAME
distribution-free guarantee: among all certified nominations, at most a fraction q are
expected to be wrong (upper bound; CORUM incompleteness). So instead of one hand-picked
nominee we emit, per complex, a certified nomination SET with a stated risk level.

DepMap co-essentiality (held-out) is used ONLY to RANK within a set, never to select
(purity firewall). Selection is purely conformal.
"""
from __future__ import annotations
import json, collections
import numpy as np
import pandas as pd
from . import config as C


def _corum():
    members, name = collections.defaultdict(set), {}
    path = C.RAW / "corum_human.txt"
    for line in open(path, encoding="utf-8", errors="ignore"):
        p = line.rstrip("\n").split("\t")
        if len(p) > 12 and p[0] != "complex_id":
            cid = p[0]; name[cid] = p[1]
            for u in p[11].split(";"):
                if u:
                    members[cid].add(u)
    return members, name


def build(q=0.1):
    df = pd.read_parquet(C.PROCESSED / "certified.parquet")
    cert = df[df[f"certified@{q}"]].copy()
    members, name = _corum()

    # optional held-out ranking signal (DepMap); load if present, else rank by conf p-value
    coess = {}
    slice_path = C.INTERIM / "depmap_slice.npz"
    if slice_path.exists():
        z = np.load(slice_path, allow_pickle=True)
        genes = list(z["genes"]); gls_p = z["gls_p"]; gls_s = z["gls_sign"]
        gi = {g: i for i, g in enumerate(genes)}
        def gene_vs_set(gene, memset_genes):
            if gene not in gi:
                return None
            vals = []
            for m in memset_genes:
                if m in gi and m != gene:
                    p = gls_p[gi[gene], gi[m]]; s = gls_s[gi[gene], gi[m]]
                    vals.append(s * -np.log10(max(p, 1e-300)))
            return float(np.mean(vals)) if vals else None
    else:
        gene_vs_set = None

    # gene symbols per uniprot (from the interactome table)
    im = pd.read_parquet(C.INTERIM / "interactome.parquet")
    sym = {}
    for _, r in im.iterrows():
        sym[r["uniprot_a"]] = r["gene_a"]; sym[r["uniprot_b"]] = r["gene_b"]

    # gene members per complex (for DepMap lookup)
    gmembers = {cid: {sym.get(u, u) for u in mem} for cid, mem in members.items()}

    nsets = {}
    for _, r in cert.iterrows():
        a, b = r["uniprot_a"], r["uniprot_b"]
        for cid, mem in members.items():
            inside_a, inside_b = a in mem, b in mem
            if inside_a ^ inside_b:
                nominee_u = b if inside_a else a
                nominee_g = sym.get(nominee_u, nominee_u)
                partner_u = a if inside_a else b
                if nominee_g in gmembers[cid]:   # already a gene-level member; skip
                    continue
                cs = gene_vs_set(nominee_g, gmembers[cid]) if gene_vs_set else None
                nsets.setdefault(cid, {"complex_id": cid, "complex_name": name.get(cid, ""),
                                       "n_members": len(mem), "nominations": []})
                nsets[cid]["nominations"].append({
                    "nominee": nominee_g, "nominee_uniprot": nominee_u,
                    "via_member": sym.get(partner_u, partner_u),
                    "conf_pvalue": round(float(r["conf_pvalue"]), 6),
                    "afm_score": round(float(r["score"]), 4),
                    "coess_signed_neglogp": None if cs is None else round(cs, 3)})

    # de-dup nominees within a complex (keep best conf_pvalue), sort by coess then conf
    out = []
    for cid, s in nsets.items():
        best = {}
        for n in s["nominations"]:
            k = n["nominee"]
            if k not in best or n["conf_pvalue"] < best[k]["conf_pvalue"]:
                best[k] = n
        noms = sorted(best.values(),
                      key=lambda n: (-(n["coess_signed_neglogp"] or -1e9), n["conf_pvalue"]))
        s["nominations"] = noms
        s["n_nominations"] = len(noms)
        out.append(s)
    out.sort(key=lambda s: -s["n_nominations"])

    result = {"q": q, "fdr_guarantee": f"<= {q} expected-false among certified nominations (upper bound; CORUM incomplete)",
              "n_complexes_with_nominations": len(out),
              "n_total_nominations": sum(s["n_nominations"] for s in out),
              "selection": "conformal-BH certified edges only (DepMap used for ranking, not selection)",
              "complexes": out}
    (C.PROCESSED / "nomination_sets.json").write_text(json.dumps(result, indent=1))
    return result


if __name__ == "__main__":
    r = build(q=0.1)
    print(f"complexes with >=1 certified nomination: {r['n_complexes_with_nominations']}")
    print(f"total certified nominations: {r['n_total_nominations']}  (FDR <= {r['q']}, upper bound)")
    print("\ntop complexes by nomination count:")
    for s in r["complexes"][:10]:
        noms = ", ".join(f"{n['nominee']}(via {n['via_member']})" for n in s["nominations"][:3])
        print(f"  [{s['complex_id']}] {s['complex_name'][:42]:42s} n={s['n_nominations']}  {noms}")
