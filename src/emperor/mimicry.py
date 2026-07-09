"""mimicry.py — second-map audit: Krogan-lab host–pathogen molecular-mimicry interactome.
Paper: "AlphaFold models of host–pathogen interactions elucidate the prevalence and structural
modes of molecular mimicry" (Krogan et al., bioRxiv 2025.06.04.657796). Data: Zenodo record
15588019. (First-author name not independently verified in-session; cited per the project's
source brief as Krogan et al.)

DAY-1 GATE finding (verified firsthand 2026-07): the Zenodo deposit ships AF-Multimer
ranked-0 PDB models + FASTA, but NO tabulated per-pair confidence scores. ipTM is NOT
recoverable (no deposited PAE matrix — that needs GPU recompute). pDockQ IS computable on
CPU directly from each ranked-0 two-chain PDB (pLDDT in the B-factor column). So the
confidence axis for THIS map is **pDockQ**, distinct from CM4AI's ipTM — which is exactly
why calibration must be strictly per-map (never pool raw scores across maps).

pDockQ (Bryant, Pozzati & Elofsson 2022, Nat Commun 13:1265):
    pDockQ = 0.724 / (1 + exp(-0.052*(x - 152.611))) + 0.018
    x = <mean interface pLDDT> * ln(number of interface residue–residue contacts)
    interface = residues with any heavy atom within 8 Å of the other chain.

PURITY FIREWALL: the validator for this map is the paper's experimental binding assays
(primary, clean positives) + mimicry-interface reuse as a confirmatory LOWER BOUND only
(explicitly non-independent — the mimicry call is itself AF-derived, circular with the score).
DepMap does NOT apply (cross-species host–pathogen). No validator enters the pDockQ score.
"""
from __future__ import annotations
import math, glob, os, json
import numpy as np
from scipy.spatial import cKDTree


def _parse_cb(path):
    """One representative atom per residue: CB, or CA for glycine / when CB absent.
    Returns per-chain list of (resid, x, y, z, bfactor=pLDDT)."""
    # first pass: collect CB and CA per (chain,resid)
    cb, ca = {}, {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom not in ("CB", "CA"):
            continue
        ch = line[21]
        resid = int(line[22:26])
        rec = (ch, resid)
        x, y, z, b = float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66])
        (cb if atom == "CB" else ca)[rec] = (x, y, z, b)
    chains = {}
    for rec in set(cb) | set(ca):
        x, y, z, b = cb.get(rec, ca.get(rec))
        ch, resid = rec
        chains.setdefault(ch, []).append((resid, x, y, z, b))
    for ch in chains:
        chains[ch].sort()
    return chains


def pdockq(path, dist=8.0):
    """CPU pDockQ from a two-chain AF-M ranked-0 PDB (Bryant 2022: Cβ–Cβ contacts within 8 Å,
    Cα for glycine). Returns dict or None (not a dimer)."""
    ch = _parse_cb(path)
    cks = sorted(ch)
    if len(cks) != 2:
        return None
    A, B = ch[cks[0]], ch[cks[1]]
    Ac = np.array([(a[1], a[2], a[3]) for a in A])
    Bc = np.array([(b[1], b[2], b[3]) for b in B])
    ta, tb = cKDTree(Ac), cKDTree(Bc)
    pairs = ta.query_ball_tree(tb, dist)          # residue(Cβ)-level neighbor lists
    if_A, if_B, n_contacts = set(), set(), 0
    for i, js in enumerate(pairs):
        if js:
            if_A.add(A[i][0])
            for j in js:
                if_B.add(B[j][0])
            n_contacts += len(js)                 # # of inter-chain Cβ–Cβ residue contacts
    if not if_A or not if_B:
        return {"pdockq": 0.0, "n_if_res": 0, "n_contacts": 0, "mean_if_plddt": 0.0,
                "chain_a": cks[0], "chain_b": cks[1]}
    pl = [a[4] for a in A if a[0] in if_A] + [b[4] for b in B if b[0] in if_B]
    mean_if_plddt = float(np.mean(pl))
    x = mean_if_plddt * math.log(max(n_contacts, 1))
    pd = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
    return {"pdockq": round(pd, 4), "n_if_res": len(if_A) + len(if_B),
            "n_contacts": int(n_contacts), "mean_if_plddt": round(mean_if_plddt, 2),
            "chain_a": cks[0], "chain_b": cks[1]}


def score_dir(pdb_glob):
    """Score every PDB matching pdb_glob. Returns list of dicts with a 'model' key."""
    out = []
    for f in sorted(glob.glob(pdb_glob)):
        r = pdockq(f)
        if r:
            r["model"] = os.path.basename(f).replace("-1.pdb", "").replace(".pdb", "")
            out.append(r)
    return out


if __name__ == "__main__":
    import sys
    g = sys.argv[1] if len(sys.argv) > 1 else "data/mimicry/models/benchmark_afmultimer_ranked0/*.pdb"
    res = score_dir(g)
    pds = np.array([r["pdockq"] for r in res])
    print(f"scored {len(res)} models")
    print("pDockQ pctiles [5,25,50,75,95]:", [round(float(np.percentile(pds, p)), 3) for p in (5, 25, 50, 75, 95)])
    print("frac >=0.23 (acceptable):", round(float((pds >= 0.23).mean()), 3),
          "| frac >=0.5 (confident):", round(float((pds >= 0.5).mean()), 3))
