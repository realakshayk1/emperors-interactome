"""labels.py — calibration labels from CORUM 5.3 + native map decoys.

Design (verified against the data, see DECISIONS):
  Positives (label 1): interactome pairs whose two proteins co-occur in the SAME
    CORUM complex — AF-M-scored AND independently CORUM-validated. Selected on
    CORUM membership ALONE (independent of the map's True/random flag), so this
    includes the 1 decoy-flagged pair that is in fact a CORUM complex: 191 total
    (190 from True candidates + 1 from the random set), all removed from negatives.
  Negatives (label 0): the map's NATIVE random decoy pairs (1,788; same AF-M
    pipeline, DepMap-independent). Only 1 is accidentally CORUM-same-complex and
    is dropped. A stricter negative set (both proteins in complexes but never the
    same one -> `diff_complex`) is tagged for the pos:neg sensitivity analysis.
  Split: complex-DISJOINT cal/test for positives (a CORUM complex's pairs never
    span both splits -> no leakage). Negatives split by seeded hash.

Output: data/interim/labels.parquet with SPEC columns
  [uniprot_a, uniprot_b, label, complex_id, split, score, neg_kind]
Purity: NO DepMap here (held-out firewall). Run: python -m emperor.labels
"""
from __future__ import annotations
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd

from . import config as C


def _corum_pairs() -> tuple[dict[frozenset, str], dict[str, set]]:
    corum = pd.read_csv(C.RAW / "corum_human.txt", sep="\t", dtype=str)
    same: dict[frozenset, str] = {}
    prot2c: dict[str, set] = defaultdict(set)
    for accs, cid in zip(corum["subunits_uniprot_id"], corum["complex_id"]):
        if not isinstance(accs, str):
            continue
        al = sorted({a.strip() for a in accs.split(";") if a.strip()})
        for a in al:
            prot2c[a].add(cid)
        for a, b in itertools.combinations(al, 2):
            same.setdefault(frozenset((a, b)), cid)
    return same, prot2c


def main() -> None:
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    same, prot2c = _corum_pairs()

    def status(a, b):
        fs = frozenset((a, b))
        if fs in same:
            return "same_complex", same[fs]
        if a in prot2c and b in prot2c:
            return "diff_complex", None
        return "unannotated", None

    st = inter.apply(lambda r: status(r.uniprot_a, r.uniprot_b), axis=1)
    inter["corum_status"] = [s[0] for s in st]
    inter["complex_id"] = [s[1] for s in st]

    # Positives: same-complex pairs (independent of the True/random flag)
    pos = inter[inter["corum_status"] == "same_complex"].copy()
    pos["label"] = 1
    pos["neg_kind"] = None

    # Negatives: native random decoys that are NOT accidentally same-complex
    neg = inter[(~inter["is_true_pair"]) & (inter["corum_status"] != "same_complex")].copy()
    neg["label"] = 0
    neg["neg_kind"] = np.where(neg["corum_status"] == "diff_complex", "strict", "random")

    # Complex-disjoint split for positives: assign whole complexes to cal/test ~50/50
    rng = np.random.default_rng(C.SEED)
    complexes = pos["complex_id"].dropna().unique().tolist()
    rng.shuffle(complexes)
    # greedily balance pair counts across splits
    counts = pos["complex_id"].value_counts().to_dict()
    cal_c, test_c, cal_n, test_n = set(), set(), 0, 0
    for cid in sorted(complexes, key=lambda c: -counts[c]):
        if cal_n <= test_n:
            cal_c.add(cid); cal_n += counts[cid]
        else:
            test_c.add(cid); test_n += counts[cid]
    pos["split"] = pos["complex_id"].map(lambda c: "cal" if c in cal_c else "test")

    # Negatives: seeded 50/50 split
    neg_split = rng.random(len(neg)) < 0.5
    neg["split"] = np.where(neg_split, "cal", "test")

    cols = ["uniprot_a", "uniprot_b", "label", "complex_id", "split",
            "score", "neg_kind", "corum_status", "gene_a", "gene_b", "is_true_pair"]
    labels = pd.concat([pos[cols], neg[cols]], ignore_index=True)
    dest = C.INTERIM / "labels.parquet"
    labels.to_parquet(dest, index=False)

    # Verification prints
    print("labels:", len(labels))
    print(pd.crosstab(labels["label"], labels["split"], margins=True))
    print("\nnegative kinds:", neg["neg_kind"].value_counts().to_dict())
    # leakage check: no complex spans cal & test
    span = pos.groupby("complex_id")["split"].nunique()
    assert (span <= 1).all(), "LEAK: a complex spans cal and test"
    # purity check: no positive appears as a negative
    posset = {frozenset((r.uniprot_a, r.uniprot_b)) for r in pos.itertuples()}
    negset = {frozenset((r.uniprot_a, r.uniprot_b)) for r in neg.itertuples()}
    assert not (posset & negset), "positive appears as negative"
    print(f"\ncomplex-disjoint split OK; pos/neg disjoint OK -> {dest}")


if __name__ == "__main__":
    main()
