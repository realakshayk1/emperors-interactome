"""idmap.py — gene-symbol <-> UniProt-accession crosswalk.

Strategy (LEARNINGS: join through UniProt accessions, never raw symbols):
  1. Local, no-network: harvest symbol<->accession pairs from CORUM's dual
     `subunits_gene_name` / `subunits_uniprot_id` columns.
  2. Network fallback: for symbols still unmapped, query the UniProt REST API
     (reviewed, human) one batch at a time.
Result cached to data/interim/idmap.parquet as columns [symbol, accession].
Public functions: symbols_to_acc(symbols) -> {symbol: accession}.
Run: python -m emperor.idmap
"""
from __future__ import annotations
import json
import time
import urllib.parse
import urllib.request

import pandas as pd

from . import config as C

UA = {"User-Agent": "emperors-interactome-audit", "Accept": "text/plain"}
_CACHE = C.INTERIM / "idmap.parquet"


def _corum_symbol_acc() -> dict[str, str]:
    """Build symbol->accession from CORUM dual columns (first-seen wins)."""
    df = pd.read_csv(C.RAW / "corum_human.txt", sep="\t", dtype=str)
    m: dict[str, str] = {}
    for accs, names in zip(df["subunits_uniprot_id"], df["subunits_gene_name"]):
        if not isinstance(accs, str) or not isinstance(names, str):
            continue
        for a, n in zip(accs.split(";"), names.split(";")):
            a, n = a.strip(), n.strip()
            if a and n and n not in m:
                m[n] = a
    return m


def _uniprot_lookup(symbols: list[str]) -> dict[str, str]:
    """Query UniProt REST for reviewed human accession of each gene symbol."""
    out: dict[str, str] = {}
    for i in range(0, len(symbols), 20):
        batch = symbols[i:i + 20]
        clause = "+OR+".join(f"gene_exact:{urllib.parse.quote(s)}" for s in batch)
        q = ("https://rest.uniprot.org/uniprotkb/search?query="
             f"({clause})+AND+organism_id:9606+AND+reviewed:true"
             "&fields=accession,gene_primary&format=tsv&size=500")
        try:
            req = urllib.request.Request(q, headers=UA)
            txt = urllib.request.urlopen(req, timeout=60).read().decode()
        except Exception as e:  # noqa: BLE001
            print(f"  [uniprot] batch {i} failed: {e!r}")
            continue
        for line in txt.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1]:
                acc, prim = parts[0], parts[1]
                if prim not in out:
                    out[prim] = acc
        time.sleep(0.2)
    # Second pass: symbols still unmapped may be aliases/renamed -> search the full
    # `gene` field (synonyms), keyed back to the ORIGINAL query symbol.
    still = [s for s in symbols if s not in out]
    for s in still:
        head = s.split("-")[0]  # readthrough like PMF1-BGLAP -> PMF1
        q = ("https://rest.uniprot.org/uniprotkb/search?query="
             f"(gene:{urllib.parse.quote(head)})+AND+organism_id:9606+AND+reviewed:true"
             "&fields=accession,gene_names&format=tsv&size=5")
        try:
            req = urllib.request.Request(q, headers=UA)
            txt = urllib.request.urlopen(req, timeout=60).read().decode()
        except Exception:  # noqa: BLE001
            continue
        for line in txt.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and (head in parts[1].split() or s in parts[1].split()):
                out[s] = parts[0]
                break
        time.sleep(0.2)
    return out


def build(symbols: list[str] | None = None) -> pd.DataFrame:
    """Build/refresh the crosswalk covering `symbols` (defaults to CORUM only)."""
    mapping = _corum_symbol_acc()
    if symbols:
        missing = sorted(set(symbols) - set(mapping))
        if missing:
            print(f"  [idmap] {len(missing)} symbols not in CORUM -> UniProt REST")
            mapping.update(_uniprot_lookup(missing))
    df = pd.DataFrame(sorted(mapping.items()), columns=["symbol", "accession"])
    C.INTERIM.mkdir(parents=True, exist_ok=True)
    df.to_parquet(_CACHE, index=False)
    return df


def _load_cache() -> dict[str, str]:
    if _CACHE.exists():
        df = pd.read_parquet(_CACHE)
        return dict(zip(df["symbol"], df["accession"]))
    return {}


def symbols_to_acc(symbols: list[str]) -> dict[str, str]:
    """Map gene symbols -> UniProt accessions, extending + caching as needed."""
    cache = _load_cache()
    missing = sorted(set(symbols) - set(cache))
    if missing:
        build(symbols)  # rebuild covering everything requested
        cache = _load_cache()
    return {s: cache[s] for s in symbols if s in cache}


def main() -> None:
    inter = C.INTERIM / "interactome.parquet"
    symbols = None
    if inter.exists():
        idf = pd.read_parquet(inter)
        symbols = sorted(set(idf["gene_a"]) | set(idf["gene_b"]))
    df = build(symbols)
    print(f"crosswalk: {len(df)} symbol->accession entries -> {_CACHE}")
    if symbols:
        cov = sum(s in set(df["symbol"]) for s in symbols)
        print(f"interactome symbol coverage: {cov}/{len(symbols)} ({100*cov/len(symbols):.1f}%)")


if __name__ == "__main__":
    main()
