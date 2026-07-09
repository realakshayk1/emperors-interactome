"""interactome.py — parse Krogan Suppl. Table 5 into the interactome data contract.

Output: data/interim/interactome.parquet with columns
  [pair, gene_a, gene_b, cluster, score, iptm_mean, ptm_mean, bench_fdr,
   is_true_pair, high_conf, high_conf_novel, uniprot_a, uniprot_b, source]

`score` is the AF-Multimer combined confidence (~mean 0.8*iptm+0.2*ptm across 5
models) and is the confidence axis for the audit (METHODS §2, single-metric path:
Table 5 has no pDockQ2 / tabulated interface-PAE). `is_true_pair`/random gives the
paper's native, same-pipeline, DepMap-independent negative set.
Run: python -m emperor.interactome
"""
from __future__ import annotations
import openpyxl
import pandas as pd

from . import config as C
from .idmap import symbols_to_acc


def load_table5() -> pd.DataFrame:
    xlsx = C.RAW / C.KROGAN_XLSX_NAME
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = wb[C.KROGAN_SHEET]
    rows = list(ws.iter_rows(values_only=True))
    header = rows[1]  # row 0 is the "Supplementary Table 5" title
    df = pd.DataFrame(rows[2:], columns=header)
    return df


def main() -> None:
    df = load_table5()
    iptm_cols = [f"iptm_{k}" for k in range(5)]
    ptm_cols = [f"ptm_{k}" for k in range(5)]
    out = pd.DataFrame({
        "pair": df["Pair"],
        "gene_a": df["geneA"].astype(str),
        "gene_b": df["geneB"].astype(str),
        "cluster": df["Cluster"],
        "score": pd.to_numeric(df[C.KROGAN_CONF_COL], errors="coerce"),
        "iptm_mean": df[iptm_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1),
        "ptm_mean": df[ptm_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1),
        "bench_fdr": pd.to_numeric(df[C.KROGAN_BENCHFDR_COL], errors="coerce"),
        "is_true_pair": df[C.KROGAN_TRUEFLAG_COL].astype(bool),
        "high_conf": df[C.KROGAN_HICONF_COL].astype(bool),
        "high_conf_novel": df[C.KROGAN_NOVEL_COL].astype(bool),
        "source": "krogan_nature2025_table5",
    })
    # Map gene symbols -> UniProt accessions (join key; METHODS/LEARNINGS: join on acc)
    genes = sorted(set(out["gene_a"]) | set(out["gene_b"]))
    acc = symbols_to_acc(genes)
    out["uniprot_a"] = out["gene_a"].map(acc)
    out["uniprot_b"] = out["gene_b"].map(acc)

    C.INTERIM.mkdir(parents=True, exist_ok=True)
    dest = C.INTERIM / "interactome.parquet"
    out.to_parquet(dest, index=False)

    n = len(out)
    mapped = out[["uniprot_a", "uniprot_b"]].notna().all(axis=1).sum()
    print(f"interactome: {n} pairs  ({out['is_true_pair'].sum()} true / "
          f"{(~out['is_true_pair']).sum()} random)")
    print(f"high_conf={out['high_conf'].sum()}  high_conf_novel={out['high_conf_novel'].sum()}")
    print(f"UniProt-mapped pairs: {mapped}/{n} ({100*mapped/n:.1f}%)")
    print(f"-> {dest}")


if __name__ == "__main__":
    main()
