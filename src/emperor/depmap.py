"""depmap.py — RAM-safe slice of the DepMap co-essentiality GLS matrix.

The GLS_p / GLS_sign arrays are 17,634 x 17,634 float64 (2.49 GB each). We never
download or materialize them fully (LEARNINGS: OOM risk on a small machine).
Instead we HTTP-range-fetch only the columns for interactome genes. The .npy is
Fortran-ordered, so column j (all rows) is contiguous on disk — one range read
per needed gene -> a ~1,049 x 1,049 submatrix (a few MB) into RAM.

Output: data/interim/depmap_slice.npz with keys [genes, gls_p, gls_sign].
DepMap is the HELD-OUT referee: it appears ONLY here and in validate.py, never in
labels/nonconformity/calibration. Run: python -m emperor.depmap
"""
from __future__ import annotations
import struct
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from . import config as C

_UA = {"User-Agent": "emperors-interactome-audit"}
N = 17634  # verified matrix dimension
ITEM = 8   # float64


def _npy_data_offset(url: str) -> tuple[int, tuple[int, int], bool]:
    """Parse the .npy header via a range read; return (data_offset, shape, fortran)."""
    req = urllib.request.Request(url, headers={**_UA, "Range": "bytes=0-255"})
    hdr = urllib.request.urlopen(req, timeout=60).read()
    assert hdr[:6] == b"\x93NUMPY", "not a .npy"
    hlen = struct.unpack("<H", hdr[8:10])[0]
    meta = eval(hdr[10:10 + hlen].decode("latin1"))  # trusted numpy header dict
    return 10 + hlen, meta["shape"], meta["fortran_order"]


def _fetch_columns(url: str, offset: int, cols: list[int], nrows: int,
                   workers: int = 16) -> np.ndarray:
    """Fetch full columns (Fortran-order: column-contiguous) -> (nrows, len(cols)).
    Parallelized across `workers` threads (I/O-bound HTTP range reads)."""
    out = np.empty((nrows, len(cols)), dtype="<f8")
    colbytes = nrows * ITEM

    def one(k_j):
        k, j = k_j
        start = offset + j * colbytes
        end = start + colbytes - 1
        req = urllib.request.Request(url, headers={**_UA, "Range": f"bytes={start}-{end}"})
        raw = urllib.request.urlopen(req, timeout=120).read()
        out[:, k] = np.frombuffer(raw, dtype="<f8", count=nrows)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one, enumerate(cols)))
    return out


def run():
    genes = [l.strip() for l in open(C.RAW / "depmap_genes.txt")]
    gidx = {g: i for i, g in enumerate(genes)}
    inter = pd.read_parquet(C.INTERIM / "interactome.parquet")
    inter_genes = sorted(set(inter.gene_a) | set(inter.gene_b))
    keep_genes = [g for g in inter_genes if g in gidx]
    cols = [gidx[g] for g in keep_genes]

    p_url, s_url = C.URLS["depmap_gls_p"], C.URLS["depmap_gls_sign"]
    off_p, shape, fort = _npy_data_offset(p_url)
    assert shape == (N, N) and fort, f"unexpected header {shape} fortran={fort}"
    off_s, _, _ = _npy_data_offset(s_url)

    # Fetch full columns, then subset to the same rows -> square submatrix
    p_cols = _fetch_columns(p_url, off_p, cols, N)          # (N, k)
    s_cols = _fetch_columns(s_url, off_s, cols, N)
    rows = np.array(cols)
    gls_p = p_cols[rows, :]                                  # (k, k)
    gls_sign = s_cols[rows, :]

    dest = C.INTERIM / "depmap_slice.npz"
    np.savez_compressed(dest, genes=np.array(keep_genes), gls_p=gls_p, gls_sign=gls_sign)
    print(f"DepMap slice: {len(keep_genes)} genes ({len(inter_genes)} interactome, "
          f"{len(inter_genes)-len(keep_genes)} not in DepMap)")
    print(f"submatrix {gls_p.shape}  p-range[{np.nanmin(gls_p):.2g},{np.nanmax(gls_p):.2g}] -> {dest}")
    return dest


if __name__ == "__main__":
    run()
