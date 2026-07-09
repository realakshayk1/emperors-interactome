"""download.py — fetch raw datasets, checksum, write provenance rows.

VERIFIED firsthand 2026-07. Notable deviations from the original planning docs:
  * Krogan Suppl. Table 5 is inside the Suppl. Tables 1-10 .xlsx, fetched from
    Europe PMC (PMC12137143) because nature.com is not reachable from the sandbox.
  * CORUM moved to a FastAPI backend serving release 5.3 (not the old .txt.zip).
  * DepMap GLS matrices are 2.49 GB each — we only fetch genes.txt here; the big
    arrays are streamed + mmap-sliced in validate.py, never fully downloaded to RAM.
Run: python -m emperor.download
"""
from __future__ import annotations
import hashlib
import json
import ssl
import time
import urllib.request
import zipfile
from pathlib import Path

from . import config as C

UA = {"User-Agent": "emperors-interactome-audit"}
# The CORUM host (mips.helmholtz-muenchen.de) presented a broken TLS chain when
# fetched from this sandbox (urllib raised CERTIFICATE_VERIFY_FAILED, observed
# firsthand 2026-07). Verification is disabled for that single public CC-BY file
# only; the payload is size/schema-checked after download. See LEARNINGS.md.
_NOVERIFY = ssl.create_default_context()
_NOVERIFY.check_hostname = False
_NOVERIFY.verify_mode = ssl.CERT_NONE


def _fetch(url: str, dest: Path, verify: bool = True, timeout: int = 180) -> bytes:
    ctx = None if verify else _NOVERIFY
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        data = r.read()
    dest.write_bytes(data)
    return data


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _prov_row(name, fname, url, sha, size, license_, note=""):
    return dict(dataset=name, file=fname, url=url, sha256=sha,
                size_bytes=size, license=license_, date=time.strftime("%Y-%m-%d"),
                note=note)


def main() -> None:
    C.RAW.mkdir(parents=True, exist_ok=True)
    prov = []

    # 1) Krogan Suppl. Tables workbook (Europe PMC supplementary zip) -> extract xlsx
    supp_zip = C.RAW / "krogan_supp.zip"
    if not (C.RAW / C.KROGAN_XLSX_NAME).exists():
        data = _fetch(C.URLS["krogan_supp"], supp_zip)
        with zipfile.ZipFile(supp_zip) as zf:
            member = next(n for n in zf.namelist() if n.endswith(".xlsx"))
            zf.extract(member, C.RAW)
            (C.RAW / member).rename(C.RAW / C.KROGAN_XLSX_NAME)
    xb = (C.RAW / C.KROGAN_XLSX_NAME).read_bytes()
    prov.append(_prov_row("Krogan/Ideker Nature 2025", C.KROGAN_XLSX_NAME,
                          "PMC12137143 (Europe PMC supplementaryFiles)",
                          _sha256(xb), len(xb), "CC BY 4.0",
                          "Suppl. Tables 1-10 workbook; Table5 = primary interactome"))

    # 2) CORUM 5.3 human complexes (FastAPI backend; TLS-unverified host)
    corum = C.RAW / "corum_human.txt"
    if not corum.exists():
        _fetch(C.URLS["corum_human"], corum, verify=False)
    cb = corum.read_bytes()
    prov.append(_prov_row("CORUM 5.3 (human)", "corum_human.txt",
                          C.URLS["corum_human"], _sha256(cb), len(cb),
                          "CC BY 4.0", "release 5.3, 2026-04-14; TSV"))

    # 3) DepMap co-essentiality — genes.txt only (big GLS arrays streamed in validate.py)
    genes = C.RAW / "depmap_genes.txt"
    if not genes.exists():
        _fetch(C.URLS["depmap_genes"], genes)
    gb = genes.read_bytes()
    prov.append(_prov_row("DepMap co-essentiality (Wainberg 2021)", "depmap_genes.txt",
                          C.URLS["depmap_genes"], _sha256(gb), len(gb),
                          "academic-open",
                          "17,634 HGNC symbols; GLS_p/GLS_sign 2.49GB each, mmap-sliced on use"))

    (C.RAW / "provenance.json").write_text(json.dumps(prov, indent=2))
    for p in prov:
        print(f"{p['dataset']:38s} {p['size_bytes']:>10,d}b  {p['sha256'][:12]}  {p['file']}")
    print(f"\nprovenance -> {C.RAW/'provenance.json'}")


if __name__ == "__main__":
    main()
