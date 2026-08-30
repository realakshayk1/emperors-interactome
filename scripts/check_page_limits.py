"""Measure where each manuscript's body actually ends.

The obvious check -- find the page carrying "References" and call the one before it the last
body page -- is wrong, because the bibliography usually starts partway down a page that still
carries body text above it. That check reported "body ends p4" for a paper whose Discussion and
Limitations were sitting on page 5.

This locates the References heading within the page's text and asks whether anything precedes
it, which is the question the page limit actually turns on.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# body page limits, excluding references and appendices
LIMITS = {"evalues": 4, "icbinb": 8, "gem": 5}


def _page_text(pdf: Path, page: int) -> str:
    return subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), str(pdf), "-"],
                          capture_output=True, text=True).stdout


def _n_pages(pdf: Path) -> int:
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    return int(out.split("Pages:")[1].split()[0])


def main() -> int:
    fails = []
    for name, limit in LIMITS.items():
        pdf = ROOT / "paper" / name / "main.pdf"
        if not pdf.exists():
            print(f"  {name}: no main.pdf")
            continue
        n = _n_pages(pdf)
        last_body, split_page = None, None
        for pg in range(1, n + 1):
            txt = _page_text(pdf, pg)
            idx = txt.find("References")
            if idx < 0:
                last_body = pg
                continue
            # strip the lineno margin numbers before judging what precedes the heading
            before = "".join(c for c in txt[:idx] if not c.isdigit()).strip()
            if len(before) > 40:
                last_body, split_page = pg, pg
            break
        status = "ok" if last_body and last_body <= limit else "OVER"
        note = f" (body shares p{split_page} with References)" if split_page else ""
        print(f"  {name:9s} {n} pages total | body ends p{last_body} | limit {limit}  {status}{note}")
        if not last_body or last_body > limit:
            fails.append(f"{name} body ends on p{last_body}, limit {limit}")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("all bodies within their page limits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
