"""Check that the three manuscripts agree on the quantities they share.

They audit one dataset and quote overlapping numbers, so a correction applied to one paper and
not the others leaves a contradiction that no single-paper check can see. That happened: the
tier's share of the calibration-to-candidate shift is 91.65%, which the e-values paper rounds
to 92% and the other two rounded to 91% -- the same sentence, three papers, two answers.

This reads the .tex sources rather than the artifacts, because the failure mode is a number
that reproduces perfectly against its artifact in each paper separately while the papers
disagree with each other.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPERS = ("evalues", "icbinb", "gem")

# A shared claim: a regex that captures the number, and the value all papers must agree on.
# Every entry must capture the value, or the comparison is between the string "present" and
# itself. A paper missing the quantity is a finding, not something to skip.
SHARED = {
    "tier share of the 0.60 sigma shift": (r"\$(9\d)\\%\$ of (?:it|that|that gap) comes from the 161", "92"),
    "measured shift (sigma)": (r"(0\.60)\\sigma|\$(0\.60)\$ standard deviations", None),
    "degree KS": (r"degree[^.]{0,40}?KS \$?(0\.25\d)", None),
    "score KS": (r"score[^.]{0,30}?(0\.12\d)", None),
    "certified at q=0.10": (r"certifies (13\d) of the|\b(13\d) of 1\{,\}666 candidates certify", None),
}


def main() -> int:
    fails = []
    for label, (pat, expect) in SHARED.items():
        found = {}
        for name in PAPERS:
            src = ROOT / "paper" / name / "main.tex"
            if not src.exists():
                continue
            t = src.read_text(encoding="utf-8")
            m = re.search(pat, t)
            if not m:
                found[name] = "ABSENT"
                continue
            groups = [g for g in m.groups() if g is not None]
            found[name] = groups[0] if groups else "present"
        if not found:
            continue
        vals = set(found.values())
        # Three ABSENTs agree with each other, so equality alone would pass a quantity that
        # has vanished from every paper -- which is how a shared number gets silently dropped.
        if "ABSENT" in vals:
            status = "MISSING"
        elif len(vals) == 1:
            status = "ok"
        else:
            status = "DISAGREE"
        detail = ", ".join(f"{k}={v}" for k, v in found.items())
        print(f"  {label:36s} {status:9s} {detail}")
        if len(vals) > 1 or "ABSENT" in vals:
            fails.append(f"{label}: {detail}")
        if expect and vals and next(iter(vals)) != expect and len(vals) == 1:
            fails.append(f"{label}: all papers say {next(iter(vals))}, expected {expect}")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("the manuscripts agree on every shared quantity checked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
