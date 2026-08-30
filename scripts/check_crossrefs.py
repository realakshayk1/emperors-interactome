"""Find labels nothing points at, and references that point at the wrong kind of thing.

LaTeX warns about undefined references but says nothing about a defined label with no \\ref --
so a figure can be orphaned by an edit that deletes the only sentence citing it, and the build
stays clean. That happened here: collapsing a section left the shift figure unreferenced, and
left two appendix pointers aimed at appendices that no longer held the material.

The kind check is a heuristic on the label prefix: fig: should be cited as "Figure~\\ref",
tab: as "Table~\\ref", app: and sec: as "Section~\\ref" or "Appendix~\\ref". It catches a
figure cited as a section, but NOT a pointer aimed at the wrong appendix -- that reference is
kind-correct and only a reader can see it is aimed at the wrong place. The orphan arm is the
one with real power here.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPERS = ("evalues", "icbinb", "gem")

KIND = {
    "fig": ("Figure",),
    "tab": ("Table",),
    "app": ("Appendix", "Section"),
    "sec": ("Section", "Appendix"),
    "eq": ("Equation", "Eq"),
    "prop": ("Proposition",),
}


def main() -> int:
    fails = []
    for name in PAPERS:
        src = ROOT / "paper" / name / "main.tex"
        if not src.exists():
            continue
        t = re.sub(r"(?m)^\s*%.*$", "", src.read_text(encoding="utf-8"))

        labels = set(re.findall(r"\\label\{([^}]+)\}", t))
        refs = re.findall(r"(\w+)~?\\(?:ref|eqref)\{([^}]+)\}", t)
        refd = {r[1] for r in refs}

        orphans = sorted(labels - refd)
        wrong = []
        for word, lab in refs:
            pre = lab.split(":")[0]
            if word.lower() in {"and", "or", "to", "in", "of", "see", "from"}:
                continue  # second element of a list: "Propositions 1 and 3"
            allowed = KIND.get(pre)
            if allowed and word.rstrip("s") not in tuple(a.rstrip("s") for a in allowed):
                wrong.append(f"{word}~\\ref{{{lab}}}")

        print(f"  {name:9s} {len(labels):2d} labels, {len(refd):2d} referenced")
        if orphans:
            print(f"    never referenced: {', '.join(orphans)}")
            fails.append(f"{name}: orphaned {', '.join(orphans)}")
        if wrong:
            print(f"    wrong kind: {', '.join(sorted(set(wrong)))}")
            fails.append(f"{name}: mis-typed reference {sorted(set(wrong))[0]}")

    print()
    if fails:
        print("FAIL: " + "; ".join(fails))
        return 1
    print("every label is referenced and every reference names the right kind")
    return 0


if __name__ == "__main__":
    sys.exit(main())
