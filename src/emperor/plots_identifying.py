"""Figure for the §A identifying experiment. Reads identifying_experiment.json
(run `python -m emperor.identifying` first) and writes identifying_experiment.png."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import config as C

FOCAL, COMPARE, LIGHT = "#1f6feb", "#9aa0a6", "#c7d7f0"


def main():
    src = C.PROCESSED / "identifying_experiment.json"
    o = json.load(open(src))
    qp = o["q_primary"]; qs = o["q_sweep"]; rows = o["rows"]
    d = np.array([r["delta"] for r in rows])
    pal = {0.05: COMPARE, 0.1: FOCAL, 0.2: LIGHT}

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for q in sorted(qs):
        y = np.array([r[f"fdr@{q}"] for r in rows]); se = np.array([r[f"se@{q}"] for r in rows])
        focal = abs(q - qp) < 1e-9
        ax.plot(d, y, "-o", color=pal.get(q, COMPARE), lw=2.4 if focal else 1.3,
                ms=6 if focal else 4, zorder=3 if focal else 2,
                label=f"realized FDR (q={q})" + (" — primary" if focal else ""))
        if focal:
            ax.fill_between(d, y - se, y + se, color=pal[q], alpha=0.2, zorder=1)
    ax.axhline(qp, ls="--", color="#444", lw=1.2)
    ax.text(0.0, qp + 0.015, f"nominal q = {qp}", fontsize=6.5, color="#444", ha="left")
    first_off = next(r["delta"] for r in rows if not r["control_holds@primary"])
    ax.axvspan(-0.05, first_off - 0.001, color="#5aa469", alpha=0.10, zorder=0)
    ax.text(0.0, 0.72, "exchangeable null\ncontrol ON\n(FDR \u2264 q)", fontsize=7, color="#3a7a48", ha="left", va="top")
    ax.text(1.28, 0.60, "non-exchangeable null\ncontrol OFF\n(FDR \u226b q)", fontsize=7, color="#b23a4b", ha="left")
    ax.set_xlabel("\u03b4  \u2014  shift of the true null toward the interaction signal (\u03c3)\n"
                  "(the ONLY variable toggled; everything else held fixed)")
    ax.set_ylabel("realized held-out FDR")
    ax.set_title("Distribution-free control holds if and only if the null is exchangeable\n"
                 "(\u00a7A identifying experiment: toggle exchangeability, hold all else fixed)", fontsize=8.5)
    ax.set_ylim(-0.02, 0.82); ax.set_xlim(-0.08, d.max() + 0.08)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right", bbox_to_anchor=(1.0, 0.0))
    fig.tight_layout()
    out = C.FIGURES / "identifying_experiment.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
