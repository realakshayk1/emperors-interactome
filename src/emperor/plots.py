"""plots.py — all figures for the audit (reliability, FDR curve, prevalence
shift, held-out enrichment). Each figure is a function writing to results/figures/.
Run: python -m emperor.plots
"""
from __future__ import annotations
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as C

try:  # publication figure style (Claude Science skill helper), optional at runtime
    from figure_style import apply_figure_style  # type: ignore
    apply_figure_style()
except Exception:  # noqa: BLE001
    plt.rcParams.update({"figure.dpi": 140, "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})

FOCAL = "#1f6feb"     # certified / calibrated
COMPARE = "#9aa0a6"   # raw / comparator
DROP = "#d1495b"      # dropped / alarm


def fig_reliability():
    r = json.loads((C.PROCESSED / "calibration.json").read_text())
    raw = pd.DataFrame(r["reliability_raw"])
    iso = pd.DataFrame(r["reliability_isotonic"])
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.plot([0, 1], [0, 1], ls="--", lw=1, color="#888", zorder=0)
    ax.plot(raw["conf"], raw["acc"], "o-", color=COMPARE, label=f"raw AF-M score (ECE {r['ece_raw']:.2f})")
    ax.plot(iso["conf"], iso["acc"], "s-", color=FOCAL, label=f"isotonic-calibrated (ECE {r['ece_isotonic']:.2f})")
    ax.set_xlabel("predicted confidence")
    ax.set_ylabel("observed CORUM precision")
    ax.set_title("Raw AF-Multimer confidence is miscalibrated\nvs held-out CORUM labels")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.text(0.97, 0.03, "diagonal = perfect calibration", ha="right", va="bottom",
            fontsize=6.5, color="#666", transform=ax.transAxes)
    fig.tight_layout()
    dest = C.FIGURES / "reliability_raw.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dest


def fig_fdr_curve():
    s = json.loads((C.PROCESSED / "audit_summary.json").read_text())
    qs = sorted(float(q) for q in s["certified_by_q"])
    cert = [s["certified_by_q"][str(q)] if str(q) in s["certified_by_q"]
            else s["certified_by_q"][q] for q in qs]
    # keys may be str or float depending on json; normalize
    cb = {float(k): v for k, v in s["certified_by_q"].items()}
    db = {float(k): v for k, v in s["high_conf_dropped_by_q"].items()}
    qs = sorted(cb)
    cert = [cb[q] for q in qs]
    dropped = [db[q] for q in qs]
    n_hc = s["n_high_conf"]

    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    ax.plot(qs, cert, "o-", color=FOCAL, label="certified candidates")
    ax.plot(qs, dropped, "s-", color=DROP,
            label=f"paper 'high-conf' edges dropped\n(of {n_hc})")
    for q, c, d in zip(qs, cert, dropped):
        ax.annotate(f"{c}", (q, c), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=7, color=FOCAL)
        ax.annotate(f"{d}", (q, d), textcoords="offset points", xytext=(0, -12),
                    ha="center", fontsize=7, color=DROP)
    ax.set_xlabel("conformal FDR level  q")
    ax.set_ylabel("number of edges")
    ax.set_title("Honest error control drops a fifth of the\n'high-confidence' complexes at q=0.10")
    ax.set_xticks(qs)
    ax.legend(frameon=False, fontsize=7, loc="center right")
    ax.margins(0.08)
    fig.tight_layout()
    dest = C.FIGURES / "fdr_curve.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dest


def fig_prevalence_shift():
    d = json.loads((C.PROCESSED / "prevalence_shift.json").read_text())
    pv = d["prevalences"]
    bench = d["benchmark_fdr"]
    conf = d["conformal_fdr"]
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    ax.plot(pv, bench, "o-", color=DROP, label="fixed benchmark cutoff\n(realized error)")
    ax.plot(pv, conf, "s-", color=FOCAL, label="conformal + BH\n(realized error)")
    ax.axhline(d["q"], ls="--", lw=1, color="#888")
    ax.text(0.011, d["q"] + 0.02, f"claimed FDR q={d['q']}", va="bottom", ha="left",
            fontsize=6.5, color="#666")
    ax.annotate("conformal certifies nothing\nrather than exceed q",
                (pv[4], conf[4]), textcoords="offset points", xytext=(10, 34),
                fontsize=6, color=FOCAL, ha="left",
                arrowprops=dict(arrowstyle="-", color=FOCAL, lw=0.6))
    ax.set_xlabel("true interaction prevalence in test pool")
    ax.set_ylabel("realized FDR")
    ax.set_title("A benchmark-tuned cutoff breaks when interactions\nare rare; conformal control holds")
    ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.margins(y=0.12)
    fig.tight_layout()
    dest = C.FIGURES / "prevalence_shift.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dest


def fig_heldout_enrichment():
    v = json.loads((C.PROCESSED / "validation.json").read_text())
    order = ["dropped", "raw_high_conf", "certified"]
    labels = {"dropped": "dropped by\nconformal FDR", "raw_high_conf": "paper\n'high-confidence'",
              "certified": "conformally\ncertified"}
    colors = {"dropped": DROP, "raw_high_conf": COMPARE, "certified": FOCAL}
    fracs = [v["sets"][k]["frac_coess"] for k in order]
    ns = [v["sets"][k]["n_edges"] for k in order]
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    bars = ax.bar([labels[k] for k in order], fracs,
                  color=[colors[k] for k in order], width=0.62)
    for b, f, n in zip(bars, fracs, ns):
        ax.text(b.get_x() + b.get_width() / 2, f + 0.012, f"{f:.0%}\n(n={n})",
                ha="center", va="bottom", fontsize=7.5)
    ax.set_ylabel("fraction with held-out DepMap\nco-essentiality support")
    ax.set_title("Edges dropped by conformal FDR lack\nindependent co-essentiality support")
    ax.set_ylim(0, max(fracs) * 1.28)
    p = v.get("permutation_vs_dropped", {})
    if p:
        ax.text(0.5, 0.94, f"certified vs dropped: +{p['obs_frac_diff']:.0%},  "
                           f"permutation p = {p['p_frac']:.3f}",
                transform=ax.transAxes, ha="center", fontsize=7, color="#333")
    fig.tight_layout()
    dest = C.FIGURES / "heldout_enrichment.png"
    fig.savefig(dest, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dest


def main():
    for fn in [fig_reliability, fig_fdr_curve, fig_prevalence_shift, fig_heldout_enrichment]:
        try:
            print("wrote", fn())
        except FileNotFoundError as e:
            print(f"skip {fn.__name__}: {e}")


if __name__ == "__main__":
    main()
