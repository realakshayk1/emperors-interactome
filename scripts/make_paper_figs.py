#!/usr/bin/env python
"""Regenerate the four manuscript figures from the committed result artifacts.

Every value plotted is read from data/processed/*.json at draw time -- the same artifacts
scripts/verify_claims.py checks the prose against -- so a figure cannot drift away from
the number in the text without one of the two failing.

Design follows the convention in this literature: figures carry the *data*, captions carry
the *interpretation*. Concretely that means no panel titles stating a conclusion, no value
labels repeating what an axis already shows, and no inline statistics. A panel gets a
letter, axis labels, and at most a small frameless legend. Every effect size, p-value, n,
and takeaway belongs in the caption, where a reader can also see the qualifications.

Layout rules, because the previous hand-made figures broke all of them: constrained_layout
throughout (labels were being clipped mid-character); legends never inside a plot area
where they can land on the data; explicit ylim so nothing collides with the axes frame.
Series stay distinguishable in greyscale via marker and line style, not colour alone.

Figures are drawn at ~2x the printed width and included at \\linewidth, so font sizes here
are ~2x their on-page size. Run: python scripts/make_paper_figs.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
OUT_DIRS = [ROOT / "results" / "figures",
            ROOT / "paper" / "icbinb" / "figs",
            ROOT / "paper" / "gem" / "figs",
            ROOT / "paper" / "evalues" / "figs"]

BLUE, ORANGE, GREY = "#1f77b4", "#e8590c", "#7f7f7f"
LIGHT, DARK, PALE = "#9ecae1", "#08519c", "#d9d9d9"

plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 17,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 15,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 200,
    "savefig.dpi": 200,
})


def load(name: str):
    return json.loads((PROCESSED / name).read_text(encoding="utf-8"))


def panel_tags(fig, n: int):
    """Panel letters in a strip reserved above the axes.

    Two earlier attempts collided: axes-relative placement landed on long y-axis labels,
    and figure-relative placement landed on the axes themselves once the titles were
    removed and the plots grew to fill the top. Shrinking the layout rect first
    guarantees the strip is empty.
    """
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.93))
    xs = {2: (0.005, 0.515), 3: (0.005, 0.345, 0.675)}[n]
    for x, letter in zip(xs, "abc"):
        fig.text(x, 0.99, letter, fontsize=21, fontweight="bold", va="top", ha="left")


def save(fig, stem: str):
    for d in OUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / f"{stem}.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {stem}.png")


# --------------------------------------------------------------------------------------
def fig1_audit_wedge():
    bench, val = load("benchmark_synth.json"), load("validation.json")
    q = 0.1
    rows = sorted((r for r in bench["rows"] if r["q"] == q), key=lambda r: -r["prevalence"])
    pis = [r["prevalence"] for r in rows]
    x = list(range(len(pis)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)

    ax1.axhline(q, ls="--", color=GREY, lw=1.5, zorder=1)
    ax1.plot(x, [r["benchmark_fdr"] for r in rows], "s-", color=ORANGE, lw=2.5, ms=9,
             label="benchmark cutoff", zorder=3)
    ax1.plot(x, [r["conformal_fdr"] for r in rows], "o-", color=BLUE, lw=2.5, ms=9,
             label="conformal + BH", zorder=3)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{p:g}" for p in pis])
    ax1.set_xlabel("prevalence $\\pi$")
    ax1.set_ylabel("realized FDR")
    ax1.set_ylim(0, 1.0)
    ax1.set_xlim(-0.3, len(x) - 0.7)
    ax1.legend(frameon=False, loc="upper left", handlelength=1.8)

    sets = val["sets"]
    order = [("certified", "certified", BLUE), ("raw_high_conf", "all high-conf", GREY),
             ("dropped", "dropped", ORANGE)]
    ys = list(range(len(order)))[::-1]
    for y, (key, _lab, colour) in zip(ys, order):
        frac = sets[key]["frac_coess"]
        ax2.barh(y, frac, height=0.5, color=colour)
    ax2.set_yticks(ys)
    ax2.set_yticklabels([lab for _, lab, _ in order])
    ax2.set_xlim(0, 0.48)
    ax2.set_ylim(-0.55, len(order) - 0.45)
    ax2.set_xlabel("fraction co-essential")

    panel_tags(fig, 2)
    save(fig, "fig1_audit_wedge")


# --------------------------------------------------------------------------------------
def fig2_guarantee():
    sa, sen = load("shift_attribution.json"), load("sensitivity.json")
    q = "0.1"

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 4.2), constrained_layout=True)

    ks = sa["per_covariate_ks"]
    labels = [("degree", "endpoint degree"), ("iptm_mean", "ipTM"), ("score", "score"),
              ("iptm_ptm_gap", "ipTM$-$pTM gap"), ("ptm_mean", "pTM")]
    ys = list(range(len(labels)))[::-1]
    for y, (key, _lab) in zip(ys, labels):
        ax1.barh(y, ks[key]["ks"], height=0.5,
                 color=BLUE if key == "degree" else PALE)
    ax1.set_yticks(ys)
    ax1.set_yticklabels([lab for _, lab in labels])
    ax1.set_xlim(0, 0.28)
    ax1.set_ylim(-0.55, len(labels) - 0.45)
    ax1.set_xlabel("KS(decoy, candidate)")

    pop = sen["realized_fdr_population"]
    names = ["nominal", "identifying\ncurve", "protein-disjoint\nempirical"]
    vals = [float(q), pop["identifying_curve_fdr"][q], pop["empirical_node_disjoint_fdr"][q]]
    ax2.bar(names, vals, color=[PALE, LIGHT, DARK], width=0.6)
    ax2.axhline(float(q), ls="--", color=GREY, lw=1.5)
    ax2.set_ylabel("realized FDR")
    ax2.set_ylim(0, 0.36)

    tilt = sen["analytic_tilt"]
    gammas = tilt["gammas"]
    styles = (("0.2", DARK, "-"), ("0.1", BLUE, "--"), ("0.05", LIGHT, ":"))
    for qq, colour, ls in styles:
        ax3.plot(gammas, tilt["worst_case_fdr_curve"][qq], color=colour, ls=ls, lw=2.5,
                 label=f"q = {float(qq):g}")
    ax3.axvline(sen["realized_fdr_population"]["equivalent_gamma"],
                color=ORANGE, ls="-", lw=2, label="measured $\\Gamma$")
    ax3.set_xlabel("sensitivity $\\Gamma$")
    ax3.set_ylabel("worst-case FDR")
    ax3.set_ylim(0, 0.34)
    ax3.set_xlim(1, max(gammas))
    ax3.legend(frameon=False, loc="upper left", handlelength=2.2)

    panel_tags(fig, 3)
    save(fig, "fig2_guarantee")


# --------------------------------------------------------------------------------------
def fig3_validation():
    ppi, bench = load("experimental_ppi_referee.json"), load("benchmark_synth.json")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)

    groups = [("certified", BLUE), ("dropped", ORANGE)]
    width = 0.33
    for i, (key, colour) in enumerate(groups):
        r = ppi["results"][key]
        ax1.bar(i - width / 2, r["obs_rate"], width, color=colour,
                label="observed" if i == 0 else None)
        ax1.bar(i + width / 2, r["matched_null_mean"], width, color=PALE,
                yerr=1.96 * r["matched_null_sd"], capsize=5,
                error_kw=dict(ecolor=GREY, lw=1.5),
                label="degree-matched null" if i == 0 else None)
    ax1.set_xticks(range(len(groups)))
    ax1.set_xticklabels([k for k, _ in groups])
    ax1.set_ylabel("IntAct evidence")   # short: long labels run into the panel letter
    ax1.set_ylim(0, 0.95)
    ax1.legend(frameon=False, loc="upper right", handlelength=1.4)

    pis = sorted(bench["prevalences"], reverse=True)
    x = list(range(len(pis)))
    for qq, colour, mk, ls in ((0.2, DARK, "s", "-"), (0.1, BLUE, "o", "--"),
                               (0.05, LIGHT, "^", ":")):
        ys = [next(r["conformal_fdr"] for r in bench["rows"]
                   if r["q"] == qq and r["prevalence"] == p) for p in pis]
        ax2.axhline(qq, ls=":", color=colour, lw=1.2)
        ax2.plot(x, ys, marker=mk, ls=ls, color=colour, lw=2.5, ms=9, label=f"q = {qq:g}")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{p:g}" for p in pis])
    ax2.set_xlabel("prevalence $\\pi$")
    ax2.set_ylabel("realized FDR")
    ax2.set_ylim(0, 0.26)
    ax2.set_xlim(-0.3, len(x) - 0.7)
    ax2.legend(frameon=False, loc="upper left", ncol=3, handlelength=2.0,
               columnspacing=1.1)

    panel_tags(fig, 2)
    save(fig, "fig3_validation")


# --------------------------------------------------------------------------------------
def fig4_secondmap():
    sm = load("secondmap_audit.json")
    tier = sm["n_high_conf"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)

    qs = ["0.05", "0.1", "0.2"]
    cert = [100 * sm["certified_by_q"][q] / tier for q in qs]
    drop = [100 * sm["dropped_by_q"][q] / tier for q in qs]
    x = list(range(len(qs)))
    ax1.bar(x, cert, 0.55, color=BLUE, label="certified")
    ax1.bar(x, drop, 0.55, bottom=cert, color=ORANGE, label="dropped")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{float(q):g}" for q in qs])
    ax1.set_xlabel("target FDR $q$")
    ax1.set_ylabel("% of tier")
    ax1.set_ylim(0, 108)
    ax1.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)

    grad = sm["referee"]["referee_alive_spoc_gradient"]
    bins = sorted(grad, key=lambda b: float(b.split("-")[0]))
    x = list(range(len(bins)))
    ax2.plot(x, [grad[b]["coess"] for b in bins], "o-", color=BLUE, lw=2.5, ms=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels([b.replace("-", "–") for b in bins], rotation=30, ha="right")
    ax2.set_xlabel("SPOC bin")
    ax2.set_ylabel("fraction co-essential")
    ax2.set_ylim(0, 0.5)

    panel_tags(fig, 2)
    save(fig, "fig4_secondmap")


def main() -> int:
    print("regenerating manuscript figures from data/processed/ ...")
    fig1_audit_wedge()
    fig2_guarantee()
    fig3_validation()
    fig4_secondmap()
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
