import json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from emperor.plots import apply_figure_style, FOCAL, COMPARE, DROP
    apply_figure_style()
except Exception:
    FOCAL, COMPARE, DROP = "#1f6feb", "#9aa0a6", "#d1495b"

cm = pd.read_parquet("data/interim/interactome.parquet")
cand = cm[cm["is_true_pair"]]["score"].values
decoy = cm[~cm["is_true_pair"]]["score"].values
hp = np.array([r["pdockq"] for r in json.load(open("data/mimicry/apms_pdockq.json"))])

fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.2))

# Panel A — CM4AI ships a null
bins = np.linspace(0, 1, 41)
axa.hist(decoy, bins=bins, density=True, color=DROP, alpha=0.55, label=f"native random decoys (n={len(decoy)})")
axa.hist(cand, bins=bins, density=True, histtype="step", lw=2.0, color=FOCAL, label=f"candidate pairs (n={len(cand)})")
axa.set_title("Map 1 (CM4AI): ships an exchangeable null\n\u2192 distribution-free FDR audit is possible", fontsize=11)
axa.set_xlabel("AF-Multimer confidence score (ipTM)")
axa.set_ylabel("density")
axa.legend(fontsize=8, loc="upper right")
axa.set_xlim(0, 1)

# Panel B — host-pathogen positives-only, saturated
axb.hist(hp, bins=bins, density=True, color=COMPARE, alpha=0.7, label=f"predicted PPIs (n={len(hp)})")
axb.axvline(0.23, ls="--", color="0.4", lw=1.2)
axb.text(0.215, axb.get_ylim()[1]*0.95, "pDockQ=0.23\n(acceptable)", fontsize=7.5, color="0.35", ha="right", va="top")
axb.axvline(0.5, ls=":", color="0.55", lw=1.1)
axb.text(0.5, axb.get_ylim()[1]*0.55, " 0.5", fontsize=7.5, color="0.45", ha="left")
axb.set_title("Map 2 (host\u2013pathogen): positives-only, saturated\n\u2192 no null \u21d2 post-hoc FDR audit not possible", fontsize=11)
axb.set_xlabel("pDockQ (CPU, from ranked-0 model)")
axb.set_ylabel("density")
axb.annotate(f"93% \u2265 0.5; only 3.9% < 0.23\n75th\u201399th pctile span = 0.000",
             xy=(0.735, axb.get_ylim()[1]*0.55), xytext=(0.27, axb.get_ylim()[1]*0.42),
             fontsize=8, color=DROP,
             arrowprops=dict(arrowstyle="->", color=DROP, lw=1.2))
axb.set_xlim(0, 1)

fig.suptitle("Auditability is a property of the data release, not the map", fontsize=12.5, y=1.02)
fig.tight_layout()
fig.savefig("results/figures/auditability_boundary.png", dpi=150, bbox_inches="tight")
print("saved results/figures/auditability_boundary.png")
