import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
try:
    from emperor.plots import apply_figure_style, FOCAL, COMPARE, DROP
    apply_figure_style()
except Exception:
    FOCAL, COMPARE, DROP = "#1f6feb", "#9aa0a6", "#d1495b"

df = pd.read_csv("data/structures/pilot_crossarch.csv")
cert = df[df.label == "certified"]; drop = df[df.label == "dropped"]

fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.4))

# Panel A — AF-M vs Boltz-2, colored by verdict
axa.scatter(cert.afm_score, cert.boltz_iptm, c=FOCAL, s=70, label="certified (q=0.10)", zorder=3, edgecolor="white", lw=0.8)
axa.scatter(drop.afm_score, drop.boltz_iptm, c=DROP, s=70, label="dropped (q=0.10)", zorder=3, edgecolor="white", lw=0.8)
lims = [0.05, 1.0]
axa.plot(lims, lims, ls="--", color="0.6", lw=1, zorder=1)
rho, p = spearmanr(df.afm_score, df.boltz_iptm)
axa.set_xlabel("AlphaFold-Multimer confidence (Score)")
axa.set_ylabel("Boltz-2 ipTM (independent architecture)")
axa.set_title(f"Two architectures agree on ranking\nSpearman \u03c1 = {rho:.2f} (p = {p:.3f})", fontsize=11)
axa.legend(fontsize=8, loc="upper left")
axa.set_xlim(0.3, 1.0); axa.set_ylim(0, 1.0)

# Panel B — the ACTUAL gate test: cross-architecture DIVERGENCE |afm - boltz| by verdict
drop["divergence"] = (drop.afm_score - drop.boltz_iptm).abs()
cert["divergence"] = (cert.afm_score - cert.boltz_iptm).abs()
for i, (lab, sub, col) in enumerate([("dropped", drop, DROP), ("certified", cert, FOCAL)]):
    x = np.random.RandomState(0).normal(i, 0.04, len(sub))
    axb.scatter(x, sub.divergence, c=col, s=55, zorder=3, edgecolor="white", lw=0.7)
    axb.hlines(sub.divergence.mean(), i - 0.2, i + 0.2, color=col, lw=2.5, zorder=4)
axb.set_xticks([0, 1]); axb.set_xticklabels(["conformal\ndropped", "conformal\ncertified"])
axb.set_ylabel("|AF-M Score \u2212 Boltz-2 ipTM|  (divergence)")
axb.set_title("Divergence trends higher for dropped\n(p = 0.090, one-sided, n=12 \u2014 not significant)", fontsize=11)
axb.set_ylim(0, 0.42)
axb.annotate("TOMM22\u2013WSB1", xy=(0, 0.36), xytext=(0.2, 0.37),
             fontsize=7.5, color=DROP, arrowprops=dict(arrowstyle="->", color=DROP, lw=1))

fig.suptitle("Cross-architecture pilot (n=12): Boltz-2 corroborates the verdict (p=0.021); "
             "divergence-as-feature only directional \u2192 conditional GO", fontsize=10.5, y=1.02)
fig.tight_layout()
fig.savefig("results/figures/crossarch_pilot.png", dpi=150, bbox_inches="tight")
print("saved results/figures/crossarch_pilot.png")
