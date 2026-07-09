# The Emperor's Interactome — submission summary

Published AI interactomes label thousands of protein complexes "high-confidence," but their
error rates are estimated on balanced benchmarks and break when true interactions are rare. We
re-audit the CM4AI multimodal cell map (Schaffer et al., *Nature* 2025) with **distribution-free
conformal FDR control** and referee it with **DepMap co-essentiality** — evidence the structure
model never saw (strict purity firewall).

Raw AF-Multimer confidence is miscalibrated against held-out CORUM (ECE 0.18→0.02; AUROC 0.70),
so we locked the "artifacts likely" branch on Day 1. At q=0.10, **35 of 161 (22%) paper
high-confidence edges fail** honest FDR control. A prevalence-shift experiment shows a
benchmark-tuned cutoff's realized FDR climbing to 0.90 as interactions grow rare, while conformal
+ BH stays bounded near q. The independent referee agrees: certified edges are **41% co-essential
vs 17% for dropped edges** (permutation p=0.016). Finally, we nominate **KANSL3** as a missing
member of the leukemia-associated MLL1-WDR5 complex (certified risk 0.007; co-essential with the
NSL submodule KANSL1/MCRS1/PHF20) — independently confirmed as an NSL-complex subunit held out of
the pipeline. An independent Boltz-2 prediction recovers a localized KANSL3–KANSL1 interface
(ipTM 0.47; interface pLDDT 0.70), corroborating the certified edge. Fully reproducible via
`make reproduce`.

## Program extension

Five further analyses stress-test and generalize the audit. **(1) Overconfidence is structured:**
Mondrian (per-complex-size) conformal shows the pooled audit holds in medium complexes (held-out FDR
0.08) but breaks in large 9-member complexes (0.65, >3× the q=0.10 target). **(2) Covariate-shift
robustness is an honest negative:** the calibration→wild-interactome shift is real (KS 0.134,
p=1.4×10⁻⁹) and neither plain nor weighted (WCS) conformal survives it on node-disjoint splits — the
finite-sample guarantee does not transfer; we report the diagnosis rather than a rescue. **(3)
Baselines:** at the same nominal FDR=0.10, conformal-BH selects the fewest edges (132) yet has the
highest held-out DepMap support (41%), beating the fixed cutoff (204, 32%) and the published flag
(161, 37%). **(4) Auditability is a property of the data release:** of three public deposits, only
CM4AI — which shipped native decoys — is auditable; the Krogan host–pathogen deposit is positives-only
with a saturated pDockQ axis and Predictomes exposes only the DepMap-contaminated SPOC composite. The
recommendation for the field is concrete: ship the null. **(5) Cross-architecture pilot (conditional
GO):** on 12 pilot dimers, Boltz-2 corroborates the conformal verdict (certified ipTM 0.71 vs dropped
0.46, p=0.021), but the specific hypothesis that cross-architecture *divergence* predicts correctness
is only directional (one-sided p=0.090, n=12) — developed as a candidate feature, not integrated.
Finally, the single nomination generalizes to an **FDR-controlled selection**: 393 complexes receive
≥1 certified missing-member nomination (810 total, all inheriting FDR ≤ q=0.10, upper bound).
