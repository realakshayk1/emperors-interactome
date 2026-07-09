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
