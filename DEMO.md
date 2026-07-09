# DEMO — 3-minute video script (notebook + recorded narration)

Vehicle: run `notebooks/1.0-demo.ipynb` and narrate over it. Hard cap **3:00**. Required by the hackathon.
Goal: land the contrarian thesis, prove it with held-out evidence, end on a concrete nomination.
Judging leverage: Demo 30% (visual + live), plus it showcases Claude Use, Depth, Impact in one arc.

## Through-line (say this in one breath at the top)
"Published AI interactomes call thousands of protein complexes 'high-confidence.' I audited one with
distribution-free error control and let evidence the model never saw — DepMap co-essentiality — be the judge.
Here's what survives, and a cancer complex member it lets me nominate."

## Frame-by-frame (≈ 6 beats, ~30s each)

**0:00–0:25 — The claim.** Title card. One line: raw AlphaFold-Multimer confidence is overconfident and the
field's FDRs aren't distribution-free. Show the target map (CM4AI cell map (Schaffer et al., *Nature* 2025), U2OS cancer cells).

**0:25–0:55 — Frame 1: overconfidence.** Show `reliability_raw.png` — raw ipTM/pDockQ2 vs actual precision on
CORUM. Narrate: "well above the diagonal — the scores claim more certainty than they earn." (ECE number on screen.)

**0:55–1:30 — Frame 2: honest re-scoring.** Show `fdr_curve.png` + the certified-vs-dropped count at q=0.10.
Narrate the headline (Branch A: "N of the 'high-confidence' complexes don't survive honest FDR control" /
Branch B: "the certified core tightens to N with a guaranteed error rate"). Then `prevalence_shift.png`:
"and unlike benchmark-estimated FDR, this guarantee holds when real interactions are rare."

**1:30–2:05 — Frame 3: the referee agrees.** Show `heldout_enrichment.png`. Narrate: "the certified set is
[2–3]× more enriched for DepMap co-essentiality — cancer-cell-line dependency data the structure model never
saw. Honest error control isn't just conservative; it's more biologically real."

**2:05–2:45 — Frame 4: the nomination.** Name the target cancer complex. Show the nominee lighting up on
held-out co-essentiality, and (Stretch 1) its Boltz-2/AF3 predicted interface — `nominee_structure.png` — with
a physical-validity check. Narrate: "so I nominate [PROTEIN] as a missing member of [COMPLEX], with a certified
error rate and independent corroboration." (Optionally: one line of the Claude-Science literature rationale.)

**2:45–3:00 — Close.** "One map, honestly audited: fewer false complexes, a cleaner core, and a corroborated
new target — the reliability layer structural interactomics is missing." Card: repo link + "Built with Claude Science."

## Recording tips
- Pre-run the notebook so cells are cached; narrate over a clean scroll (don't wait on compute on camera).
- Put the key number on each figure (ECE, #dropped, enrichment fold, certified risk) as a title — readable at 720p.
- If Stretch 1 didn't finish, replace Frame 4's structure with the co-essentiality nomination alone — still lands.
- Keep energy high in the first 25s; judges decide fast. Lead with the contrarian claim, not the methods.

## Assets checklist (all in results/figures/)
- [ ] reliability_raw.png  - [ ] fdr_curve.png  - [ ] prevalence_shift.png
- [ ] heldout_enrichment.png  - [ ] nominee_structure.png (Stretch 1)  - [ ] (opt) predictomes_robustness.png
