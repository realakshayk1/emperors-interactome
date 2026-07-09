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
Narrate the headline: **"35 of the 161 'high-confidence' complexes — 22% — don't survive honest FDR
control."** Then `prevalence_shift.png`: "and unlike benchmark-estimated FDR, this guarantee holds when real
interactions are rare — the benchmark cutoff's error climbs to 90% as interactions get sparse; conformal stays near 10%."

**1:30–2:05 — Frame 3: the referee agrees.** Show `heldout_enrichment.png`. Narrate: "the certified set is
**41% co-essential vs 17% for the dropped edges** (permutation p=0.016) — cancer-cell-line dependency data the
structure model never saw. Honest error control isn't just conservative; it's more biologically real."

**2:05–2:45 — Frame 4: the nomination.** Target = leukemia-associated **MLL1-WDR5 complex (CORUM 5386)**. Show
the nominee lighting up on held-out co-essentiality, and its Boltz-2 predicted interface — `nominee_structure.png`
— with the physical-validity check. Narrate: "so I nominate **KANSL3** as a missing member of the MLL1-WDR5
complex, certified risk 0.007, independently confirmed as an NSL-complex subunit held out of the pipeline and
corroborated by a Boltz-2 interface." (Optional extension beat: the audit generalizes — 393 complexes get
FDR-controlled missing-member nominations; and auditability itself is a property of the data release.)

**2:45–3:00 — Close.** "One map, honestly audited: fewer false complexes, a cleaner core, and a corroborated
new target — the reliability layer structural interactomics is missing." Card: repo link + "Built with Claude Science."

## Recording tips
- Pre-run the notebook so cells are cached; narrate over a clean scroll (don't wait on compute on camera).
- Put the key number on each figure (ECE, #dropped, enrichment fold, certified risk) as a title — readable at 720p.
- If Stretch 1 didn't finish, replace Frame 4's structure with the co-essentiality nomination alone — still lands.
- Keep energy high in the first 25s; judges decide fast. Lead with the contrarian claim, not the methods.

## Assets checklist (all in results/figures/)
Core arc: [x] reliability_raw.png  [x] fdr_curve.png  [x] prevalence_shift.png
[x] heldout_enrichment.png  [x] nominee_structure.png
Extension (optional deeper-cut beats): [x] overconfidence_structured.png (Mondrian)
[x] shift_diagnostic.png + wcs_vs_plain_fdr.png (honest negative)  [x] baseline_headtohead.png
[x] auditability_boundary.png (ship-the-null)  [x] crossarch_pilot.png (Boltz-2 corroboration)
