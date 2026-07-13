# DEMO — walkthrough script

A short narrated walkthrough of `notebooks/1.0-demo.ipynb`. One thesis, proven with held-out
evidence, ending on the reliability standard. Target ~3 minutes.

## Through-line (say this at the top)

"The CM4AI cell map — a Nature 2025 AI interactome from the Bridge2AI consortium co-led by Nevan Krogan (Gladstone) and Trey Ideker — is
one of the few you can actually audit, because it ships a native random-pair null. I audited its
'high-confidence' complexes with distribution-free error control and let evidence the model never
saw — DepMap co-essentiality — be the judge. Here's what survives, and why it matters."

## Beats (~5, ~30s each)

**0:00–0:25 — The claim.** Title card. One line: a "high-confidence" cutoff tuned on a balanced
benchmark hides a much larger false-discovery rate in the real regime, where true interactions are
rare. Name the map (CM4AI, Schaffer et al. / Ideker lab, Bridge2AI Krogan/Ideker consortium, *Nature* 2025).

**0:25–1:00 — The prevalence wedge.** Show `fig1_audit_wedge.png` (panel A). Narrate: "on
known-truth data, a fixed benchmark cutoff's realized FDR climbs toward 0.9 as interactions grow
rare — while conformal + BH stays bounded at the target q. That gap is the problem." This is the
single headline motivation.

**1:00–1:35 — The audit.** Show the certified-vs-dropped count at q=0.10. Narrate the headline:
"**35 of 161 published high-confidence edges — 22% — don't survive honest FDR control.** The audit
has teeth precisely because this map ships the null that makes distribution-free control possible."

**1:35–2:15 — The referee agrees.** Show `fig1_audit_wedge.png` (panel B) / `heldout_enrichment.png`.
Narrate: "the certified set is **41% co-essential versus 17% for the dropped edges** (p=0.016) —
cancer-cell dependency data the structure model never saw. Honest error control isn't just
conservative; the surviving edges are more biologically real."

**2:15–2:45 — Reliability, made precise.** One line that the guarantee was stress-tested: the
strongly-certified core stays controlled under exchangeability violations far larger than the one
actually present, and it is independently enriched for wet-lab physical interactions. (Depth cue;
don't belabor the statistics on camera.)

**2:45–3:00 — Close.** "One map, honestly audited: fewer false complexes, a cleaner core, and a
concrete reliability standard for AI interactomes — ship a raw confidence axis and a native null,
so your map can be audited too." Card: repo link.

## Recording tips
- Pre-run the notebook so cells are cached; narrate over a clean scroll.
- Put the key number on each figure (wedge 0.9 vs q; 35/161; 41% vs 17%) as a title — readable at 720p.
- Lead with the wedge and the map lineage, not the methods. One headline number: 22%. One positive
  result: 41% vs 17% co-essential.

## Assets (results/figures/)
- [ ] fig1_audit_wedge.png (prevalence wedge + co-essentiality)
- [ ] fig2_guarantee.png (reliability stress-test — supporting)
- [ ] fig3_validation.png (independent physical validation — supporting)
- [ ] heldout_enrichment.png (referee enrichment)
