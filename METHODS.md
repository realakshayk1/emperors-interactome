# METHODS — The scientific core

This is the load-bearing document. It specifies the nonconformity score, the conformal FDR
procedure, the physical-validity add-on, the held-out validation design, and the nomination.
Cite the machinery; the novelty is the *first structural-PPI-confidence instantiation + held-out audit*.

---

## 1. Problem setup
We have candidate protein pairs (edges) each with structural confidence metrics from AlphaFold-Multimer.
We want a **distribution-free, FDR-controlled** subset of edges we certify as true interactions, plus a
**calibrated probability** per edge, and an honest demonstration that the certified set reflects real
biology using a channel the model never saw.

Two prediction targets, both handled:
- **Calibration (probabilities):** map raw confidence → P(true interaction). Diagnostic; shows overconfidence.
- **Selection (FDR control):** output a set S of certified edges with FDR ≤ q. The headline result.

---

## 2. Nonconformity score
For edge i, define a scalar nonconformity score `s_i` where **lower = more confident it's a true interaction**.
Base it on raw AF-M interface confidence plus a physical-validity term:

```
s_i = w1 * (1 - iptm_i)
    + w2 * (1 - pdockq2_i)
    + w3 * normalized_interface_PAE_i
    + w4 * phys_penalty_i
```

- Start with equal weights or, better, fit a small monotone/logistic model on the calibration split so
  `s_i` is a learned "improbability of true interaction." Keep it simple and seeded.
- If Table 5 has only one metric (e.g., only pDockQ2), use it alone + the physical-validity term. Don't fabricate metrics.
- `phys_penalty_i` comes from §4.

**Design rule:** `s_i` may use ONLY structural/AF quantities and physical validity. It must NOT use DepMap,
co-expression, GO, or any channel reserved for validation (§5). This is the held-out-purity firewall.

---

## 3. Calibration labels (positives / negatives)
- **Positives:** all within-complex protein pairs from CORUM human core complexes (VERIFIED release 5.3,
  2026-04-14; see DATA.md §2 for the version/host correction) — explode each complex's
  semicolon-separated subunit list into all pairwise combinations). Optionally union with hu.MAP 3.0 high-confidence complexes.
- **Negatives (decoys):** protein pairs not co-occurring in any known complex. Construction matters:
  - Start with **random non-co-complex pairs** matched on degree/abundance where possible, ratio 1:1.
  - Report sensitivity to the positive:negative ratio (1:1, 1:5) — the PPI negative-set problem is real; be explicit.
  - A stricter negative set: pairs whose proteins are each in *some* complex but never the *same* one.
- **Split:** partition labels into a **calibration** set (used to compute conformal p-values) and a held-out
  **CORUM test** set (used only to check empirical FDR ≤ q). Keep complexes disjoint across splits (no leakage of a complex's pairs across cal/test).

---

## 4. Physical-validity add-on (your signature contribution)
Reuse the physical-validity reflex from your `flow-match`/`calibrated-dock` work, transposed to interfaces.
From each deposited AF-M complex structure (the CM4AI map ships 111 high-confidence models; for others
recompute lightly or use interface metrics already tabulated), compute an interface-quality vector:

- **Steric clashes** at the interface (count of atom pairs below VdW threshold; Biopython/`freesasa`).
- **Buried surface area (BSA)** on complex formation — implausibly small BSA ⇒ likely spurious interface.
- **Interface pLDDT** (mean pLDDT of interface residues) and **interface PAE** (mean predicted aligned error
  across the interface) — low-confidence interfaces should be penalized.
- Optional: a strain/energy proxy on interface side chains.

Combine into `phys_penalty_i ∈ [0,1]` (clashes ↑, tiny BSA ↑, high interface-PAE ↑ → larger penalty). This term
is what lets a "high-ipTM but physically impossible" edge be correctly demoted — the crux of the audit.

Tooling: Biopython (`Bio.PDB`), `freesasa` for BSA, numpy for PAE aggregation. PoseBusters is ligand-oriented;
for protein–protein use the clash+BSA+interface-PAE composite above. Metric computation is CPU.

**BioNeMo upgrade (Modal-backed):** if the CM4AI deposited models or their tabulated interface metrics are
incomplete, *generate* the interface structures with the **Boltz-2 / OpenFold3 / AlphaFold2-Multimer BioNeMo
skills** and compute clash/BSA/interface-PAE from those. Budget-limited ($30 Modal ≈ few hundred predictions) —
apply only to the high-leverage subset (the ~111 high-confidence complex interfaces + decision-boundary edges),
NOT the whole map. See DATA.md §5 for the budget. This turns the physical-validity term from a fallback into a
first-class, real-structure signal.

---

## 5. Conformal FDR control (the headline method)
Use **conformal p-values + Benjamini–Hochberg**, following conformal selection / conformal link-prediction FDR
(Marandon 2023, arXiv:2306.14693; Blanchard et al. 2024, arXiv:2404.02542; conformal selection, Jin & Candès 2023).

For each candidate edge j (from the interactome, the "test" pool), and using the calibration **negatives**
{s_k : k in cal, label=0} (size n), compute the conformal p-value:

```
p_j = ( 1 + #{ k in cal_negatives : s_k <= s_j } ) / ( n + 1 )
```

Intuition: p_j is small when edge j looks *more like a true interaction* than the calibration negatives do
(its nonconformity is lower than most negatives). Under exchangeability, p_j is a valid p-value for the null
"j is a non-interaction."

Then apply **Benjamini–Hochberg** to {p_j} at level q → the rejected set S is the **certified interactome**
with FDR ≤ q (distribution-free, marginal). Report:
- the certified set size at q=0.10, plus a sweep q∈{0.05,0.10,0.20};
- how many previously "high-confidence" edges are **dropped** (Branch A headline);
- a comparison to the interactome's own published cutoff / SPOC benchmark-FDR (show ours holds under prevalence shift, theirs doesn't).

**Prevalence-shift centerpiece (the wedge vs SPOC):** benchmark-estimated FDR assumes the test prevalence
equals the benchmark's positive:negative ratio. Re-weight the test pool to a realistic (low) interaction
prevalence and show the benchmark-FDR undershoots true errors while the conformal FDR guarantee still holds
on the held-out CORUM test split. This is the strongest "why conformal, not SPOC" evidence.

**Validity check (unit test):** on synthetic data with known nulls, verify empirical FDR ≤ q across seeds.

---

## 6. Held-out validation design (the referee)
The certified set is only convincing if an *independent* signal agrees. Use **DepMap co-essentiality**
(Wainberg 2021 GLS p-value/sign matrix). Purity firewall: DepMap is used ONLY here, never in §2–§5.

For each edge set (certified@q, raw-high-confidence, dropped-artifacts), compute co-essentiality strength
for its protein pairs (look up GLS value; higher |GLS| / significant = stronger co-essentiality). Then:
- **Enrichment metric:** mean co-essentiality, and AUROC/odds-ratio of "is a strong co-essential pair" for
  certified vs raw vs dropped.
- **Permutation p-value:** shuffle edge-set membership to get a null; report p.
- **Expected result:** certified > raw high-confidence > dropped (certified 2–3× enriched). If the map is
  well-calibrated (Branch B), certified ≈ raw but both > dropped, and the *certified core* is the cleanest —
  still a positive, quotable result.

Purity checklist (put in the notebook): DepMap not in `nonconformity.py`; primary interactome is AP-MS+AF
(DepMap-independent); on Predictomes we audit raw pDockQ2, not SPOC (SPOC ate DepMap). See LEARNINGS.md.

Optional second held-out channel for extra depth: evolutionary coupling / co-conservation (independent of DepMap).

---

## 7. Calibration pre-check (Day-1, decides the headline branch)
Before building the full pipeline, quickly: fit isotonic regression mapping raw ipTM (or pDockQ2) → CORUM label
on the calibration split; plot the reliability diagram on the test split; compute ECE.
- If badly miscalibrated / overconfident → Branch A likely (artifacts exist). 
- If already near-diagonal → Branch B (lead with certified-core enrichment).
Either way keep the reliability diagram — it's demo frame 1.

---

## 8. Missing-member nomination
Pick a **target cancer complex** from the CM4AI map 275 assemblies — ideally one flagged as recurrently mutated
in cancer, with a plausible missing subunit. For each non-member candidate protein c:

```
score(c) = certified_confidence(c ↔ complex members)  ×  held_out_coessentiality(c, complex)
```

- `certified_confidence`: does c have edges to complex members that pass conformal certification (low conf p-value)?
- `held_out_coessentiality`: mean DepMap co-essentiality of c with complex members (the independent vote).
Nominate the top c with certified risk < q AND positive held-out co-essentiality. Output `nomination.json`
and have Claude Science write a literature-grounded rationale (PubMed/bioRxiv connector): what the complex does,
why c is plausible, what's known/unknown. This is the "discrete finding" the research track wants.

---

## 9. Nominee structural validation (now CORE — BioNeMo skill, was Stretch 1)
Because Boltz-2/OpenFold3 are callable BioNeMo skills in Claude Science (Modal-backed), predicting the
{nominated member + top complex partner} structure is a single reliable skill call — promote it from stretch
to **core**. Report ipTM/ligand_iptm + run the §4 physical-validity checks on the predicted interface. A
physically valid, confident predicted interface for a nominee *also* corroborated by held-out co-essentiality
is the wow moment. ~minutes for one complex; ~5–15 predictions well within the $30 Modal budget (DATA.md §5).

## 10. Stretch 2 — second-interactome robustness (Predictomes)
Repeat §2–§6 on Predictomes (~1.6M human pairs; predictomes.org/downloads). **Audit the raw pDockQ2 metric,
NOT SPOC** (SPOC ingests DepMap/co-expression → would break held-out purity). Show the audit conclusion
(overconfidence + certified-core enrichment) generalizes to a second, independently-built map. If the easy
Predictomes CSV lacks pDockQ2, use the small genome-maintenance subset
---

# Program extension (NEXT_DIRECTIONS, executed 2026-07-09)

## 11. Named audit
For interpretability, the 35 high-confidence edges dropped at q=0.10 are named and cross-referenced
against CORUM: 7 pair proteins co-occur in the same annotated complex (R2TP 8735, RNA exosome 788/789,
spliceosome pre-B 8370, NSL 7221, PRC1.5 8387/8388, Golgi SNARE 875/877). Upper-bound wording is
enforced repo-wide: a dropped edge is "not supported at level q", never "false". Output:
`reports/proposed_named_audit.md`, `results/tables/dropped_high_conf.csv`.

## 12. Mondrian (stratified) conformal
Conformal calibration is repeated within complex-size strata (small 3–5, medium 6–8, large 9),
matching each test edge to same-stratum decoys ("Mondrian" conformal, Vovk). Held-out FDR per stratum:
small 0.12 (wide CI), medium 0.08, large **0.65** (>3× q). Pooled calibration certifies 132 high-tertile
edges vs 78 under score-matched Mondrian (1.7× inflation) — the pooled audit's overconfidence is
*structured*, concentrated in large complexes. Output: `results/tables/mondrian_fdr.csv`,
`results/figures/overconfidence_structured.png`.

## 13. Covariate-shift robustness (weighted conformal)
The exchangeability assumption (VERIFY §method) is interrogated directly. The calibration decoys vs the
wild candidate pool differ significantly in score (KS D=0.134, p=1.4e-9; a benign same-distribution
sanity split gives KS 0.026, p=0.92). We test whether weighted conformal selection (WCS; Jin & Candès
2023, likelihood-ratio reweighting, ESS fraction 0.92) restores control on node-disjoint
(protein-disjoint) cal/test splits. Result: **both plain and WCS conformal exceed nominal q where they
certify**, and WCS mainly collapses power (certifies 0/2/16 edges at q=0.05/0.10/0.20 vs plain
14/22/27) without an FDR gain. The finite-sample guarantee does not survive the real shift — reported
as a limitation. The rigorous guarantee this work stands on remains the synthetic-null exchangeability
unit test (exact control); held-out CORUM/DepMap are corroboration. Output:
`results/figures/shift_diagnostic.png`, `wcs_vs_plain_fdr.png`, `data/processed/wcs_results.json`.

## 14. Baselines (head-to-head)
Conformal-BH is compared against three alternatives at the same nominal FDR=0.10: the map's own
benchmark-tuned fixed cutoff (Score≥0.394), the published high-confidence flag, and the dropped set.
Metric = held-out DepMap co-essentiality support (the referee), which none of the methods used.
Conformal-BH selects the fewest edges (132) yet has the highest support (41% co-essential, mean
−log10 GLS p 4.36) vs fixed cutoff (204, 32%, 3.14) and published (161, 37%, 3.84). SPOC is NOT
included as a baseline (it ingests DepMap → circular against this referee). Output:
`results/tables/baseline_headtohead.csv`, `results/figures/baseline_headtohead.png`.

## 15. Second-map audit → auditability boundary
The audit was attempted on two further public deposits. A distribution-free conformal audit requires
(i) per-pair raw confidence with ranking information intact and (ii) an exchangeable null (matched
decoys through the same pipeline). The Krogan host–pathogen deposit (Zenodo 15588019) ships
AF-Multimer ranked-0 models with NO tabulated scores; pDockQ is recoverable on CPU from the PDBs
(Bryant 2022: Cβ–Cβ interface contacts × interface pLDDT from B-factors — `src/emperor/mimicry.py`,
validated on 452 benchmark dimers), but the released set is positives-only and pDockQ saturates
(median 0.74, 93%≥0.5, no null). Predictomes exposes only the DepMap-contaminated SPOC composite plus
a contact count — no raw pDockQ, no null. **Finding: post-hoc auditability is a property of the data
RELEASE, not the map.** Only CM4AI (native decoys shipped) is auditable. Recommendation: interactome
releases must ship per-pair raw scores + a matched decoy null. Output:
`reports/second_map_auditability.md`, `results/figures/auditability_boundary.png`,
`data/mimicry/{apms_pdockq,auditability_boundary}.json`.

## 16. Cross-architecture pilot (conditional GO)
Whether an architecturally independent predictor corroborates the conformal verdict, and whether
cross-architecture *divergence* predicts correctness, are tested on 12 pilot dimers (6 certified + 6
dropped high-conf) folded on Boltz-2 (v2.2.1, A100-80GB, free ColabFold MSA). Two distinct tests,
reported separately: **(i) corroboration** — certified Boltz-2 ipTM mean 0.71 vs dropped 0.46
(Mann-Whitney p=0.021), AF-M vs Boltz-2 Spearman ρ=0.80 (significant → drops are not AF-M-specific
artifacts); **(ii) divergence-as-feature** (the gate) — |Score − Boltz ipTM| larger for dropped (0.13)
than certified (0.03), one-sided p=0.090 at n=12 (directional, NOT significant). Verdict: develop
cross-architecture divergence as a *candidate* nonconformity feature, do not integrate on this
evidence; a powered score-matched test on a map where certification and raw score are not separable is
required. Output: `results/figures/crossarch_pilot.png`, `data/structures/crossarch_verdict.json`,
`pilot_crossarch.csv`.

## 17. Nomination as FDR-controlled selection
The single §8 nomination generalizes. Every conformally certified edge (q=0.10) that bridges a CORUM
complex member to a non-member is a missing-member nomination; because the certified set is
BH-FDR-controlled at q, the collection of nominations inherits the same guarantee (≤q expected wrong,
upper bound). Selection is purely conformal; DepMap co-essentiality only RANKS within each set (purity
firewall). Result: 393 complexes with ≥1 certified nomination, 810 total; the §8 flagship
KANSL3→MLL1-WDR5 (5386) is present and top-ranked by co-essentiality within its set. Module:
`src/emperor/nominate_sets.py` (`make nominate_sets`), output `data/processed/nomination_sets.json`.

## 18. Motivating proposition
The prevalence-shift result (§ Results) is stated formally as a boxed two-part proposition
(`reports/proposition.md`): (a) conformal+BH is prevalence-invariant (Bates 2023 / BH-PRDS); (b) a
benchmark-tuned cutoff's realized FDR = (1−π)FPR / [(1−π)FPR + π·TPR] → 1 as prevalence π → 0. Boxed as
framing/motivation, explicitly NOT claimed as a theorem or contribution.
