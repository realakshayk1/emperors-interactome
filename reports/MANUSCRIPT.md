# When can you audit an AI interactome? A conditional distribution-free guarantee and the null it requires

*Draft manuscript — Emperor's Interactome. Target venue: Cell Systems / Bioinformatics / MLSB. Every numeric claim traces to a reproducible artifact in `data/processed/` (regenerable via `make audit-self`); citations to project artifacts are given as `[→ file.json]`.*

---

## One-sentence contribution (G-NOVELTY)

**We show that post-hoc, distribution-free FDR auditing of an AI-predicted interactome is only valid to the extent the map ships a calibration null that is exchangeable with its candidate edges — we quantify, on a real map, exactly how much exchangeability violation the guarantee tolerates before it breaks (Γ\*), and we argue that shipping a raw (untrained) confidence axis plus a native random-pair null should be a release standard for AI interactomes — a standard we test by transferring the audit to a second real map.**

This is a *conditional-guarantee* result, not a claim that we repaired unconditional control. We are explicit about where the marginal guarantee fails, and we convert that failure into a measurable robustness statement and a community standard.

---

## Abstract

Published AI interactomes label thousands of protein pairs "high-confidence," but those error rates are estimated on class-balanced benchmarks and inflate when true interactions are rare — the regime of a genome-scale screen. We re-audit the Krogan/Ideker AlphaFold-Multimer cell map (*Nature* 2025) with distribution-free conformal FDR control, calibrated against the map's own native random-pair decoys and refereed by held-out DepMap co-essentiality under a strict purity firewall (the referee never enters any score, label, or calibration set). At a target FDR of q=0.10, 35 of 161 (22%) of the map's high-confidence edges fail conformal FDR control, and certified edges are enriched for independent co-essentiality relative to dropped edges (41% vs 17%; permutation p=0.016). We then interrogate the guarantee itself. The audit's validity rests on exchangeability between the calibration decoys and the candidate edges; we measure a real calibration-to-candidate distribution shift of 0.60σ in nonconformity, driven predominantly by an unmodeled hub/degree-selection axis rather than by the confidence score, and we show the marginal FDR guarantee provably breaks at a shift of only δ*≈0.08σ. Under the measured shift the realized FDR is ≈0.30 at nominal q=0.10, corroborated by two independent routes (an identifying δ-curve and a protein-disjoint empirical estimate). Neither a hard-negative-enriched null nor degree-conditioned (Mondrian) calibration restores unconditional control. We therefore report a *conditional* result: the strongly-certified core of the map is robust to a bounded exchangeability violation up to Γ\*≈31 (Rosenbaum-style sensitivity), and it is independently enriched for physical interactions in wet-lab PPI data (77% vs a degree-matched null of 34%, odds ratio 6.2, permutation p<10⁻⁴). We confirm the audit transfers to a second, independent real map (genome-scale Predictomes), which also surfaces a sharpening of the standard: the auditable axis must be a *raw* structural-confidence readout, since a trained interaction classifier separates candidates from the null too perfectly to audit. Of current AI interactomes, CM4AI is — to our knowledge — the only one shipping such a raw auditable axis in bulk. We provide a semi-synthetic benchmark demonstrating that a fixed benchmark-tuned cutoff's realized FDR climbs to 0.84–0.92 as prevalence falls while conformal control holds at ≤q, and we release the full audit as a reproducible standard.

---

## 1. Introduction

AI structure prediction has moved from single structures to genome-scale interaction screens: AlphaFold-Multimer and successors are now used to nominate thousands of protein-protein interactions at once. A recurring feature of these releases is a "high-confidence" tier defined by a confidence-score cutoff, with an error rate estimated on a curated, roughly class-balanced benchmark (e.g. known complexes vs. random pairs at ~1:1).

Two problems follow. First, a threshold with a given precision on a balanced benchmark has a very different realized false-discovery rate in the true operating regime, where genuine interactions are a small minority of all scored pairs — the *prevalence shift* that inflates FDR precisely where confident-looking errors are most costly. Second, the natural fix — post-hoc, distribution-free error control via conformal p-values plus Benjamini-Hochberg — carries a guarantee that is only as good as its calibration null: the certified FDR is controlled if the calibration negatives are exchangeable with the candidate edges being tested.

This paper takes both problems seriously and, crucially, does not claim to have solved the second one by engineering a better null. Instead we (i) perform the audit and report what it removes; (ii) *measure* the exchangeability violation on a real map and show the marginal guarantee genuinely fails; (iii) convert that failure into a conditional robustness statement and independent physical validation; and (iv) argue that the ingredient that made this map auditable at all — a native, released null — should become a standard.

---

## 2. Data and firewall

**Primary map.** The Krogan/Ideker CM4AI multimodal cell map (Schaffer et al., *Nature* 2025, doi:10.1038/s41586-025-08878-3), AlphaFold-Multimer interface confidence (ipTM) as the confidence axis. The map ships **1,788 native random-pair decoys**, which we use as calibration negatives — the feature that makes the map auditable.

**Positive labels.** CORUM 5.3 human complexes, with a complex-disjoint calibration/test split so that no complex contributes to both calibration and evaluation.

**Referee (held out).** DepMap co-essentiality (Wainberg et al. 2021). The referee is subject to a strict **purity firewall**: it never enters any nonconformity score, label, calibration set, or hard-negative set. It is used only to (a) evaluate certified vs. dropped edges after the fact and (b) *rank* (never select) within an already-certified nomination set. The firewall is enforced by adversarial unit tests (§6).

---

## 3. The audit and what it removes

We define a nonconformity score s = (1 − ipTM) and compute conformal p-values for each of the 1,666 candidate edges against the calibration decoys, then apply Benjamini-Hochberg at q. At q=0.10, **132 candidate edges are certified**. Restricting attention to the map's own published high-confidence tier (161 edges), **35 of those 161 (22%) fail certification and are dropped** [→ audit_summary.json]. (The two counts have different denominators: 132/1,666 certified across all candidates; 35/161 dropped within the published high-confidence tier.) The dropped set includes edges to biologically plausible partners that nonetheless fail distribution-free control at the map's own decoy-implied null.

The independent referee agrees with the audit's direction: **certified edges are 41% co-essential vs 17% for dropped edges** (permutation p=0.016), while certified-vs-raw-high-confidence is not significantly different (p=0.51) — i.e. the audit's value is in *removing* the co-essentiality-poor tail, not in re-ranking the whole set [→ validation.json].

**Dependence.** Conformal p-values computed against a common calibration set are mutually dependent across test edges. Benjamini-Hochberg is nonetheless valid here because these p-values are provably PRDS across test points (Bates, Candès, Lei, Romano, Sesia, *Ann. Statist.* 51(1):149–178, 2023, doi:10.1214/22-AOS2244), and BH controls FDR under PRDS (Benjamini & Yekutieli 2001). For completeness we also report the arbitrary-dependence floor: under Benjamini-Yekutieli and e-BH, **zero edges certify at any q** — the harmonic penalty is unpayable at this calibration size (n_cal≈900, m=1,666 tests) [→ dependence_robustness.json]. We therefore rest the headline on PRDS, and report BY/e-BH as the honest worst case.

---

## 4. Interrogating the guarantee (the core result)

The audit above is only valid if the calibration decoys are exchangeable with the candidate edges. They are not, and we quantify the gap.

**4.1 The shift is real and is not about the score.** We measure a calibration-to-candidate shift of **0.60σ** in nonconformity units [→ gamma_seed.json]. Decomposing it across local covariates, the largest single divergence is not the confidence score (KS 0.129) but the **union-graph degree** of the edge's endpoints (KS 0.256) — a hub/selection axis. A density-ratio classifier separates calibration from candidate edges at AUC 0.64 on the full covariate vector but only 0.56 on the score alone: **62% of the divergence is invisible to one-dimensional score reweighting** [→ shift_attribution.json]. This is why an earlier weighted-conformal-shift attempt (reweighting on the score) could not close the gap: it was correcting the wrong axis.

**4.2 The marginal guarantee provably breaks — and by a lot.** Using an identifying experiment that toggles a synthetic shift δ and reads off realized FDR, control breaks (realized FDR exceeds q) at a shift of only **δ*≈0.08σ** [→ gamma_seed.json]. The measured shift of 0.60σ is roughly **7× past that breaking point**. Two independent routes put the realized FDR at nominal q=0.10 at **≈0.30**: the identifying-curve interpolation gives 0.295, and a protein-disjoint empirical held-out estimate gives 0.319 [→ sensitivity.json]. The unconditional guarantee does not hold on this map.

**4.3 Attempts to restore unconditional control fail honestly.** We tried the two natural fixes. (a) A *hard-negative* null enriched for high-degree decoys (pre-registered rule, degree in the top tertile of the candidate distribution) — its KS to the candidate distribution is if anything slightly higher (0.156 vs 0.129), an honest partial-negative, because degree and score are nearly independent [→ hard_negatives.json]. (b) *Degree-conditioned Mondrian* calibration — it collapses power (near-zero certifications). On protein-disjoint splits, neither fix moves the held-out FDR toward q [→ shift_control.json]. We report this as an attempted-and-insufficient negative, not a solved problem.

**4.4 What we *can* guarantee: conditional robustness (Path B headline).** The failure of the marginal guarantee does not make the certified edges worthless — it means their validity is conditional. We make that conditionality precise with a Rosenbaum-style sensitivity model: a Γ-bounded tilt on the certified set's conformal p-values. The specific certified edges have very small conformal p-values, so their aggregate worst-case FDR stays below q up to a **large tilt Γ\*≈31** (a genuine grid-crossing, not a censored bound), and at the measured shift's equivalent Γ≈1.8 their worst-case FDR is ≈0.006 [→ sensitivity.json].¹

To be explicit about what this does and does not claim — because the two numbers (realized FDR ≈0.30, robustness to Γ\*≈31) invite a "having-it-both-ways" reading — they are answers to two *different* questions about two *different* populations, and both are reported as caveats, not as a rescue:
- **(Q1) Is the map's high-confidence tier, as a whole, FDR-controlled under the real shift? No.** The realized FDR is ≈0.30 at nominal q=0.10 — a 3× miss. This is the population-level verdict and it is negative; the marginal guarantee is void.
- **(Q2) Are the *individual edges we certified* likely to be real? For the strongly-certified ones, yes, robustly.** Their conformal p-values are so small that even a tilt far larger than the measured one cannot lift their worst-case FDR above q. This is an edge-set property, not a rescue of the tier-level guarantee, and it says nothing about the edges Q1 counts as errors.

The two do not, and should not, coincide: a set can contain a controlled core and an uncontrolled tail simultaneously. The defensible claim is therefore narrow: *the map's overall marginal FDR guarantee is void under the measured shift; separately, the strongly-certified core is robust to bounded exchangeability violations far larger than the one actually present.* We do not claim the second statement repairs the first.

¹ Exact grid-crossing values from `sensitivity.json`: Γ\* = 29.0 (q=0.05), 31.0 (q=0.10), 31.0 (q=0.20); "≈31" throughout refers to the headline q=0.10 value. Measured-shift equivalent Γ = 1.82.

---

## 5. Independent validation of the certified set

**5.1 Orthogonal physical channel.** Every channel used so far is either the AF score or functional (co-essentiality). We add a *physical* channel independent of both: IntAct experimental protein-protein interactions. Certified edges are strongly enriched — **101 of 132 (77%) have independent physical evidence vs a degree-matched null of 34% (permutation p<10⁻⁴, matched-null odds ratio 6.2)** [→ experimental_ppi_referee.json]. (The naive odds ratio against the *unmatched* IntAct background rate of 13% is 21.9, but that figure is inflated by the hub confound; the degree-matched OR of 6.2 is the confound-controlled effect we report.) Dropped edges are **not** significantly depleted (40% vs 33%, p=0.85). This is an honest boundary: the audit flags edges that fail *statistical control*, not edges proven *physically false* — some dropped edges do have physical support. The physical evidence for the certified edges is not circular with the AP-MS screen that generated the candidates: over a 40-edge provenance subset it spans 5 source databases, 13 detection methods (co-IP, yeast two-hybrid, pull-down, size-exclusion, cross-linking — several predating AlphaFold), and 35 distinct publications.

**5.2 Method validation on known truth.** On a semi-synthetic benchmark (Beta laws fit to the real positive/negative score distributions), conformal+BH holds realized FDR ≤ q at every controlled prevalence π∈{0.30,0.10,0.05,0.02}, while a fixed benchmark-tuned cutoff's realized FDR inflates from 0.20 (π=0.30) to **0.84 at π=0.02 (q=0.10), and to 0.92 at q=0.20** [→ benchmark_synth.json]. This demonstrates the prevalence wedge on data where the truth is known, and confirms the conformal machinery is correctly calibrated in the exchangeable regime — isolating the real-data non-exchangeability of §4 as the sole source of the marginal-guarantee failure.

**5.3b Held-out member recovery (LOCO).** The nomination procedure's central promise is that certification recovers true complex members. We test this directly with a leave-one-complex-member-out design: certification depends only on ipTM versus the calibration decoys and never on complex membership, so we can ask, for each CORUM complex with ≥2 members in the candidate universe (73 complexes, 216 member trials), whether a member has a certified edge to another member of its own complex — a genuinely held-out recovery. The fair comparison group is *eligible impostors*: non-members that nonetheless have a candidate edge into the complex (276 impostor trials), which removes the trivial confound that a random protein simply has no edge to test. True members are recovered at **49.5% vs 23.2% for eligible impostors (odds ratio 3.3, within-complex label-permutation p<10⁻⁴)** [→ loco_validation.json]. This is the primary procedure-validation result; the temporal split originally planned for this role is infeasible on this data (§8, limitation 3), and the substitution is logged as a dated pre-registration deviation.

**5.3 A validated positive control (not a discovery headline).** The pipeline recovers KANSL3 as a member of the NSL/MLL1-WDR5 chromatin-modifying assembly (CORUM 5386/610/7221) — an annotation held out of the pipeline, where KANSL3→MSL carries the highest co-essentiality of any nomination. We present this as a labeled sanity check that the nomination procedure surfaces true biology, and we release a dated, pre-registered prospective shortlist of *novel* nominations (KANSL3 excluded), each with explicit confirm/refute criteria [→ prospective_nominations.md].

---

## 6. Why this map, and the auditability standard (paired contribution)

The single ingredient that made this audit possible is a **native random-pair null on a raw, untrained confidence axis**. Without a calibration negative that a user can access, there is no distribution-free p-value to compute and nothing to control — and, as we show below, the null must live on a *raw* score, not a trained interaction classifier.

**A second real map: the audit transfers (§6.1).** We tested generality on an independent, already-published AF-Multimer interactome: the genome-scale human **Predictomes** (Schmid & Walter, *Mol Cell* 2025; 20,196 proteins, 1,614,047 pairs — a different lab, pipeline, and protein universe, at no compute cost) [→ secondmap_audit.json]. Two findings result:

1. **The audit axis must be a raw structural readout, not a trained classifier.** Predictomes' headline score, SPOC, is a machine-learning classifier trained to separate interactors from random pairs. On that axis the high-confidence tier is *perfectly* separated from a random-pair null, so every conformal p-value collapses to the 1/(n+1) floor and the audit is degenerate (100% certification, nothing to discriminate). We therefore audit on `num_unique_contacts`, a raw interface-size readout with genuine candidate/null overlap — the analog of CM4AI's raw ipTM. On this axis the audit is well-behaved (SPOC≥0.9 tier, 12,767 edges; certified/dropped 11,642/1,125, 12,420/347, 12,765/2 at q=0.05/0.10/0.20). This explains *why* CM4AI's ipTM was auditable and sharpens the standard: ship a **raw** confidence axis, not only a curated score.

2. **A classifier-curated map has less for the audit to remove.** The DepMap co-essentiality referee is alive on Predictomes (co-essential fraction rises 0.08→0.43 across SPOC bins) but does *not* separate certified from dropped within the SPOC≥0.9 tier (certified 0.43 vs dropped 0.50, background 0.43). The reason is structural: Predictomes is *classifier-curated* — SPOC has already pre-filtered the tier, so it is uniformly co-essentiality-rich and the audit finds little junk to drop. CM4AI's high-confidence tier is an *uncurated* raw-ipTM threshold, which is exactly why the audit removes 22% of it (certified 0.41 vs dropped 0.17). The contrast is itself a result: **the audit's yield scales with how uncurated the map's confidence tier is.**

**Claim scope, corrected.** An earlier version of this work claimed CM4AI was "the only externally auditable AI interactome." The Predictomes audit shows that is too strong: the correct statement is that CM4AI is, to our knowledge, the only current map shipping a raw auditable confidence axis *in bulk* — Predictomes ships a full pair space, but its bulk-downloadable axis is a trained classifier (degenerate for auditing), while its raw AF-M ipTM resides only in ~1.5 TB of per-structure archives. An apples-to-apples ipTM-vs-ipTM cross-map referee comparison is left as future work needing that access. The auditability standard is therefore refined: release a **raw** confidence axis alongside any curated score.

We turn this into a positive recommendation: **AI-interactome releases should ship, or enable regeneration of, a calibration null exchangeable with their candidate edges — and should characterize the dominant selection axes (e.g. degree) along which that exchangeability is likely to fail.** Our shift-attribution procedure (§4.1) is exactly the diagnostic a release could run on itself. This is the resource/standard contribution, complementary to the conditional guarantee.

Supporting analyses (cross-architecture divergence; Mondrian per-stratum breakdown showing overconfidence concentrated in large complexes) are reported as supporting evidence, not headline.

---

## 6b. Related work and positioning

We position the contribution against four literatures; the novelty is their *intersection applied to a real AI interactome with a measured exchangeability budget*, not any single ingredient.

**Distribution-free error control.** Conformal prediction (Vovk et al. 2022; Angelopoulos & Bates 2023) provides finite-sample validity under exchangeability. We use the conformal *p-value* / outlier-detection formulation of Bates et al. (2023), which also supplies the result we lean on for multiplicity: conformal p-values against a common calibration set are PRDS across test points, so Benjamini-Hochberg (1995) controls FDR via the Benjamini-Yekutieli (2001) PRDS theorem. Our contribution here is not the machinery but its *audit* use — turning a published confidence tier into a rejection decision with an explicit q — and reporting the arbitrary-dependence e-BH / BY floor (Wang & Ramdas 2022; Benjamini-Yekutieli 2001) as the honest worst case, which for this calibration size is empty.

**The exchangeability assumption and shift.** That conformal validity degrades under distribution shift is known (Quiñonero-Candela et al. 2008; Barber et al. 2023, "conformal prediction beyond exchangeability"). Prior applied work typically either assumes exchangeability or corrects for a *known* covariate shift via reweighting. Our distinct move is diagnostic and quantitative: we *measure* the calibration-to-candidate shift on a real map (0.60σ), attribute it across covariates (finding degree, not the confidence score, dominant), and locate the exact breaking point δ*≈0.08σ — showing the shift is both real and 1-D-reweighting-invisible. To our knowledge this exchangeability *budget* has not been reported for an AI interactome.

**Sensitivity analysis.** The conditional guarantee borrows the logic of Rosenbaum-style sensitivity analysis from observational causal inference (Rosenbaum 2002): rather than assume no unmeasured confounding, bound how large a hidden bias Γ would have to be to overturn a conclusion. We transplant this to conformal auditing — bounding how large an exchangeability violation Γ the certified set tolerates before its worst-case FDR exceeds q (Γ*≈31). This framing is, to our knowledge, novel in the conformal-audit setting.

**AI interactomes and their evaluation.** Structure-based interactome prediction (AlphaFold-Multimer, Evans et al. 2021; Humphreys et al. 2021; Burke et al. 2023; the CM4AI map, Schaffer et al. 2025; the classifier-curated genome-scale Predictomes, Schmid & Walter 2025) is typically evaluated against curated positives with random-pair negatives at balanced ratios — the benchmark regime whose prevalence-dependence we make explicit (cf. Hart et al. 2014 on gold-standard error rates; Davis & Goadrich 2006 on prevalence and precision). Our referee (DepMap co-essentiality, Wainberg et al. 2021; CORUM positives, Tsitsiridis et al. 2022) is standard; what is new is (i) using the map's *own* native null for calibration, (ii) an orthogonal *physical* referee to break circularity, and (iii) the argument that shipping such a null should be a release standard — the "auditability boundary."

**In one sentence:** conformal outlier p-values (Bates 2023) + PRDS-justified BH + a Rosenbaum-Γ sensitivity bound (Rosenbaum 2002), applied to a real AI interactome with a *measured* shift budget and an orthogonal physical referee, is the combination this paper contributes.

---

## 7. Reproducibility and integrity

The full hardened audit runs via `make audit-self` (~1 min): dependence-robust FDR → shift attribution → hard negatives → shift-control gate → Γ seed → sensitivity bound → semi-synthetic benchmark → experimental-PPI referee. The purity firewall is enforced by 21 unit tests including 5 adversarial firewall checks (no DepMap token in any score source; certification invariant to injected/permuted DepMap columns) and a no-signal negative control. Every success threshold was pre-registered in `DECISIONS.md` before the run that tested it; deviations (e.g. the infeasible temporal split; the corrected sensitivity model) are logged with their corrections rather than hidden.

---

## 8. Limitations

1. **Cross-map referee replication.** The audit *method* transfers to a second real map (Predictomes, §6.1), but the CM4AI-strength referee separation (certified ≫ dropped co-essentiality) does not replicate there — because Predictomes' only bulk-accessible raw axis (interface contacts) is weaker than ipTM and its tier is classifier-curated. An apples-to-apples ipTM-vs-ipTM cross-map comparison needs the ~1.5 TB per-structure archives and is left as future work; the empirical referee numbers remain strongest on CM4AI.
2. **Conditional, not unconditional.** The headline guarantee is conditional on bounded exchangeability violation; the marginal FDR guarantee is void under the measured shift, and we say so.
3. **Temporal validation instrument infeasible here.** The map's 690-protein screen does not intersect recent CORUM additions, so a prospective time-split has zero positives [→ timesplit_validation.json]. The validation *construct* it was meant to serve — does the procedure recover held-out true members? — is instead measured by leave-one-complex-member-out recovery (§5.3b, OR 3.3, p<10⁻⁴), corroborated by the semi-synthetic benchmark and the orthogonal physical referee. What remains genuinely untested is *forward-in-time* generalization specifically (recovery of complexes curated after the map was built), which the LOCO design does not address.
4. **CORUM incompleteness.** Positive labels are incomplete, so certified-vs-dropped co-essentiality contrasts are lower bounds on the true separation.

---

## Figure plan

- **Fig 1 — The audit and the wedge** (`fig1_audit_wedge.png`, rendered). (A) FDR-vs-prevalence: benchmark cutoff climbing to 0.84 vs conformal held at q (semi-synthetic, q=0.10). (B) certified vs dropped co-essentiality (41% vs 17%; p=0.016). *Establishes problem + audit.*
- **Fig 2 — Interrogating the guarantee (core)** (`fig2_guarantee.png`, rendered as a single 3-panel composite). (A) Per-covariate KS: degree (0.26) is the dominant shift axis, not score (0.13). (B) Three-route realized-FDR reconciliation at nominal q=0.10 — BH marginal 0.10 vs identifying-curve 0.29 vs empirical node-disjoint 0.32, with control breaking at δ*≈0.08σ. (C) Rosenbaum-Γ worst-case FDR curve: certified core stays ≤ q up to Γ\*≈31, far past the measured Γ≈1.8. *This is the paper's spine.* (Standalone `fig2c_sensitivity_gamma.png` retained as the panel-C-only asset.)
- **Fig 3 — Independent validation** (`fig3_validation.png`, rendered). (A) IntAct enrichment: certified observed 77% vs degree-matched null 34% (matched-null OR 6.2, p<10⁻⁴); dropped not depleted (p=0.85). (B) Semi-synthetic conformal FDR ≤ q across π at all three q. 
- **Fig 4 — Second real map (Predictomes)** (`fig4_secondmap.png`, rendered). (A) The audit transfers on a raw axis: the SPOC≥0.9 high-confidence tier (12,767 edges) certified/dropped across q on the raw interface-contact axis — the trained SPOC classifier axis is itself degenerate (all p-values pin to the floor; not shown). (B) The DepMap referee is alive on this map (co-essential fraction climbs 0.08→0.43 across SPOC bins) but the SPOC≥0.9 tier is already uniformly co-essentiality-rich — classifier-curated, so little for the audit to drop, unlike CM4AI's uncurated ipTM tier. *Establishes generality + the raw-axis standard.*
- **Supp. Fig — KANSL3 positive control + Boltz-2 interface** (supporting; asset exists in repo, not re-rendered here).

**Figure assets (all rendered from `data/processed/*.json`; regenerable):**
- `fig1_audit_wedge.png` → Fig 1
- `fig2_guarantee.png` → Fig 2 (3-panel composite: KS attribution / three-route FDR / Γ* curve)
- `fig3_validation.png` → Fig 3
- `fig4_secondmap.png` → Fig 4 (second real map: raw-axis transfer + referee-alive gradient)
- `fig2c_sensitivity_gamma.png` → standalone Γ* panel (Fig 2C alone, e.g. for slides)
- `ws1_ws2_shift_control.png` (earlier composite): superseded by the cleaner Fig 1/Fig 2/Fig 3 assets; retained for provenance only.

All four manuscript figures are rendered and verified (geometric §9.1 + perceptual §9.2 vision check).

---

## References

- Angelopoulos, Bates (2023). *Conformal Prediction: A Gentle Introduction*. Foundations and Trends in ML 16(4). doi:10.1561/2200000101
- Barber, Candès, Ramdas, Tibshirani (2023). *Conformal prediction beyond exchangeability*. Annals of Statistics 51(2):816–845. doi:10.1214/23-AOS2276
- Bates, Candès, Lei, Romano, Sesia (2023). *Testing for outliers with conformal p-values*. Annals of Statistics 51(1):149–178. doi:10.1214/22-AOS2244
- Benjamini, Hochberg (1995). *Controlling the false discovery rate*. JRSS-B 57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Benjamini, Yekutieli (2001). *FDR control in multiple testing under dependency*. Annals of Statistics 29(4):1165–1188. doi:10.1214/aos/1013699998
- Burke et al. (2023). *Towards a structurally resolved human protein interaction network*. Nature Structural & Molecular Biology 30:216–225. doi:10.1038/s41594-022-00910-8
- Davis, Goadrich (2006). *The relationship between precision-recall and ROC curves*. ICML. doi:10.1145/1143844.1143874
- Evans et al. (2021). *Protein complex prediction with AlphaFold-Multimer*. bioRxiv. doi:10.1101/2021.10.04.463034
- Hart, Brown, Moffat (2014). *Measuring error rates in genomic perturbation screens*. Molecular Systems Biology 10:733. doi:10.15252/msb.20145216
- Humphreys et al. (2021). *Computed structures of core eukaryotic protein complexes*. Science 374(6573). doi:10.1126/science.abm4805
- Quiñonero-Candela et al. (2008). *Dataset Shift in Machine Learning*. MIT Press. doi:10.7551/mitpress/9780262170055.001.0001
- Rosenbaum (2002). *Observational Studies (2nd ed.)*. Springer Series in Statistics. doi:10.1007/978-1-4757-3692-2
- Schaffer et al. (2025). *A multimodal cell map of human cells (CM4AI)*. Nature. doi:10.1038/s41586-025-08878-3
- Schmid, Walter (2025). *Predictomes, a classifier-curated database of AlphaFold-modeled protein-protein interactions*. Molecular Cell 85(6):1216–1232.e5. doi:10.1016/j.molcel.2025.01.034
- Tsitsiridis et al. (2022). *CORUM: comprehensive resource of mammalian protein complexes – 2022*. Nucleic Acids Research 51(D1). doi:10.1093/nar/gkac1015
- Vovk, Gammerman, Shafer (2022). *Algorithmic Learning in a Random World (2nd ed.)*. Springer. doi:10.1007/978-3-031-06649-8
- Wainberg et al. (2021). *A genome-wide atlas of co-essential modules*. Nature Genetics 53:638–649. doi:10.1038/s41588-021-00840-z
- Wang, Ramdas (2022). *False discovery rate control with e-values*. JRSS-B 84(3):822–852. doi:10.1111/rssb.12489
