# When Can You Audit an AI Interactome? A Conditional Distribution-Free Guarantee and the Null It Requires

**Authors:** [Author Name(s)]¹ *(to be completed)*
**Affiliations:** ¹[Affiliation] *(to be completed)*
**Correspondence:** [email] *(to be completed)*

*Every numeric claim traces to a reproducible artifact in* `data/processed/` *(regenerable via* `make audit-self`*). Target venues: Cell Systems / Bioinformatics / MLSB.*

**Keywords:** conformal prediction; false discovery rate; distribution-free inference; protein–protein interactions; AlphaFold-Multimer; exchangeability; sensitivity analysis; AI interactome reliability

---

## Abstract

Published AI interactomes label thousands of protein pairs "high-confidence," but those error rates are estimated on class-balanced benchmarks and inflate when true interactions are rare — the regime of a genome-scale screen. We re-audit the CM4AI AlphaFold-Multimer cell map (Schaffer et al., *Nature* 2025; NIH Bridge2AI Cell Maps for AI project; contact PI Ideker, co-PI Krogan) with distribution-free conformal FDR control, calibrated against the map's own native random-pair decoys and refereed by held-out DepMap co-essentiality under a strict purity firewall (the referee never enters any score, label, or calibration set). At a target FDR of q=0.10, 35 of 161 (22%) of the map's high-confidence edges fail conformal FDR control, and certified edges are enriched for independent co-essentiality relative to dropped edges (41% vs 17%; permutation p=0.016). We then interrogate the guarantee itself. The audit's validity rests on exchangeability between the calibration decoys and the candidate edges; we measure a real calibration-to-candidate distribution shift of 0.60σ in nonconformity, driven predominantly by an unmodeled hub/degree-selection axis rather than by the confidence score, and we show the marginal FDR guarantee provably breaks at a shift of only δ*≈0.08σ. Under the measured shift the realized FDR is ≈0.30 at nominal q=0.10, corroborated by two independent routes (an identifying δ-curve and a protein-disjoint empirical estimate). Neither a hard-negative-enriched null nor degree-conditioned (Mondrian) calibration restores unconditional control. We therefore report a *conditional* result: the strongly-certified core of the map is robust to a bounded exchangeability violation up to Γ\*≈31 (Rosenbaum-style sensitivity), and it is independently enriched for physical interactions in wet-lab PPI data (77% vs a degree-matched null of 34%, odds ratio 6.2, permutation p<10⁻⁴). We confirm the audit transfers to a second, independent real map (genome-scale Predictomes), which also surfaces a sharpening of the standard: the auditable axis must be a *raw* structural-confidence readout, since a trained interaction classifier separates candidates from the null too perfectly to audit. Of current AI interactomes, CM4AI is — to our knowledge — the only one shipping such a raw auditable axis in bulk. We provide a semi-synthetic benchmark demonstrating that a fixed benchmark-tuned cutoff's realized FDR climbs to 0.84–0.92 as prevalence falls while conformal control holds at ≤q, and we release the full audit as a reproducible standard.

---

## 1. Introduction

AI structure prediction has moved from single structures to genome-scale interaction screens: AlphaFold-Multimer and successors are now used to nominate thousands of protein–protein interactions at once. A recurring feature of these releases is a "high-confidence" tier defined by a confidence-score cutoff, with an error rate estimated on a curated, roughly class-balanced benchmark (e.g. known complexes vs. random pairs at ~1:1).

Two problems follow. First, a threshold with a given precision on a balanced benchmark has a very different realized false-discovery rate in the true operating regime, where genuine interactions are a small minority of all scored pairs — the *prevalence shift* that inflates FDR precisely where confident-looking errors are most costly. Second, the natural fix — post-hoc, distribution-free error control via conformal p-values plus Benjamini–Hochberg — carries a guarantee that is only as good as its calibration null: the certified FDR is controlled if the calibration negatives are exchangeable with the candidate edges being tested.

This paper takes both problems seriously and, crucially, does not claim to have solved the second one by engineering a better null. Instead we (i) perform the audit and report what it removes; (ii) *measure* the exchangeability violation on a real map and show the marginal guarantee genuinely fails; (iii) convert that failure into a conditional robustness statement and independent physical validation; and (iv) argue that the ingredient that made this map auditable at all — a native, released null on a raw confidence axis — should become a release standard.

**Contributions.** Post-hoc, distribution-free FDR auditing of an AI-predicted interactome is only valid to the extent the map ships a calibration null that is exchangeable with its candidate edges. We (1) audit a published map and quantify what honest FDR control removes; (2) measure, on a real map, exactly how much exchangeability violation the guarantee tolerates before it breaks (Γ\*); (3) convert the failure of the marginal guarantee into a conditional robustness statement plus orthogonal physical and held-out-recovery validation; and (4) argue that shipping a raw (untrained) confidence axis plus a native random-pair null should be a release standard for AI interactomes — a standard we test by transferring the audit to a second real map. This is a *conditional-guarantee* result, not a claim that we repaired unconditional control: we are explicit about where the marginal guarantee fails, and convert that failure into a measurable robustness statement and a community standard.

---

## 2. Results

### 2.1 The audit removes 22% of the published high-confidence tier

We define a nonconformity score s = (1 − ipTM) and compute conformal p-values for each of the 1,666 candidate edges against the map's native random-pair calibration decoys, then apply Benjamini–Hochberg at q. At q=0.10, **132 candidate edges are certified**. Restricting attention to the map's own published high-confidence tier (161 edges), **35 of those 161 (22%) fail certification and are dropped** (the two counts have different denominators: 132/1,666 certified across all candidates; 35/161 dropped within the published high-confidence tier). The dropped set includes edges to biologically plausible partners that nonetheless fail distribution-free control at the map's own decoy-implied null.

The independent referee agrees with the audit's direction: certified edges are 41% co-essential versus 17% for dropped edges (permutation p=0.016), while certified-vs-raw-high-confidence is not significantly different (p=0.51) — i.e. the audit's value is in *removing* the co-essentiality-poor tail, not in re-ranking the whole set (Figure 1B).

Conformal p-values computed against a common calibration set are mutually dependent across test edges. Benjamini–Hochberg remains valid here because these p-values are provably PRDS across test points (Bates et al. 2023), and BH controls FDR under PRDS (Benjamini & Yekutieli 2001). For completeness we also report the arbitrary-dependence floor: under Benjamini–Yekutieli and e-BH, zero edges certify at any q — the harmonic penalty is unpayable at this calibration size (n_cal≈900, m=1,666 tests). We therefore rest the headline on PRDS and report BY/e-BH as the honest worst case (Methods 5.3).



![Figure 1](../results/figures/fig1_audit_wedge.png)

**Figure 1. The prevalence wedge and the audit.** (A) On a semi-synthetic benchmark with known truth, a fixed benchmark cutoff's realized FDR inflates to 0.84 as prevalence falls, while conformal + BH holds at the target q=0.10. (B) On the real map, certified edges are 41% co-essential versus 17% for dropped edges (held-out DepMap referee; permutation p=0.016); the raw high-confidence tier sits between (37%).

### 2.2 The exchangeability violation is real and is not about the score

The audit above is valid only if the calibration decoys are exchangeable with the candidate edges. They are not, and we quantify the gap. We measure a calibration-to-candidate shift of 0.60σ in nonconformity units. Decomposing it across local covariates, the largest single divergence is not the confidence score (KS 0.129) but the union-graph degree of the edge's endpoints (KS 0.256) — a hub/selection axis (Figure 2A). A density-ratio classifier separates calibration from candidate edges at AUC 0.64 on the full covariate vector but only 0.56 on the score alone: 62% of the divergence is invisible to one-dimensional score reweighting. This is why a weighted-conformal-shift correction on the score cannot close the gap — it corrects the wrong axis.

### 2.3 The marginal guarantee provably breaks, and attempts to restore it fail honestly

Using an identifying experiment that toggles a synthetic shift δ and reads off realized FDR, control breaks (realized FDR exceeds q) at a shift of only δ*≈0.08σ. The measured shift of 0.60σ is roughly 7× past that breaking point. Two independent routes put the realized FDR at nominal q=0.10 at ≈0.30: the identifying-curve interpolation gives 0.295, and a protein-disjoint empirical held-out estimate gives 0.319 (Figure 2B). The unconditional guarantee does not hold on this map.

We tried the two natural fixes. (a) A hard-negative null enriched for high-degree decoys (pre-registered rule: degree in the top tertile of the candidate distribution) has a KS to the candidate distribution that is if anything slightly higher (0.156 vs 0.129), an honest partial-negative, because degree and score are nearly independent. (b) Degree-conditioned (Mondrian) calibration collapses power (near-zero certifications). On protein-disjoint splits, neither fix moves the held-out FDR toward q. We report this as an attempted-and-insufficient negative, not a solved problem.

### 2.4 What we can guarantee: conditional robustness of the certified core

The failure of the marginal guarantee does not make the certified edges worthless — it means their validity is conditional. We make that conditionality precise with a Rosenbaum-style sensitivity model: a Γ-bounded tilt on the certified set's conformal p-values. The specific certified edges have very small conformal p-values, so their aggregate worst-case FDR stays below q up to a large tilt Γ\*≈31 (a genuine grid-crossing, not a censored bound), and at the measured shift's equivalent Γ≈1.8 their worst-case FDR is ≈0.006 (Figure 2C).

To be explicit about what this does and does not claim — because the two numbers (realized FDR ≈0.30, robustness to Γ\*≈31) invite a "having-it-both-ways" reading — they answer two *different* questions about two *different* populations, and both are reported as caveats, not as a rescue:

- **(Q1) Is the map's high-confidence tier, as a whole, FDR-controlled under the real shift? No.** The realized FDR is ≈0.30 at nominal q=0.10 — a 3× miss. This is the population-level verdict and it is negative; the marginal guarantee is void.
- **(Q2) Are the individual edges we certified likely to be real? For the strongly-certified ones, yes, robustly.** Their conformal p-values are so small that even a tilt far larger than the measured one cannot lift their worst-case FDR above q. This is an edge-set property, not a rescue of the tier-level guarantee.

The two do not, and should not, coincide: a set can contain a controlled core and an uncontrolled tail simultaneously. The defensible claim is therefore narrow — *the map's overall marginal FDR guarantee is void under the measured shift; separately, the strongly-certified core is robust to bounded exchangeability violations far larger than the one actually present* — and we do not claim the second statement repairs the first. (Exact grid-crossing values: Γ\* = 29.0, 31.0, 31.0 at q = 0.05, 0.10, 0.20; "≈31" refers to the headline q=0.10 value; measured-shift equivalent Γ = 1.82.)



![Figure 2](../results/figures/fig2_guarantee.png)

**Figure 2. Interrogating the guarantee.** (A) Per-covariate shift: union-graph degree (KS 0.256) is the dominant calibration-to-candidate axis, not the confidence score (KS 0.129). (B) Three routes agree the realized FDR under the measured shift is ≈0.30 at nominal q=0.10, with control breaking at δ*≈0.08σ. (C) Rosenbaum-Γ worst-case FDR: the certified core stays ≤ q up to Γ\*≈31, far past the measured Γ≈1.8.

### 2.5 Independent validation of the certified set

**Orthogonal physical channel.** Every channel used so far is either the AF score or functional (co-essentiality). We add a *physical* channel independent of both: IntAct experimental protein–protein interactions. Certified edges are strongly enriched — 101 of 132 (77%) have independent physical evidence versus a degree-matched null of 34% (permutation p<10⁻⁴, matched-null odds ratio 6.2; Figure 3A). The naive odds ratio against the unmatched IntAct background rate of 13% is 21.9, but that figure is inflated by the hub confound; the degree-matched OR of 6.2 is the confound-controlled effect we report. Dropped edges are not significantly depleted (40% vs 33%, p=0.85). This is an honest boundary: the audit flags edges that fail *statistical control*, not edges proven *physically false*. The physical evidence is not circular with the AP-MS screen that generated the candidates: over a 40-edge provenance subset it spans 5 source databases, 13 detection methods (co-IP, yeast two-hybrid, pull-down, size-exclusion, cross-linking — several predating AlphaFold), and 35 distinct publications.

**Method validation on known truth.** On a semi-synthetic benchmark (Beta laws fit to the real positive/negative score distributions), conformal+BH holds realized FDR ≤ q at every controlled prevalence π ∈ {0.30, 0.10, 0.05, 0.02}, while a fixed benchmark-tuned cutoff's realized FDR inflates from 0.20 (π=0.30) to 0.84 at π=0.02 (q=0.10), and to 0.92 at q=0.20 (Figure 3B). This demonstrates the prevalence wedge on data where the truth is known and confirms the conformal machinery is correctly calibrated in the exchangeable regime — isolating the real-data non-exchangeability of Section 2.2 as the sole source of the marginal-guarantee failure.

**Held-out member recovery.** The nomination procedure's central promise is that certification recovers true complex members. We test this directly with a leave-one-complex-member-out (LOCO) design: certification depends only on ipTM versus the calibration decoys and never on complex membership, so for each CORUM complex with ≥2 members in the candidate universe (73 complexes, 216 member trials) we ask whether a member has a certified edge to another member of its own complex — a genuinely held-out recovery. The fair comparison group is eligible impostors: non-members that nonetheless have a candidate edge into the complex (276 impostor trials), removing the trivial confound that a random protein simply has no edge to test. True members are recovered at 49.5% versus 23.2% for eligible impostors (odds ratio 3.3, within-complex label-permutation p<10⁻⁴). This is the primary procedure-validation result; the temporal split originally planned for this role is infeasible on this data (Section 4), and the substitution is logged as a dated pre-registration deviation.

**A validated positive control.** The pipeline recovers KANSL3 as a member of the NSL/MLL1-WDR5 chromatin-modifying assembly (CORUM 5386/610/7221) — an annotation held out of the pipeline, where KANSL3→MSL carries the highest co-essentiality of any nomination. We present this as a labeled sanity check that the nomination procedure surfaces true biology (not as a discovery headline), and release a dated, pre-registered prospective shortlist of *novel* nominations (KANSL3 excluded), each with explicit confirm/refute criteria.



![Figure 3](../results/figures/fig3_validation.png)

**Figure 3. Independent validation.** (A) Certified edges are enriched for independent IntAct physical evidence (77% vs a degree-matched null of 34%, matched-null OR 6.2, p<10⁻⁴); dropped edges are not significantly depleted (40% vs 33%, p=0.85). (B) Conformal control holds realized FDR ≤ q across prevalences on known-truth data.

### 2.6 Generality: a second real map and the raw-axis refinement

We tested generality on an independent, already-published AF-Multimer interactome: the genome-scale human Predictomes (Schmid & Walter 2025; 20,196 proteins, 1,614,047 pairs — a different lab, pipeline, and protein universe, at no compute cost). Two findings result.

First, **the audit axis must be a raw structural readout, not a trained classifier.** Predictomes' headline score, SPOC, is a machine-learning classifier trained to separate interactors from random pairs. On that axis the high-confidence tier is perfectly separated from a random-pair null, so every conformal p-value collapses to the 1/(n+1) floor and the audit is degenerate (100% certification, nothing to discriminate). We therefore audit on num_unique_contacts, a raw interface-size readout with genuine candidate/null overlap — the analog of CM4AI's raw ipTM. On this axis the audit is well-behaved (SPOC≥0.9 tier, 12,767 edges; certified/dropped 11,642/1,125, 12,420/347, 12,765/2 at q=0.05/0.10/0.20; Figure 4A). This explains *why* CM4AI's ipTM was auditable and sharpens the standard: ship a raw confidence axis, not only a curated score.

Second, **a classifier-curated map has less for the audit to remove.** The DepMap co-essentiality referee is alive on Predictomes (co-essential fraction rises 0.08→0.43 across SPOC bins; Figure 4B) but does not separate certified from dropped within the SPOC≥0.9 tier (certified 0.43 vs dropped 0.50, background 0.43). The reason is structural: Predictomes is classifier-curated — SPOC has already pre-filtered the tier, so it is uniformly co-essentiality-rich and the audit finds little junk to drop. CM4AI's high-confidence tier is an uncurated raw-ipTM threshold, which is exactly why the audit removes 22% of it (certified 0.41 vs dropped 0.17). The contrast is itself a result: the audit's yield scales with how uncurated the map's confidence tier is.

![Figure 4](../results/figures/fig4_secondmap.png)

**Figure 4. A second real map (Predictomes).** (A) The audit transfers on a raw interface-contact axis (SPOC≥0.9 tier, 12,767 edges; certified/dropped across q). (B) The DepMap referee is alive on this map (co-essential fraction climbs 0.08→0.43 across SPOC bins) but the tier is already uniformly co-essentiality-rich — classifier-curated, so little for the audit to drop, unlike CM4AI's uncurated ipTM tier.

---

## 3. Discussion

The practical claim is narrow and, we believe, robust: on this map, 22% of the published high-confidence tier does not survive honest distribution-free FDR control, and an orthogonal, held-out biological referee agrees with the direction of that removal. The methodological claim is more general and more important. Post-hoc auditing of an AI interactome inherits a guarantee only through its calibration null, and that null is almost never exchangeable with the candidate set in a genome-scale screen — here the dominant violation is a hub/degree-selection axis orthogonal to the confidence score, so the standard weighted-conformal shift correction, which operates on the score, cannot repair it.

Rather than hide this, we quantify it three ways — a breaking-point δ*, a realized-FDR estimate from two routes, and a Rosenbaum-Γ robustness bound on the certified core — and separate the tier-level verdict (marginal control is void) from the edge-level one (the strongly-certified core is robust). The generality experiment then converts the enabling ingredient into a prescription: an auditable map must ship a *raw* confidence axis plus a native random-pair null, because a trained classifier score is un-auditable by construction (it separates candidates from the null too perfectly), and the audit's yield scales with how uncurated the tier is.

## 4. Limitations

1. **The marginal FDR guarantee does not hold on the real map.** The realized FDR at nominal q=0.10 is ≈0.30. Our positive statements are conditional (Rosenbaum-Γ) or refer to the certified core, not to tier-level marginal control.
2. **The hard-negative and Mondrian fixes did not restore unconditional control.** We report them as attempted-and-insufficient, not as solutions.
3. **The temporal-holdout validation is infeasible on this data.** CORUM's served release is a single current snapshot (the release parameter is ignored server-side), and zero candidate pairs are newly-appearing post-cutoff complex pairs, so a literature-dated split yields no positives. The LOCO recovery test is the pre-registered substitute (a logged deviation); it is equal-or-harder and measures the same construct.
4. **The second-map referee does not separate certified from dropped within the SPOC≥0.9 tier.** This is a consequence of classifier curation, not of audit failure; we report it as an honest non-replication and explain the mechanism.
5. **The "only raw auditable axis in bulk" claim rests on our own survey**, hedged accordingly; it is a call for a standard, not a proof of uniqueness.
6. **Co-essentiality and physical-interaction referees are enrichment signals, not ground truth.** Some dropped edges have physical support; the audit flags statistical-control failures, not proven-false edges.

## 5. Methods

**5.1 Data.** Primary map: the CM4AI multimodal cell map (Schaffer et al., *Nature* 2025, doi:10.1038/s41586-025-08878-3; NIH Bridge2AI *Cell Maps for AI* project — contact PI Trey Ideker, UC San Diego; co-PI Nevan Krogan, Gladstone Institutes / UCSF, whose lab supplied the interaction mapping), using AlphaFold-Multimer interface confidence (ipTM) as the confidence axis. The map ships 1,788 native random-pair decoys, used here as calibration negatives — the feature that makes the map auditable. Positive labels: CORUM 5.3 human complexes, with a complex-disjoint calibration/test split so no complex contributes to both. Held-out referee: DepMap co-essentiality (Wainberg et al. 2021). Second map: human Predictomes (Schmid & Walter 2025).

**5.2 The purity firewall.** The held-out referee (DepMap co-essentiality) never enters any score, label, or calibration set used by the audit; it is consulted only after certification, for validation. This is enforced in code and covered by adversarial unit tests that fail if any referee-derived quantity reaches the scoring path.

**5.3 Conformal FDR control and its dependence justification.** Nonconformity s = (1 − ipTM); conformal p-values are computed for each candidate edge against the calibration decoys; Benjamini–Hochberg is applied at q ∈ {0.05, 0.10, 0.20}. These p-values are mutually dependent (shared calibration set); BH is valid because they are PRDS across test points (Bates et al. 2023) and BH controls FDR under PRDS (Benjamini & Yekutieli 2001). The arbitrary-dependence floor (Benjamini–Yekutieli, e-BH) is reported as the honest worst case: zero certifications at this calibration size.

**5.4 Shift attribution and breaking point.** The calibration-to-candidate shift is measured in nonconformity σ-units and decomposed per covariate by two-sample KS; a density-ratio classifier quantifies the fraction of divergence invisible to score-only reweighting. An identifying experiment toggles a synthetic shift δ and reads off realized FDR to locate the breaking point δ*.

**5.5 Sensitivity analysis.** A Rosenbaum-style Γ-bounded tilt is applied to the certified set's conformal p-values; the worst-case FDR is computed on a grid extended to Γ=50, with grid-crossings (Γ\*) reported and censoring flagged when no crossing occurs.

**5.6 Independent validation.** IntAct physical interactions provide an orthogonal channel; enrichment is tested against a degree-matched null by 10,000-permutation test, with a provenance/independence audit over source databases, detection methods, and publications. The semi-synthetic benchmark fits Beta laws to the real positive/negative score distributions and sweeps prevalence π over 100 seeds. LOCO recovery uses within-complex label permutation (10,000 permutations) against an eligible-impostor null. All randomness uses SEED=42.

**5.7 Reproducibility.** Every numeric claim traces to a result JSON in `data/processed/`, regenerable via `make audit-self` (the full hardened suite, ~1 min) from a bare clone; interim tables regenerate from raw data and the second-map score table auto-fetches from its public source. Unit, adversarial-firewall, and held-out-recovery tests (21 tests) accompany the pipeline.

## Data and Code Availability

All code, configuration, and result artifacts are in the project repository under an MIT license. Primary map: doi:10.1038/s41586-025-08878-3 (CM4AI). Second map: Predictomes (predictomes.org; doi:10.1016/j.molcel.2025.01.034). CORUM 5.3, DepMap co-essentiality (Wainberg et al. 2021), and IntAct are publicly available from their respective providers. Each numeric claim in this paper corresponds to a file in `data/processed/`; the pipeline regenerates all of them from a bare clone via `make audit-self`.

## Author Contributions

*(To be completed.)*

## Acknowledgments

*(To be completed — funding, compute, and data providers.)*

## References

- Angelopoulos, Bates (2023). *Conformal Prediction: A Gentle Introduction*. Foundations and Trends in Machine Learning 16(4). doi:10.1561/2200000101
- Barber, Candès, Ramdas, Tibshirani (2023). *Conformal prediction beyond exchangeability*. Annals of Statistics 51(2):816–845. doi:10.1214/23-AOS2276
- Bates, Candès, Lei, Romano, Sesia (2023). *Testing for outliers with conformal p-values*. Annals of Statistics 51(1):149–178. doi:10.1214/22-AOS2244
- Benjamini, Hochberg (1995). *Controlling the false discovery rate*. Journal of the Royal Statistical Society B 57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Benjamini, Yekutieli (2001). *The control of the false discovery rate in multiple testing under dependency*. Annals of Statistics 29(4):1165–1188. doi:10.1214/aos/1013699998
- Burke et al. (2023). *Towards a structurally resolved human protein interaction network*. Nature Structural & Molecular Biology 30:216–225. doi:10.1038/s41594-022-00910-8
- Davis, Goadrich (2006). *The relationship between precision-recall and ROC curves*. ICML. doi:10.1145/1143844.1143874
- Evans et al. (2021). *Protein complex prediction with AlphaFold-Multimer*. bioRxiv. doi:10.1101/2021.10.04.463034
- Hart, Brown, Moffat (2014). *Measuring error rates in genomic perturbation screens*. Molecular Systems Biology 10:733. doi:10.15252/msb.20145216
- Humphreys et al. (2021). *Computed structures of core eukaryotic protein complexes*. Science 374(6573). doi:10.1126/science.abm4805
- Quiñonero-Candela et al. (2008). *Dataset Shift in Machine Learning*. MIT Press. doi:10.7551/mitpress/9780262170055.001.0001
- Rosenbaum (2002). *Observational Studies* (2nd ed.). Springer Series in Statistics. doi:10.1007/978-1-4757-3692-2
- Schaffer et al. (2025). *Multimodal cell maps as a foundation for structural and functional genomics* (CM4AI). Nature. doi:10.1038/s41586-025-08878-3
- Schmid, Walter (2025). *Predictomes, a classifier-curated database of AlphaFold-modeled protein–protein interactions*. Molecular Cell 85(6):1216–1232.e5. doi:10.1016/j.molcel.2025.01.034
- Tsitsiridis et al. (2022). *CORUM: the comprehensive resource of mammalian protein complexes — 2022*. Nucleic Acids Research 51(D1). doi:10.1093/nar/gkac1015
- Vovk, Gammerman, Shafer (2022). *Algorithmic Learning in a Random World* (2nd ed.). Springer. doi:10.1007/978-3-031-06649-8
- Wainberg et al. (2021). *A genome-wide atlas of co-essential modules identifies functional gene interactions*. Nature Genetics 53:638–649. doi:10.1038/s41588-021-00840-z
- Wang, Ramdas (2022). *False discovery rate control with e-values*. Journal of the Royal Statistical Society B 84(3):822–852. doi:10.1111/rssb.12489

---

## Appendix A. Supplementary figures

- **Supp. Fig. S1 — KANSL3 positive control and Boltz-2 interface.** Recovery of the KANSL3/NSL-family positive control and its predicted interface (asset in repository; not re-rendered here).
- Additional diagnostic figures (`fdr_curve.png`, `prevalence_shift.png`, `reliability_raw.png`, `heldout_enrichment.png`, `nominee_structure.png`, `fig2c_sensitivity_gamma.png`) are available in `results/figures/`.

## Appendix B. Related work and positioning

*Conformal prediction and distribution-free inference.* We build on conformal p-values and outlier testing (Bates et al. 2023; Vovk et al. 2022; Angelopoulos & Bates 2023), FDR control (Benjamini & Hochberg 1995; Benjamini & Yekutieli 2001; Wang & Ramdas 2022), and conformal prediction beyond exchangeability (Barber et al. 2023), whose weighted-conformal correction we show is insufficient here because the dominant shift is off-score. *Dataset shift* framing follows Quiñonero-Candela et al. (2008). *Sensitivity analysis* follows Rosenbaum (2002). *AI interactomes and complexes* — AlphaFold-Multimer (Evans et al. 2021), computed eukaryotic complexes (Humphreys et al. 2021), structurally resolved human PPIs (Burke et al. 2023), CM4AI (Schaffer et al. 2025), and Predictomes (Schmid & Walter 2025). *Reference resources* — CORUM (Tsitsiridis et al. 2022), DepMap co-essentiality (Wainberg et al. 2021; Hart et al. 2014). To our knowledge this is the first post-hoc, distribution-free FDR audit of a published AI interactome that measures the exchangeability violation on the real map and converts it into a conditional guarantee plus a release standard.
