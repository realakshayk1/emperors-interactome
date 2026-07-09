# A distribution-free reliability audit of AI-predicted protein interactomes

*Draft skeleton — MLSB workshop / bioRxiv preprint. Target length ~5 pages. Numbers are pulled
verbatim from `data/processed/*.json` and the phase artifacts; citations keyed to `NEXT_DIRECTIONS §9`.*

---

## Abstract (~150 words)

AI structure predictors now generate protein interactomes at scale, published as sets of
"high-confidence" complexes gated on an interface confidence score (ipTM / pDockQ) with a
benchmark-estimated false-discovery rate. We show these self-reported error rates are
prevalence-fragile and that the underlying confidence axis is overconfident in a *structured* way.
We apply distribution-free **conformal selection** (conformal p-values + Benjamini–Hochberg) to the
CM4AI multimodal cell map (Schaffer et al., Nature 2025), calibrated against the map's own native
random-pair decoys and refereed by held-out DepMap co-essentiality — which never enters the score.
**35 of 161 (22%) published high-confidence edges fail FDR control at q = 0.10.** The dropped edges
are depleted of independent co-essentiality support (17% vs 41%, permutation p = 0.016); a second,
independent structure predictor (Boltz-2) corroborates the verdict. We further show that whether such
an audit is even *possible* is a property of the data release, not the map: of three public deposits,
only the one shipping a native null is auditable.

## 1. Introduction / motivation

- The field ships AF-Multimer/pDockQ interactomes with author-set cutoffs and benchmark-estimated
  FDRs (importance hook: 2026 AFDB complex expansion, millions of confidence-scored complexes).
- Benchmark FDR is a prevalence-fragile lookup (boxed **Proposition**, §`reports/proposition.md`):
  positive half — conformal+BH is prevalence-invariant; negative half — a benchmark-tuned cutoff's
  realized FDR → 1 as interactions become rare. Elementary (Bates 2023 + BH/PRDS); stated as
  motivation, not a contribution.
- **Contribution (reframed, guardrail-compliant):** the first distribution-free, exchangeability-
  *audited* conformal reliability audit of published AF-M interactomes — FDR-controlled *selection*
  applied to AF-M PPI scores, with the exchangeability assumption interrogated head-on. We do not
  claim "certified rejection" as a new capability (it is the complement of selection), nor a new
  network-dependence theorem (Marandon 2023 / Du 2025 own that; we cite/apply).

## 2. Methods

- **Data.** CM4AI map (1,666 candidate pairs + 1,788 native random decoys; axis = ipTM). CORUM 5.3
  for calibration labels + ID crosswalk. DepMap co-essentiality (GLS matrix) as held-out referee.
  **Purity firewall:** DepMap/co-essentiality never enters the nonconformity score or calibration
  labels; SPOC is never audited (it ingests DepMap → circular).
- **Conformal selection.** Nonconformity `s = (1 − score) + w·phys_penalty`; conformal p-values vs the
  decoy null; BH at q ∈ {0.05, 0.10, 0.20}. Complex-disjoint calibration/test split. Held-out FDR
  reported as an UPPER bound (CORUM incompleteness), both pooled and mean-split estimators.
- **Structure/physics feature.** Boltz-2 ipTM as an independent-architecture signal + interface
  physical-validity (one feature in a multi-signal score; AdaDetect legitimises a learned score —
  Marandon et al. 2024). Not a headline (SPOC owns the physics-classifier framing).

## 3. Results

### 3.1 The confidence axis is miscalibrated → Branch A
AUROC 0.70; ECE 0.176 → 0.016 after isotonic recalibration. (Fig: `reliability_raw.png`.)

### 3.2 22% of high-confidence edges fail honest FDR control
Certified-by-q {0.05: 78, 0.10: 132, 0.20: 177}; high-conf dropped {0.05: 83, 0.10: 35, 0.20: 12}.
**Headline: 35/161 (22%) fail at q = 0.10.** Named casualties in same-CORUM-complex pairs (R2TP,
RNA exosome, spliceosome pre-B, NSL, PRC1.5, Golgi SNARE). (Figs: `fdr_curve.png`,
`dropped_high_conf.csv`; §`proposed_named_audit.md`. Upper-bound wording throughout.)

### 3.3 The prevalence wedge
Benchmark-tuned cutoff realized FDR climbs 0.08 → 0.90 as prevalence 0.5 → 0.01; conformal+BH stays
near q, certifying nothing rather than exceeding it at low prevalence. (Fig: `prevalence_shift.png`.)

### 3.4 Overconfidence is structured (Mondrian)
Pooled conformal holds in medium complexes (held-out FDR 0.08) but is marginally above target in small
(0.12, wide CI) and **breaks in large 9-member complexes (0.65, >3× target)**. (Fig:
`overconfidence_structured.png`.)

### 3.5 Held-out referee: dropped edges are depleted of support
Certified 41% co-essential vs dropped 17% (permutation p = 0.016); certified-vs-dropped is the
discriminating contrast (certified-vs-raw is confounded by subset relation, p = 0.51). (Fig:
`heldout_enrichment.png`.)

### 3.6 Baselines
At the same nominal FDR = 0.10, conformal-BH selects fewest edges (132) yet highest held-out support
(41%, mean −log10 p 4.36) vs fixed cutoff (204, 32%) and published flag (161, 37%). (Fig:
`baseline_headtohead.png`.)

### 3.7 Covariate-shift robustness (honest negative)
The calibration→wild shift is real (KS 0.134, p = 1.4e-9). On node-disjoint splits **both plain and
weighted (WCS) conformal exceed nominal q where they certify**; WCS mainly collapses power without an
FDR gain. The guarantee does **not** survive the real shift — reported as a limitation, not a rescue.
(Figs: `shift_diagnostic.png`, `wcs_vs_plain_fdr.png`.) *This is the paper's core defensibility: the
shift is diagnosed and quantified, not hidden.*

### 3.8 Cross-architecture pilot (Pillar B, GO with caveat)
Boltz-2 corroborates the conformal verdict: certified ipTM 0.71 vs dropped 0.46 (Mann–Whitney
p = 0.021); AF-M vs Boltz-2 Spearman ρ = 0.80. Underpowered (n = 12), and certification is
near-separable from AF-M score on this map, so this shows *agreement* (drops aren't AF-M-specific
artifacts), not yet independent signal beyond score. Funded scale-up is the path to a calibrated
feature. (Fig: `crossarch_pilot.png`.)

### 3.9 Nomination as FDR-controlled selection
393 complexes receive ≥1 certified missing-member nomination (810 total, FDR ≤ q = 0.10 upper bound).
Flagship worked example: **KANSL3 → MLL1-WDR5 (CORUM 5386)**, certified risk 0.0066, top-ranked by
DepMap co-essentiality (3.07), independently confirmed as an annotated NSL-complex member (CORUM
7221) absent from the target entry, and structurally corroborated by Boltz-2 (localized interface
KANSL3[279–448]–KANSL1[624–690]). (§`RATIONALE.md`, `nomination_sets.json`.)

### 3.10 Auditability is a property of the release
Of three public AF-M deposits, only CM4AI (ships native decoys) is auditable. The Krogan host–pathogen
map is positives-only with a saturated pDockQ axis (median 0.74, 93% ≥ 0.5); Predictomes exposes only
the DepMap-contaminated SPOC composite. **Recommendation: interactome releases must ship their nulls.**
(Fig: `auditability_boundary.png`; §`second_map_auditability.md`.)

## 4. Limitations

- Held-out FDR is an upper bound (CORUM incomplete); the rigorous guarantee is the synthetic-null unit
  test, held-out CORUM/DepMap are corroboration.
- The covariate-shift guarantee does not survive the real shift (§3.7) — the contribution there is the
  diagnosis.
- Cross-architecture signal is a GO-to-develop, not a calibrated feature (§3.8).
- Single fully-auditable map; generalization is carried by the auditability-boundary analysis, not a
  second full audit.

## 5. Conclusion

A distribution-free audit turns "high-confidence" from an author-set label into a statement with a
stated, prevalence-robust error rate — and reveals that 22% of one flagship map's high-confidence
edges do not survive it. The broader lesson for the field is a release standard: ship the null, or the
reported FDR cannot be independently checked.

## References
See `NEXT_DIRECTIONS.md §9` (Jin & Candès 2022; Bates 2023; Marandon 2023; Blanchard 2024; Du 2025;
Marandon-AdaDetect 2024; Schaffer et al. 2025; SPOC/Predictomes 2025).
