# DECISIONS — append-only decision log

ADR-style. Newest at bottom. Mark superseded entries rather than deleting. Captures *why* the project
is shaped this way so a fresh agent doesn't re-litigate settled choices.

## 2026-07-08 — Project = falsification audit, not a trust-layer tool
- **Status:** accepted
- **Context:** The obvious idea (recalibrate AF-Multimer confidence) is crowded and reads as tool-building.
- **Decision:** Reframe as an *audit* of a published interactome — certify what survives honest FDR, refereed by held-out biology.
- **Rationale:** Contrarian + falsifiable is more novel and more memorable; wins on either result branch.
- **Alternatives rejected:** Boltz-2 affinity trust layer (methods-not-biology, GPU-bound, no Gladstone data); Krogan "obvious" PPI trust layer (the idea anyone would generate).

## 2026-07-08 — Primary interactome = Krogan/Ideker Nature 2025
- **Status:** accepted
- **Context:** Need a published, high-profile structural interactome that is DepMap-independent and cancer-relevant.
- **Decision:** Use the Nature 2025 "Multimodal cell maps" (U2OS) Suppl. Table 5 as primary.
- **Rationale:** Small/clean (1,666 pairs, 111 complexes), cancer context matches DepMap's cancer cell lines, and it's an AP-MS+AF map so DepMap stays genuinely held-out. Gladstone/Krogan data → Award-eligible.
- **Alternatives rejected:** Predictomes as primary (SPOC ingests DepMap → held-out purity broken); demoted to Stretch 2 with raw-pDockQ2 auditing.

## 2026-07-08 — Confidence axis = raw AF-M metrics + physical validity (never SPOC on primary)
- **Status:** accepted
- **Rationale:** Held-out purity: any score feeding the audit must exclude DepMap. SPOC ate DepMap/co-expression, so auditing SPOC and validating on DepMap would be circular.
- **Consequence:** On Predictomes (Stretch 2) audit raw pDockQ2, not SPOC.

## 2026-07-08 — FDR machinery = conformal p-values + Benjamini-Hochberg
- **Status:** accepted
- **Rationale:** Distribution-free, tiny code, provable FDR, and the exact wedge vs SPOC's prevalence-fragile benchmark-FDR.
- **Alternatives rejected:** Bespoke Bayesian model (overkill, not distribution-free), plain Platt calibration only (no selection guarantee).
- **Guardrail:** Do NOT claim to invent conformal FDR — cite Marandon 2023 / Blanchard 2024 / Jin-Candès 2023 as machinery; our contribution is the first structural-PPI-confidence instantiation + held-out audit.

## 2026-07-08 — Held-out referee = DepMap co-essentiality
- **Status:** accepted
- **Rationale:** Orthogonal to structure, cancer-cell-line-derived (matches U2OS), precomputed matrix available with no login. Purity firewall enforced in AGENTS.md/LEARNINGS.md.
- **Optional second channel:** evolutionary coupling (independent of DepMap) for extra depth if time allows.

## 2026-07-08 — Scope = core + 2 stretches; two headline branches both count as wins
- **Status:** accepted
- **Decision:** Core (audit + validation + 1 nomination) + Stretch 1 (Boltz/AF3 structural validation of nominee) + Stretch 2 (Predictomes robustness). Headline is Branch A (artifacts found) or Branch B (certified-core enrichment) — decided by the Day-1 pre-check; both are strong.
- **Rationale:** De-risks the "what if the map is well-calibrated" failure mode. Cut order defined in TASKS.md.

## 2026-07-08 — Compute = local/Claude-Science except Boltz (16 GB RAM)
- **Status:** accepted
- **Rationale:** DepMap matrix sliced to interactome genes → fits in RAM; all audit/validation is CPU. Only Boltz-2/AF3 (Stretch 1) needs GPU → Colab/Modal.

## 2026-07-08 — BioNeMo skills + $30 Modal: structure prediction promoted to CORE
- **Status:** accepted (v3); **EXECUTED 2026-07-09** — Modal connected, Boltz-2 run on A100-80GB.
- **Decision:** Physical-validity of the nominee's interface (Boltz-2 / OpenFold3 / AF2-Multimer via
  BioNeMo skills, Modal-backed, ~$30 credit ≈ 12 A100-hrs) moves from optional Stretch 1 to core.
- **Realized:** Boltz-2 v2.2.1 (image `im-HUmTG4YGPRXCeHPD6lJUXV`, weights hydrated on the user's
  Modal account) predicted the full-length KANSL3–KANSL1 dimer (2,009 aa) on an A100-80GB. Result:
  global ipTM 0.47, with a localized confidently-placed interface KANSL3[279–448]–KANSL1[624–690]
  (90 inter-chain contacts, interface pLDDT 0.70/0.59, best-contact PAE 4.8 Å). This is an
  **independent** structure model (architecturally distinct from the audited AlphaFold-Multimer),
  so it corroborates the certified edge rather than restating it. `src/emperor/structure.py` documents
  the exact remote invocation and defines `phys_penalty(iptm)=max(0, 0.5−iptm)` (=0.031 here). The
  shipped audit still uses the single-metric ipTM path (`phys_penalty`=0 for the map-wide FDR run);
  the structure term is applied to the nominee as explicit corroboration, saved in
  `data/processed/nomination.json → physical_validity` and `data/structures/`.
- **GPU-tier note:** the user's Modal plan gates every GPU ≥40 GB behind a payment method; the 24 GB
  tiers (A10G/L4/T4) OOM on this 2,009-residue dimer. A payment method was added to unlock A100-80GB.
  The ColabFold MSA is cached in the `claude-science-boltz-cache` volume, so re-runs are inference-only
  (~5 min). Cut-order (TASKS.md) still keeps the nominee valid without a predicted structure.

## 2026-07-09 — Confidence axis = ipTM (verified), single-metric path
- **Status:** accepted. Krogan Table 5 exposes ipTM (`iptm_0..4`), pTM (`ptm_0..4`), combined `Score`
  (≈ mean AF `0.8·ipTM+0.2·pTM`, verified: mean|err|=0.012), the paper's own `FDR`, and boolean flags.
  **No pDockQ2, no tabulated interface-PAE** — so the nonconformity score is single-metric on `Score`.
- **Consequence:** METHODS' multi-metric nonconformity reduces to `(1-score)+w_phys·phys_penalty`.

## 2026-07-09 — Calibration negatives = 1,788 native random pairs
- **Status:** accepted. Table 5 ships 1,788 "random" decoy pairs scored by the SAME AF-M pipeline,
  DepMap-independent. Used as conformal calibration negatives (cleaner exchangeability than ad-hoc
  decoys). Positives = CORUM same-complex pairs selected on membership alone (191 total), removed from
  negatives. Cal/test split is complex-DISJOINT for positives.

## 2026-07-09 — Branch A (artifacts likely) LOCKED
- **Status:** accepted (Day-1 pre-check). Raw `Score` is miscalibrated vs held-out CORUM
  (ECE 0.176 raw → 0.016 isotonic; AUROC 0.701). Rule `ece_raw > 2·ece_isotonic` → Branch A.
- **Headline:** at q=0.10, 35/161 (22%) paper high-confidence edges fail conformal FDR.

## 2026-07-09 — Held-out FDR reported as an UPPER bound, both estimators
- **Status:** accepted. Monte-Carlo over 200 complex-disjoint splits. Report pooled AND mean-split
  FDR side by side (q=0.10: 0.157/0.134; q=0.20: 0.215/0.175) — both run modestly above nominal at
  q=0.10, near/at nominal at q=0.20. This is an UPPER bound (CORUM incompleteness: a "random" decoy
  certified as an edge may be a real-but-unannotated interaction). The RIGOROUS guarantee is the
  synthetic-null unit test (exact control); held-out CORUM corroborates; DepMap is the referee.

## 2026-07-09 — Nomination target = MLL1-WDR5 (CORUM 5386), nominee KANSL3
- **Status:** accepted. Target chosen (config `TARGET_COMPLEX_ID`), NOT auto-picked: a leukemia-
  defining H3K4-methyltransferase/chromatin assembly anchored by KMT2A(MLL1) and WDR5, with many
  certified edges to non-members and good DepMap coverage.
- **Nominee KANSL3:** certified risk 0.0066 < q; strongly co-essential with the NSL submodule
  (KANSL1 p=2e-20, MCRS1 p=1e-6, PHF20 p=6e-3). Independent positive control held out of the pipeline:
  KANSL3 is an annotated NSL-complex member (CORUM 7221) absent from the target entry 5386.
- **Alternatives ranked (nom_score):** TAF5 (2.07), TAF8 (1.27) — TAF paralogs, less specific than the
  NSL-scaffold recovery. Honest scope: KANSL3 is essentiality-relevant, not a Cancer Gene Census driver.

## 2026-07-10 — PLAN_V3 execution (WS0 + WS1 diagnosis)

### Pre-registered gates (written BEFORE running; PLAN_V3 loop discipline)
- **T0.4 (BY/e-BH):** gate = produce BH vs BY vs e-BH certified/dropped side-by-side; assert BY⊆BH per q. Passes on producing the comparison.
- **T1.1 (shift attribution):** gate = quantify the share of cal↔wild divergence explained by `score` alone vs the full local covariate vector (density-ratio AUC + NN coverage). Passes on producing the attribution, whatever it shows.
- **T1.2 hard-negative rule (when built):** construction rule + params fixed here BEFORE KS-to-wild is measured; KS-to-wild is a reported diagnostic, never an optimization target (V3 markup #1 anti-circularity).

### T0.4 — dependence-robust FDR (DONE)
BH certifies 132/35 (cert/dropped) at q=0.10; **BY and e-BH certify 0 at every q**. Diagnosis (not a bug): the arbitrary-dependence harmonic penalty H_m≈8 over m=1666 needs leading conformal p-values ~1e-4, but the calibration floor is 1/(n+1)=0.0011 and only 33 edges beat the entire null. **Consequence — headline justification corrected:** the 132/35 result is defended by PRDS (conformal p-values against a common calibration set are PRDS across test points, Bates et al. 2023 Ann. Statist. 51(1):149-178 doi:10.1214/22-AOS2244, verified 2026-07-10; BH valid under PRDS via Benjamini-Yekutieli 2001), NOT by independence. BY/e-BH are reported as the (unpayable-at-this-size) arbitrary-dependence floor. Output: data/processed/dependence_robustness.json; tests/test_dependence.py (4 tests incl. exact e-BH synthetic control).

### T1.1 — shift attribution (DONE) → decides WS1 path
Full-covariate density-ratio classifier AUC 0.643 vs score-only 0.555 ⇒ **62% of the cal↔wild divergence is INVISIBLE to the 1-D `score` reweighting WCS used** — confirming the plan's hypothesis that the WCS failure is a null-completeness problem, not a reweighting bug. Largest single shift axis = union-graph **degree** (KS 0.256, vs score 0.129): candidates are hub-enriched by AP-MS+AF pre-screening, a selection covariate the score-only ratio never sees. (Self-audit: degree recomputed on the union graph, not the candidate graph, so it is not definitionally tied to pool membership — this dropped the invisible share from an inflated 82% to an honest 62%.) NN coverage of candidates by decoys = 80% (decoys mostly span the region). Output: data/processed/shift_attribution.json.

### Γ* seed (V3 markup #5) — Path B is viable and cheap
Mapping the measured real shift (cal-neg vs wild mean nonconformity gap = 0.60σ) onto the §19 identifying-experiment curve gives an **estimated realized FDR ≈ 0.29 at nominal q=0.10** (0.18 at q=0.05, 0.45 at q=0.20). This independently corroborates the empirically-measured node-disjoint held-out FDR of 0.319 (WCS experiment) — two independent routes agree at ~0.30. δ* (shift at which control breaks) ≈ 0.08σ, so the measured 0.60σ shift is ~7× past the breaking point. **This makes Path B (conditional/Γ* guarantee) a co-equal headline candidate: the honest statement is "realized FDR ≈ 0.29 under the measured 0.60σ exchangeability violation," a quantified correction rather than a broken guarantee.** Output: data/processed/gamma_seed.json.

### T1.2 — hard-negative rule PRE-REGISTERED (before KS-to-wild measured; V3 markup #1)
**Rule (fixed now, not tuned to KS):** the WCS/T1.1 diagnosis is that the dominant *unmodeled* shift axis is union-graph **degree** (candidates are hub-enriched by AP-MS+AF pre-screening; KS 0.256). Therefore hard negatives = **non-CORUM-co-complex pairs whose union-graph degree is in the top tertile of the candidate degree distribution** — i.e. high-connectivity pairs that look "selected" like the wild candidates but are not annotated interactors. Source pool = the existing decoy set augmented with degree-matched non-co-complex pairs drawn from the interactome's own protein universe. Parameters fixed: tertile cutoff = candidate-degree 67th percentile; firewalled from DepMap and from CORUM positives.
**Anti-circularity commitments:** (i) KS(hard-null, wild) is REPORTED as a diagnostic, never minimized by adjusting the rule; (ii) required output = contamination sensitivity propagating ε∈{1,5,10}% true-interactors-in-null to the certified count; the direction of bias (anti-conservative) is reported regardless of ε's true value.
**Gate:** hard-null recalibration moves node-disjoint held-out FDR toward q relative to the plain-conformal baseline (0.29/0.32/0.37 at q=0.05/0.10/0.20). Pass = movement toward q; a null result is reported honestly, not tuned away.

### WS1 VERDICT — Path A fails, Path B leads (2026-07-10)
T1.2/T1.3 evaluated on protein-disjoint splits within one consistent harness (shift_control.json). Within-harness plain-vs-fix comparison (the valid one; the stored WCS baseline used a different split scheme so cross-comparison is apples-to-oranges and is NOT claimed): at q=0.10 pooled held-out FDR plain 0.176 → hard-negative 0.175 (no material movement) → degree-Mondrian collapses power (n≈0). All methods certify very little on protein-disjoint splits (mean certified per split ranges 0–15 across methods×q; three of nine method×q cells are 0). **Neither a better null (T1.2) nor degree-conditioning (T1.3) restores unconditional control.** Per the pre-registered gate and PLAN_V3 WS1 kill criterion, route to **Path B (conditional Γ* guarantee)** as the headline — NOT Path C, because the Γ* seed shows a *quantified* conditional statement is available (realized FDR ≈0.29 at q=0.10 under the measured 0.60σ shift, corroborated by the empirical node-disjoint 0.319). Honest negative preserved: the hard-negative fix is reported as attempted-and-insufficient, not hidden. Next: T1.4 full Rosenbaum-Γ sensitivity bound = the lead deliverable.

### T1.2 pre-registration DEVIATION (flagged, dated 2026-07-10, before final verdict framing)
The pre-registered rule said "source pool = existing decoys AUGMENTED with degree-matched non-co-complex pairs drawn from the protein universe." **This is infeasible as written:** newly-drawn pairs were never folded by AlphaFold-Multimer, so they have NO `score`/ipTM and cannot enter a score-based conformal null without GPU folding (the deferred WS2.1 path). Feasible substitution actually implemented: **partition the existing decoy set by the pre-registered degree cutoff and UPWEIGHT the hard (high-degree) region** as the shift-matched null. This preserves the rule's intent (make the null span the hub-selected region the candidates occupy) using only scored pairs. The degree cutoff (candidate 67th pct) and the anti-circularity commitments are UNCHANGED. Drawing+folding new hard negatives is logged as the GPU extension (WS2.1-adjacent). This deviation does not affect the WS1 verdict (Path A fails under both the partition and would require folding to test the true-augmentation variant).

### T1.4 confirmation of the WS1 verdict (post-hoc, sequencing note)
SEQUENCING CAVEAT: the WS1 Path-B verdict above was first written citing the preliminary Γ* *seed* (identifying-curve interpolation), before the formal T1.4 sensitivity module was built — a deviation from strict "confirm-then-conclude" ordering. T1.4 has now been implemented (sensitivity.py) and run, and it CONFIRMS the routing rather than revising it: (A) population realized FDR under the measured 0.60σ shift = 0.295 (identifying curve) / 0.319 (empirical node-disjoint) at q=0.10 — the marginal guarantee is void (δ*≈0.08σ, shift ~7× past breaking); (B) the certified SET's aggregate worst-case FDR under a Γ-bounded tilt stays below q up to Γ*≈30 (grid-crossing, not censored), so the strongly-certified core is robust even though the marginal guarantee is not. Path B stands, now with a two-part quantified statement (population-FDR concession + conditional-robustness result) that is stronger than "control fails." The verdict does not rest on the preliminary seed alone.

### T3.3 — orthogonal experimental-PPI referee (DONE; markup #3)
Physical channel independent of AF score AND DepMap: IntAct physical interactions (www.ebi.ac.uk findInteractions), queried for all 233 edge-set proteins (23,229 distinct interacting pairs; 21,051 physical). Degree-matched permutation baseline (10,000×, seed 42), matching on summed IntAct partner-degree to control the hub confound. **Certified edges: 77% have independent physical evidence vs degree-matched null 34% — strongly ENRICHED (perm p<0.0001, matched-null OR 6.2; the raw-background OR 21.9 is inflated by the hub confound the matched null removes).** Dropped edges: 40% vs null 33% — NOT significantly depleted (perm p=0.85). Independence audit (COMPUTED from the IntAct per-pair query, not hardcoded): over the 40-edge per-pair provenance subset (a stricter pair-filtered query than the 101/132 bulk enrichment count), certified-edge evidence spans 5 source databases (IntAct/UniProt/DIP/MINT/bhf-ucl), 13 detection methods (co-IP, Y2H, pull-down, size-exclusion, cross-link — several predating AlphaFold), and 35 distinct PMIDs — not a single CM4AI/Krogan deposition, so not circular. **Honest boundary:** the audit certifies edges disproportionately backed by orthogonal wet-lab evidence, but dropped edges fail STATISTICAL control, not physical reality — some have physical support. Output: data/processed/experimental_ppi_referee.json; module src/emperor/experimental_ppi.py (reproducible from handoff/intact_edges.json cache).

### T3.1 — temporal validation: INFEASIBLE (honest negative, markup #4)
Built a literature-dated temporal split from CORUM 5.3 (historical release files no longer served by the FastAPI backend; complexes dated by earliest supporting-publication year via NCBI E-utilities, all 3689 PMIDs dated, range 1968–2025). Cutoff 2021 (AlphaFold-Multimer release) separates pre-AF baseline (5481 complexes) from AF-era test (145). **Gate NOT MET — infeasible: ZERO of the 1,666 CM4AI candidate pairs are newly-appearing post-2021 CORUM co-complex pairs**, so the ground-truth positive set for the enrichment test is empty (only 12 of 252 proteins in newly-appearing pairs are even in the candidate universe; none of their pairs coincide). Root cause: the CM4AI cell-map screen (690 proteins) doesn't intersect what recent CORUM curation added. Reported as honest infeasibility, not forced. Markup #4's AF-provenance guard is moot (test can't run); had it run, the 145 AF-era complexes were the ready sensitivity-arm exclusion set. **Procedure validation is instead carried by the two channels that DID work: WS2.2 semi-synthetic (validates the METHOD on known truth) + WS3.3 experimental-PPI referee (validates the CERTIFIED SET, matched-null OR 6.2, p<1e-4).** Output: data/processed/timesplit_validation.json.

### T3.2 — KANSL3 demoted to positive control + prospective register frozen (DONE)
KANSL3→MLL1-WDR5/NSL/MSL (CORUM 5386/610/7221) reframed from "flagship discovery" to labeled positive control: the engine recovers KANSL3's known NSL-family membership as a sanity check; KANSL3→MSL(610) has the highest co-essentiality (signed −log10p 19.6) of any nomination. Prospective NOVEL shortlist (KANSL3 excluded) frozen and dated 2026-07-10 to reports/prospective_nominations.md: top 15 conformal-BH certified nominations ranked by orthogonal DepMap co-essentiality (ranking only; firewall intact), each with explicit confirm/refute criteria (co-IP, AF3/Boltz ipTM, future-release appearance). Top novel: YPEL5→WDR26 (CUL4B-DDB1, coess 68.3), CHUK↔IKBKB (IKK complex, 31.8). Gate met: committed dated shortlist + confirm/refute statement. Output: reports/prospective_nominations.md, data/processed/prospective_register.json.

### WS2.1 (GPU second-map null regeneration) — NOT PURSUED (cost decision, user, 2026-07-10)
User's final call: the second real-map audit via Boltz-2/Modal null regeneration is not worth the compute cost (~$50–150 to fold a whole decoy set for a second map). Deliberate scope decision, not an omission. Consequence for the paper: the GENERALITY claim rests on the two channels already validated — WS2.2 semi-synthetic benchmark (validates the METHOD on known truth across the prevalence sweep) and WS3.3 orthogonal experimental-PPI referee (validates the CERTIFIED SET against independent wet-lab data, matched-null OR 6.2). The auditability-boundary finding stands: CM4AI remains the only map shipping a native null, so it is the only externally-auditable one — a point to make positively in framing rather than treat as a coverage gap.

### G-VALID instrument substitution — temporal → leave-one-complex-member-out recovery (pre-registration DEVIATION, 2026-07-12)
**Deviation logged before building.** G-VALID was pre-registered on a *temporal* validation instrument (hold out post-2021 CORUM complexes, test recovery, require p<0.05 + effect size). That instrument is INFEASIBLE on this data — zero of 1,666 candidate pairs are newly-appearing post-2021 CORUM co-complex pairs (see T3.1 append). The *construct* G-VALID measures — "does the nomination procedure recover held-out true complex members better than chance" — is unchanged; only the instrument changes.

**New primary instrument: leave-one-complex-member-out (LOCO) recovery.** For each CORUM complex with ≥2 members in the candidate protein universe, hide one member's membership and ask whether the conformal-certified edge set (q=0.10) still connects it to a remaining member of the same complex (i.e. would nominate it back in). Recovery is compared to an **eligible-impostor null**: non-members that nonetheless have ≥1 candidate edge INTO the complex (so they are eligible to be recovered — this removes the tautological "a random protein has no edge into the complex" confound). Significance by within-complex member/impostor label permutation. (An initial degree-matched-random-protein null was tried and REJECTED as near-tautological: it gave OR≈333 against a ~0.3% null rate, i.e. random proteins essentially never have a within-complex edge, so it measured edge-existence not membership.) Reports member-vs-impostor recovery + permutation p + odds ratio — exactly the G-VALID gate form. This is an instrument substitution, NOT goalpost-moving: LOCO is equal-or-harder (dozens–hundreds of held-out trials vs. the temporal split's zero) and measures the same construct. The masked membership never enters conformal certification (which depends only on ipTM vs decoys), so recovery is genuinely held-out and the purity firewall is untouched (DepMap absent throughout). Pre-registered thresholds: recovery OR ≥ 2 vs matched null AND permutation p < 0.05.

### WS2.1' — second real map (Predictomes) audited WITHOUT GPU (G-GENERAL), 2026-07-12
Instead of regenerating a null by GPU folding (declined WS2.1), audited a second **already-published** real AF-Multimer map that ships its own pair space: the genome-scale human **Predictomes** (Schmid & Walter, Mol Cell 2025; predictomes.org). 20,196 proteins, 1,614,047 pairs — different lab, pipeline, protein universe from CM4AI. Zero GPU cost (28 MB pair-scores table).

**Findings (data/processed/secondmap_audit.json; src/emperor/secondmap.py, reproducible from data/external/predictomes_pair_scores.parquet):**
1. **G-GENERAL met:** the conformal + native-random-null + BH audit transfers to a second real map. Non-degenerate on a RAW axis: SPOC>=0.9 high-conf tier (12,767 edges), certified/dropped {q0.05: 11642/1125, q0.1: 12420/347, q0.2: 12765/2}.
2. **Audit-axis requirement (new methodological finding):** the axis MUST be a raw, untrained structural-confidence readout. Predictomes' headline SPOC is a *trained classifier* (separates interactors from random pairs by construction) → on the SPOC axis the high-conf tier is perfectly separated from a negative null, every conformal p-value pins to 1/(n+1) → degenerate (100% certification). Verified `spoc_axis_degenerate=True`. We therefore audit on `num_unique_contacts` (raw interface size, genuine candidate/null overlap). Explains why CM4AI's raw ipTM was auditable.
3. **Referee non-replication is a signature, not a failure:** DepMap co-essentiality is ALIVE on this map (co-ess fraction 0.08→0.43 across SPOC bins) but does NOT separate certified-vs-dropped within the SPOC>=0.9 tier (certified 0.43 vs dropped 0.50, bg 0.43) — because that tier is already uniformly co-ess-rich. Predictomes is **classifier-curated** (SPOC pre-filters), so its high-conf tier has the junk already removed and the audit finds little to drop. CM4AI's tier is an *uncurated* raw ipTM threshold, which is why the audit catches 22% of it (certified 0.41 vs dropped 0.17). The contrast is itself a finding.
4. **Claim correction:** the paper's "CM4AI is the only auditable AI interactome" is now too strong. Corrected to "CM4AI is the only map shipping a raw, auditable confidence axis IN BULK" — Predictomes ships a pair space but its bulk-downloadable axis is a trained classifier (SPOC); raw ipTM lives only in the ~1.5 TB per-structure archives (not pulled). The ipTM-vs-ipTM cross-map referee comparison is future work needing that access.
