# Self-Improvement Scorecard — PLAN_PUBLICATION_V3 execution

*Generated 2026-07-10. Baseline to beat (frozen before any change): plain-conformal
protein/node-disjoint held-out FDR **0.29 / 0.32 / 0.37** at q = 0.05 / 0.10 / 0.20
(from `wcs_results.node_disjoint_wild`). Every gate below was pre-registered in
`DECISIONS.md` before the run that tested it (deviations flagged where they occurred).*

## Gate outcomes

| Gate | Task | Status | Evidence |
|---|---|---|---|
| **G-CONTROL** | Unconditional FDR control survives real shift | ❌ **honestly fails** → reframed | WS1: realized FDR ≈0.30 at q=0.10 (two routes), δ*≈0.08σ; shift 0.60σ is ~7× past breaking |
| **G-CONTROL (conditional)** | A conditional guarantee stronger than "control fails" | ✅ **met (Path B)** | Certified-set worst-case FDR < q up to Γ*≈31 (grid-crossing); certified core is robust |
| **G-GENERAL (known truth)** | Method controls FDR where truth is known | ✅ **met** | WS2.2 semi-synthetic: conformal FDR ≤ q at all π∈{0.30,0.10,0.05,0.02}; at q=0.10 benchmark cutoff FDR inflates 0.20 (π=0.30) → 0.84 (π=0.02) |
| **G-GENERAL (2nd real map)** | Audit transfers to ≥2 real maps | ✅ **met** | WS2.1′: genome-scale Predictomes (20,196 proteins, 1.6M pairs) audited at $0 GPU; method transfers on a raw axis (SPOC classifier is degenerate → new finding); referee alive but tier classifier-curated |
| **G-VALID (recovery)** | Nomination procedure recovers held-out true members | ✅ **met** (LOCO) | Members recovered 49.5% vs eligible-impostor null 23.2%, OR 3.3, perm p<1e-4; temporal instrument infeasible (0 positives) so substituted with leave-one-member-out |
| **G-VALID (orthogonal)** | Certified set validated by independent data | ✅ **met** | WS3.3 IntAct: certified matched-null OR 6.2 enriched (perm p<1e-4, 101/132 have physical evidence) vs degree-matched null; provenance over a 40-edge per-pair subset spans 5 DBs / 13 methods / 35 PMIDs |
| **G-NOVELTY** | KANSL3 demoted; prospective register frozen | ✅ **met** | WS3.2: dated shortlist + confirm/refute criteria; KANSL3 → labeled positive control |
| **G-INTEGRITY** | Purity firewall + dependence robustness | ✅ **met** | 21/21 tests incl. 5 adversarial firewall + no-signal control + 3 LOCO held-out-recovery; BY/e-BH floor computed |

## What changed vs. the baseline

1. **Dependence-robustness (T0.4).** BH's 132/35 result is now defended by **PRDS**
   (Bates et al. 2023, Ann. Statist. 51(1):149-178, doi:10.1214/22-AOS2244), not
   independence. BY and e-BH (arbitrary-dependence) certify **0** at every q — the
   harmonic price is unpayable at n_cal≈900, m=1666. Reported as the honest floor.

2. **Shift diagnosis (T1.1).** 62% of the calibration→wild divergence is invisible to
   the 1-D score reweighting the earlier WCS attempt used; the dominant unmodeled axis
   is union-graph **degree** (hub selection, KS 0.256). WCS failure is a
   null-completeness problem, not a reweighting bug.

3. **Attempted fixes fail honestly (T1.2/T1.3).** Neither a hard-negative
   (degree-matched) null nor degree-Mondrian conditioning restores unconditional
   control on protein-disjoint splits. Reported, not tuned away.

4. **Sensitivity bound (T1.4).** Two distinct quantities, reported separately (never
   conflated): (A) population realized FDR ≈0.30 under the measured shift (marginal
   guarantee void); (B) the certified *set's* worst-case FDR stays < q up to a large
   tilt Γ*≈30 (the strongly-certified core is robust). An earlier version of this
   module was wrong (worst-case FDR under-counted, false three-way "agreement"); caught,
   corrected, and the false-agreement claim removed.

5. **Semi-synthetic benchmark (T2.2).** The prevalence wedge is now demonstrated on
   **known-truth** data, not assumed: at q=0.10 benchmark FDR inflates to 0.84 at π=0.02
   (and 0.92 at q=0.20) while conformal holds at ≤q throughout. (This is the exchangeable regime, consistent with the
   identifying experiment's δ=0 row; real-data non-exchangeability is the separate WS1
   finding.)

6. **Orthogonal physical referee (T3.3).** New, independent of AF score AND DepMap:
   certified edges are strongly enriched for IntAct physical evidence (matched-null OR 6.2). Honest
   boundary: dropped edges are **not** physically depleted — the audit flags
   statistical-control failures, not proven-false edges.

## Deferred / infeasible (stated honestly)

- **WS2.1 (second auditable map via GPU null regeneration):** **not pursued** — the
  Modal GPU cost (~$50-150 to fold a full decoy set for a second map) was judged not
  worth the marginal value over the channels already validated. Deliberate scope
  decision (2026-07-10), not an omission. The generality claim rests on WS2.2
  (semi-synthetic, validates the method) + WS3.3 (orthogonal referee, validates the
  certified set). The auditability-boundary point — CM4AI is the only map shipping a
  native null, hence the only externally-auditable one — is framed as a positive
  finding, not a coverage gap.
- **WS3.1 (temporal validation instrument):** infeasible on this data (zero positives).
  The G-VALID *construct* (does the procedure recover held-out true members?) is instead
  measured by leave-one-complex-member-out (LOCO) recovery — a denser, feasible
  instrument that meets the pre-registered gate (OR 3.3, p<1e-4); logged as a dated
  pre-registration instrument-substitution deviation in DECISIONS.md. The AF-provenance
  guard (markup #4) is ready but moot until a dataset with candidate↔new-complex overlap
  is available for the temporal route specifically.

## Reproducibility

`make audit-self` runs the entire hardened suite (dependence → shift attribution →
hard negatives → shift-control gate → Γ seed → sensitivity → semi-synthetic → PPI
referee → second real map → LOCO recovery) in ~1 min. `make test` = 21/21.
**Reproducibility check (2026-07-12):** from the release **tarball** (which ships interim
files), `make audit-self` + `make test` reproduce every number with no prior state. From a
bare **git clone** (interim files are git-ignored by design), `audit-self` first regenerates
them via its `make data idmap interactome labels depmap` prerequisite, and `secondmap`
auto-fetches the public Predictomes CSV — so `make reproduce && make audit-self` is the
clean-clone path. Verified the secondmap auto-fetch and the interim prerequisite wiring.
All numbers in this scorecard are emitted to
`data/processed/*.json`.
