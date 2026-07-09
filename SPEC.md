# SPEC — The Emperor's Interactome

The WHAT and WHY. Keep tech choices out of here (those live in PLAN.md).

## Summary
Published AlphaFold-Multimer protein interactomes report "high-confidence" complexes using
overconfident raw metrics (ipTM/pDockQ2) and benchmark-estimated FDRs that are not distribution-free.
We re-audit a published interactome with conformal FDR control, certify which complexes survive at a
guaranteed error rate, validate the certified set against held-out DepMap co-essentiality, and nominate
a missing member of a cancer complex with a certified confidence and a predicted structure.

## Goals (measurable)
- G1. Produce a **calibrated** map from raw AF-M confidence to P(true interaction), with a reliability diagram and ECE, on a CORUM-labeled set.
- G2. Produce a **conformal FDR-controlled** certification of the primary interactome at q=0.10 (and a sweep q∈{0.05,0.10,0.20}), splitting edges into certified / dropped.
- G3. Show the certified set is **more enriched for held-out DepMap co-essentiality** than the raw high-confidence set (report effect size + a permutation p-value).
- G4. Deliver **one missing-member nomination** for a named cancer complex, with its certified risk and held-out corroboration.
- G5. `make reproduce` regenerates all results from `data/raw/` with zero manual steps.

## Non-goals (out of scope — do not build)
- Not training a new PPI classifier or structure model (we audit existing outputs).
- Not re-running AlphaFold-Multimer at scale (we use published predictions).
- Not claiming to invent conformal prediction or conformal FDR (cite as machinery).
- Not a proteome-wide rerun — primary analysis is the CM4AI map (Schaffer et al. 2025); Predictomes is a robustness stretch only.
- Not wet-lab validation (held-out computational corroboration only).

## User stories
- As a **structural-interactomics researcher**, I want to know which "high-confidence" complexes survive honest error control, so I don't build on artifacts.
- As a **cancer-biology target-hunter**, I want a missing complex member nominated with a certified error rate and independent corroboration, so I can prioritize follow-up.
- As a **hackathon judge**, I want a live, visual audit that ends in a concrete, trustworthy nomination.

## Functional requirements (numbered, testable)
- R1. Ingest the primary interactome as a table of protein pairs with ≥1 confidence metric. `[NEEDS CLARIFICATION: exact CM4AI Table 5 columns — Day-1 gate]`
- R2. Build calibration labels: CORUM within-complex pairs = positives; matched decoy pairs = negatives; keyed by UniProt accession.
- R3. Compute a nonconformity score per pair from confidence metric(s) + a physical-validity term (interface clashes / buried surface area / interface-PAE from the deposited complex structures).
- R4. Compute conformal p-values for candidate edges using calibration negatives; apply Benjamini-Hochberg to control FDR at q.
- R5. Calibrate raw confidence → probability (isotonic/Platt); output reliability diagram + ECE.
- R6. Held-out validation: for certified vs raw-high-confidence vs dropped edges, compute DepMap co-essentiality enrichment (AUROC and/or odds ratio) + permutation p-value. DepMap must be absent from R2–R5.
- R7. Nominate a missing member: for a target cancer complex, rank non-member candidates by (certified interaction confidence to members) × (held-out co-essentiality), output top nomination with certified risk.
- R8. Stretch-1: predict the nominated member–partner complex with Boltz-2/AF3 (GPU), report ipTM + physical validity.
- R9. Stretch-2: repeat R3–R6 on Predictomes using raw pDockQ2 (not SPOC); report that the audit conclusion holds.

## Data contracts (interfaces between stages)
- `data/interim/interactome.parquet`: columns `[uniprot_a, uniprot_b, iptm, pdockq2, pae_int, phys_valid, source]`.
- `data/interim/labels.parquet`: columns `[uniprot_a, uniprot_b, label∈{0,1}, complex_id, split∈{cal,test}]`.
- `data/processed/certified.parquet`: `[uniprot_a, uniprot_b, nonconf, conf_pvalue, certified@q, prob_calibrated]`.
- `data/processed/validation.json`: enrichment stats per edge-set + permutation p-values.
- `data/processed/nomination.json`: `[complex_id, nominated_uniprot, certified_risk, coess_score, rationale]`.

## Acceptance criteria (the verification gates)
- A1 (G1): reliability diagram exists; ECE reported; raw metric shown mis-calibrated (ECE_raw > ECE_calibrated).
- A2 (G2): `certified.parquet` produced; certified-set empirical FDR on the held-out CORUM test split ≤ q (within conformal tolerance).
- A3 (G3): held-out enrichment of certified > raw high-confidence, with permutation p < 0.05 (report even if effect is modest).
- A4 (G4): `nomination.json` has one member with certified risk < q and positive held-out co-essentiality; a one-paragraph biological rationale is generated (Claude Science + PubMed).
- A5 (G5): fresh clone → `make reproduce` → all figures regenerate; CI-style smoke test passes.

## Test plan
- Unit: conformal p-value validity (coverage on synthetic data), BH FDR control on synthetic nulls, idmap round-trip, decoy generator excludes true positives.
- Integration: end-to-end on a 200-pair subset completes < 2 min CPU and produces all artifacts.
- E2E: `make reproduce` from clean `data/raw/`.
- Sanity: certified-set FDR on held-out CORUM ≤ q; enrichment monotonic-ish as q tightens.

## Open questions
- [NEEDS CLARIFICATION] CM4AI Table 5 confidence columns (Day-1 gate; fallback in DATA.md).
- [NEEDS CLARIFICATION] Target cancer complex choice (pick during Day 3 from the 275 assemblies).
