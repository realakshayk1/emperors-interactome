# DECISIONS — append-only decision log

ADR-style. Newest at bottom. Mark superseded entries rather than deleting. Captures *why* the project
is shaped this way so a fresh agent doesn't re-litigate settled choices.

## 2026-07-08 — Project = falsification audit, not a trust-layer tool
- **Status:** accepted
- **Context:** The obvious idea (recalibrate AF-Multimer confidence) is crowded and reads as tool-building.
- **Decision:** Reframe as an *audit* of a published interactome — certify what survives honest FDR, refereed by held-out biology.
- **Rationale:** Contrarian + falsifiable is more novel and more memorable; wins on either result branch.
- **Alternatives rejected:** Boltz-2 affinity trust layer (methods-not-biology, GPU-bound, no Gladstone data); Krogan "obvious" PPI trust layer (the idea anyone would generate).

## 2026-07-08 — Primary interactome = CM4AI cell map (Schaffer et al., Nature 2025)
- **Status:** accepted
- **Context:** Need a published, high-profile structural interactome that is DepMap-independent and cancer-relevant.
- **Decision:** Use the Nature 2025 "Multimodal cell maps" (U2OS) Suppl. Table 5 as primary.
- **Rationale:** Small/clean (1,666 pairs, 111 complexes), cancer context matches DepMap's cancer cell lines, and it is an AP-MS+AF map so DepMap stays genuinely held-out. CM4AI is a Bridge2AI consortium co-led by Krogan (Gladstone) → indirect Gladstone link, Award-eligible.
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
- **Status:** accepted. CM4AI Table 5 exposes ipTM (`iptm_0..4`), pTM (`ptm_0..4`), combined `Score`
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

## 2026-07-09 — Program extension (NEXT_DIRECTIONS Sessions 1–5)

- **Primary-map attribution corrected.** Krogan/Ideker → **CM4AI cell map (Schaffer et al., Nature 2025, doi:10.1038/s41586-025-08878-3)**; Ideker senior author, Krogan co-leads the CM4AI/Bridge2AI consortium (indirect Gladstone link), NOT a paper author. Renamed `KROGAN_*`→`CM4AI_*`, provenance `cm4ai_schaffer2025_table5`; pipeline reproduces identically.
- **Named audit (Pillar A).** The 22% headline made concrete: 35 dropped high-conf edges named; 7 pair proteins co-occur in the same CORUM complex. Upper-bound wording enforced repo-wide ("not supported at q", never "false").
- **Mondrian stratified conformal.** Overconfidence is STRUCTURED: pooled conformal holds in medium complexes (held-out FDR 0.08), is marginally above target in small (0.12, wide CI), and BREAKS in large 9-member complexes (0.65, >3× target). Supporting rigor, not a headline.
- **Covariate-shift robustness = HONEST NEGATIVE.** Cal→wild shift is real (KS 0.134, p=1.4e-9). On node-disjoint splits BOTH plain and weighted conformal (WCS, Jin&Candès 2023) exceed nominal q where they certify; WCS collapses power without an FDR gain. The guarantee does NOT survive the real shift — the contribution is the diagnosis, not a rescue. Rigorous guarantee remains the synthetic-null unit test.
- **Baselines.** At the same nominal FDR=0.10, conformal-BH selects fewest (132) yet has highest held-out support (41%) vs fixed cutoff (204, 32%) and published flag (161, 37%). SPOC never audited (ingests DepMap → circular).
- **Second map = AUDITABILITY BOUNDARY.** Krogan host-pathogen deposit (Zenodo 15588019) is positives-only, pDockQ saturated (median 0.74, 93%≥0.5), NO native null → not auditable post-hoc without GPU decoy folding. Predictomes exposes only SPOC (DepMap-contaminated) + contact count, no null. Finding: auditability is a property of the data RELEASE. Recommendation: interactome releases must ship their nulls. pDockQ recovered on CPU from ranked-0 PDBs (Bryant 2022; validated on 452 benchmark dimers).
- **Cross-architecture pilot (Pillar B) = GO with caveat.** Boltz-2 corroborates the conformal verdict on 12 pilot dimers: certified ipTM 0.71 vs dropped 0.46 (Mann-Whitney p=0.021); AF-M vs Boltz Spearman 0.80. Caveat: certification is near-separable from AF-M score on this map, so this shows AGREEMENT (drops aren't AF-M-specific artifacts), not yet independent signal beyond score. GO to develop as a nonconformity feature; funded scale-up is the path.
- **Nomination as FDR-controlled selection.** Generalized the single KANSL3 nomination to per-complex certified SETS: 393 complexes, 810 nominations, all inheriting FDR≤q=0.10 (upper bound). Selection purely conformal; DepMap ranks only. Third-map audit SKIPPED (candidate maps paywalled / no confirmed null); documented as future work.
- **Write-up.** MLSB/bioRxiv draft skeleton (`reports/DRAFT_mlsb.md`) + boxed motivating proposition (`reports/proposition.md`, framing not theorem) + corrected interactome released as browsable CSV (`data/processed/corrected_interactome.csv`).
