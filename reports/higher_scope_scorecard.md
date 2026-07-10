# Higher-scope VERIFY — scorecard (MASTER_DIRECTIONS Phase 2)

Four directions probed in parallel. Honest status: **D1, D2, and D3 are all KILLs; D4 is the sole GO.**
D2/D3 die at G1/G2; D1 was resumed after the user granted the CZI S3 access and ran its full G3 signal
probe on real data — it KILLs on executed evidence (no autoimmune-genetics enrichment). This matches the
master doc's expectation that most directions KILL — the discipline is to learn that from cheap probes.

| Dir | G1 novelty | G2 data | G3 signal probe | G4 feasibility | Ceiling | GO/KILL |
|---|---|---|---|---|---|---|
| **D1** Perturb-seq target discovery | ⚠ **PASS on method-per-se; finding likely not novel** — no published calibration-gated confident-hit pipeline on this dataset, but the source paper's own verified title couples its regulators to human immune traits and the release ships the cross-guide/cross-donor reproducibility columns the confident-hit rule reuses (paper full text NOT fetched — internal claims left unverified) | ✓ **reachable (access granted)** — CZI public S3; read the full `.obs` (33,983×28) out of the remote 16.8 GB `GWCD4i.DE_stats.h5ad` via h5py+fsspec HTTP range reads in ~3 s / <300 MB RAM, no download, no GPU | ✓ **RAN — KILL.** Confident-hit rule on all 33,983 target×condition tests → 1,480 certified rows → 878 regulators. Open Targets corroboration across 11 autoimmune diseases: 9.7% of confident vs 8.8% of all tested genes — **Fisher OR=1.13, p=0.32, NO enrichment**. **Zero** regulators pass novel-AND-corroborated; the 85 corroborated hits are known immune/TCR biology or housekeeping machinery | CPU-feasible (range reads) | KILL — reconstructs the paper's existing confidence + immune-trait design | **KILL** (at G3) |
| **D2** Structure-mechanism (Boltz-2) | ✗ **KILL** — method saturated (AlphaMissense, RosettaDDG, structure-based ΔΔG-of-binding); SOD1 dimer-destabilization mechanism is textbook (2004–2019) | ✓ reachable — 185 SOD1 ClinVar missense (117 clean pathogenic) but **~0 clean benign missense** (negatives only a gnomAD-rare proxy) | not run — froze writing fold inputs; even a clean signal would be a known method | fits budget (~$2.2 for WT+10) | — | **KILL** (at G1) |
| **D3** Noncoding-variant (ChromBPNet) | ✗ **KILL** — saturated: ChromBPNet + fine-mapping + motif-disruption is the published Corces 2023 (Cell) workflow; calibration angle covered (seq-to-activity uncertainty, 2023) | ✗ **KILL** — ChromBPNet 5-fold weights on **denylisted S3** (ungrantable); Borzoi substitute needs ≥24 GB VRAM (> $5 budget) and gives no cross-fold calibration | moot (G1+G2 both fail) | not feasible on this stack | — | **KILL** (G1 & G2) |
| **D4** Benchmark-gap (PoseBusters move) | ✓ **PASS** — genomic-FM calibration benchmark is threatened (ACMG variant-effect calibration mature: Pejaver 2022/ClinGen), but **RNA/nucleic-acid physical-validity benchmark ("PoseBusters for RNA"): no direct competitor found** — least-filled | ✓ reachable — RNA-Puzzles standardized dataset (24 native crystals + 1000 group predictions, all real PDBs), 927 parsed on CPU | ✓ **RAN — informative.** MolProbity-style clashscore discriminates native vs prediction (KS D=0.50, p=1.5e-4). Surprising: natives are NOT clash-free (median 9.6); predictions bimodal — **33% over-idealized (cleaner than any native), 14% physically implausible (up to 15,000+ clashes), only 53% in the native band** | CPU-feasible; ~1 min for 927 structures | high adoption / good venue | **GO** |

## Recommendation

The portfolio is **three KILLs and one GO (D4)** — the master doc's predicted outcome. D2/D3 are dead on
arrival (method-saturated / data-ungrantable). D1, after the user granted the CZI S3 access, ran its full
G3 probe and KILLs on executed evidence.

- **D4 — GO (winner).** The "PoseBusters for RNA/nucleic-acid physical validity" gap survived a hard
  novelty search (no direct competitor) AND its G3 informativeness probe returned a surprising,
  actionable result: on 927 real RNA-Puzzles structures, a MolProbity-style clashscore cleanly separates
  natives from predictions (KS D=0.50, p=1.5e-4), and reveals that predictions are bimodal — a third are
  *over-idealized* (energy-minimized cleaner than any real crystal) and a seventh are physically broken.
  A "fewer clashes = better" benchmark would be wrong; the informative benchmark is a **native-calibrated
  validity band**. This is a real, publishable methods/resource finding with a high ceiling.
- **D1 — KILL on executed evidence.** With S3 access granted, the calibration-gated confident-hit rule
  ran on all 33,983 target×condition tests (→ 878 confident regulators). Corroboration against Open
  Targets across 11 autoimmune diseases showed **no enrichment** (9.7% vs 8.8% background, Fisher OR=1.13,
  p=0.32), and **zero** regulators passed novel-AND-corroborated: the 85 corroborated confident hits are
  textbook immune/TCR biology or housekeeping machinery. The calibration "unlock" reconstructs the source
  paper's existing confidence + immune-trait design rather than exposing a gap. Evidence:
  `D1_confident_hits_corroborated.csv`, `D1_G3_probe.png`. (Caveat: paper full text was not fetched; the
  KILL rests on the executed null result, not on the unverified internal claims of the source paper.)

**Decision (per user):** run D4's probe (done — GO) AND execute the interactome §A identifying experiment
as the guaranteed elevating deliverable. Both shipped. D1 was subsequently resolved (KILL) after the user
granted data access.

**Guardrails honored:** no G3 result is claimed for any probe that did not execute one (D2/D3 KILLs rest
on G1/G2, not on unrun probes); validators never entered any score; attribution stays CM4AI / Schaffer et
al. 2025; the RNA-validity finding is stated as native-calibrated (not "fewer clashes = better").