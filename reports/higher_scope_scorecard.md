# Higher-scope VERIFY — scorecard (MASTER_DIRECTIONS Phase 2)

Four directions probed in parallel. Honest status: **D2 and D3 are decisive KILLs at G1/G2; D1 is
blocked on a data-access approval (with a novelty headwind); D4 is the one live candidate but its G3
informativeness probe did not execute.** This matches the master doc's expectation that most directions
KILL early — the discipline is to learn that from cheap probes before committing.

| Dir | G1 novelty | G2 data | G3 signal probe | G4 feasibility | Ceiling | GO/KILL |
|---|---|---|---|---|---|---|
| **D1** Perturb-seq target discovery | ⚠ **threatened** — the Marson 2025 paper *itself* already nominates novel regulators and links them to autoimmune GWAS, which is exactly D1's intended headline | ⏸ blocked — DE-stats live on `genome-scale-tcell-perturb-seq.s3.amazonaws.com` (off-allowlist, **grantable**); probe parked on your approval card | not run (data gated) | CPU-feasible once data granted | good journal *if* a non-headline corroborated regulator emerges | **PENDING** (needs approval; G1 headwind) |
| **D2** Structure-mechanism (Boltz-2) | ✗ **KILL** — method saturated (AlphaMissense, RosettaDDG, structure-based ΔΔG-of-binding); SOD1 dimer-destabilization mechanism is textbook (2004–2019) | ✓ reachable — 185 SOD1 ClinVar missense (117 clean pathogenic) but **~0 clean benign missense** (negatives only a gnomAD-rare proxy) | not run — froze writing fold inputs; even a clean signal would be a known method | fits budget (~$2.2 for WT+10) | — | **KILL** (at G1) |
| **D3** Noncoding-variant (ChromBPNet) | ✗ **KILL** — saturated: ChromBPNet + fine-mapping + motif-disruption is the published Corces 2023 (Cell) workflow; calibration angle covered (seq-to-activity uncertainty, 2023) | ✗ **KILL** — ChromBPNet 5-fold weights on **denylisted S3** (ungrantable); Borzoi substitute needs ≥24 GB VRAM (> $5 budget) and gives no cross-fold calibration | moot (G1+G2 both fail) | not feasible on this stack | — | **KILL** (G1 & G2) |
| **D4** Benchmark-gap (PoseBusters move) | ✓ **PASS** — genomic-FM calibration benchmark is threatened (ACMG variant-effect calibration mature: Pejaver 2022/ClinGen), but **RNA/nucleic-acid physical-validity benchmark ("PoseBusters for RNA"): no direct competitor found** — least-filled | ✓ reachable — RNA-Puzzles standardized dataset (24 native crystals + 1000 group predictions, all real PDBs), 927 parsed on CPU | ✓ **RAN — informative.** MolProbity-style clashscore discriminates native vs prediction (KS D=0.50, p=1.5e-4). Surprising: natives are NOT clash-free (median 9.6); predictions bimodal — **33% over-idealized (cleaner than any native), 14% physically implausible (up to 15,000+ clashes), only 53% in the native band** | CPU-feasible; ~1 min for 927 structures | high adoption / good venue | **GO** |

## Recommendation

The portfolio is **mostly KILL as predicted**, with **one GO: D4**. Two directions are dead on arrival
(D2 method-saturated, D3 data-ungrantable). D1 remains pending on a data-access approval with a novelty
headwind.

- **D4 — GO (winner).** The "PoseBusters for RNA/nucleic-acid physical validity" gap survived a hard
  novelty search (no direct competitor) AND its G3 informativeness probe returned a surprising,
  actionable result: on 927 real RNA-Puzzles structures, a MolProbity-style clashscore cleanly separates
  natives from predictions (KS D=0.50, p=1.5e-4), and reveals that predictions are bimodal — a third are
  *over-idealized* (energy-minimized cleaner than any real crystal) and a seventh are physically broken.
  A "fewer clashes = better" benchmark would be wrong; the informative benchmark is a **native-calibrated
  validity band**. This is a real, publishable methods/resource finding with a high ceiling.
- **D1 — parked, not pursued.** Given D4 is a clean GO and D1 carries a genuine novelty headwind (the
  Marson paper already nominates novel GWAS-linked regulators), D1 was not retried.

**Decision (per user):** run D4's probe (done — GO) AND execute the interactome §A identifying experiment
as the guaranteed elevating deliverable. Both shipped.

**Guardrails honored:** no G3 result is claimed for any probe that did not execute one (D2/D3 KILLs rest
on G1/G2, not on unrun probes); validators never entered any score; attribution stays CM4AI / Schaffer et
al. 2025; the RNA-validity finding is stated as native-calibrated (not "fewer clashes = better").