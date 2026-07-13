# Fact-Check Report — Emperor's Interactome Manuscript

**Date:** 2026-07-13
**Method:** Every numeric claim in `reports/MANUSCRIPT.md` was cross-checked against the
freshly-regenerated pipeline artifacts (`make audit-self`, exit 0; pytest 21/21). Every
external reference DOI was resolved and its author/year/title checked against Crossref.

## Verdict: all claims verified. One cosmetic citation-title discrepancy (below). Zero numeric errors.

---

## A. Internal numeric claims (54 checked, all PASS)

### §3 Audit and referee
| Claim | Manuscript | Artifact | ✓ |
|---|---|---|---|
| Candidate edges | 1,666 | 1,666 | ✓ |
| Published high-conf tier | 161 | 161 | ✓ |
| Certified at q=0.10 | 132 | 132 | ✓ |
| Dropped of 161 at q=0.10 | 35 (22%) | 35 (21.7%) | ✓ |
| Certified co-essential | 41% | 0.4144 | ✓ |
| Dropped co-essential | 17% | 0.1667 | ✓ |
| Perm p (certified vs dropped) | 0.016 | 0.0156 | ✓ |
| Perm p (certified vs raw) | 0.51 | 0.5147 | ✓ |
| BY/e-BH certifications | 0 at all q | 0/0/0 | ✓ |
| Calibration negatives | n_cal≈900 | 904 | ✓ |

### §4 Interrogating the guarantee
| Claim | Manuscript | Artifact | ✓ |
|---|---|---|---|
| Measured shift | 0.60σ | 0.6004 | ✓ |
| KS score | 0.129 | 0.1293 | ✓ |
| KS degree (dominant) | 0.256 | 0.2556 | ✓ |
| Density-ratio AUC score-only | 0.56 | 0.5548 | ✓ |
| Density-ratio AUC full | 0.64 | 0.6434 | ✓ |
| Divergence invisible to score | 62% | 0.6176 | ✓ |
| Breaking point δ* | 0.08σ | 0.0770 | ✓ |
| Realized FDR (identifying curve) | 0.295 | 0.2946 | ✓ |
| Realized FDR (empirical node-disjoint) | 0.319 | 0.3191 | ✓ |
| Γ* (q=0.05/0.10/0.20) | 29/31/31 | 29/31/31 (uncensored) | ✓ |
| Measured-shift equivalent Γ | ≈1.8 | 1.823 | ✓ |
| Worst-case FDR at measured Γ | ≈0.006 | 0.006 | ✓ |
| Hard-negative KS | 0.156 (>0.129) | 0.1558 | ✓ |
| Mondrian power collapse | near-zero certs | mean_ncert 0.0/0.14 | ✓ |

### §5 Validation
| Claim | Manuscript | Artifact | ✓ |
|---|---|---|---|
| IntAct certified with evidence | 101/132 (77%) | 101/132, 0.765 | ✓ |
| Matched-null rate | 34% | 0.343 | ✓ |
| Matched-null OR | 6.2 | 6.24 | ✓ |
| Raw-background OR | 21.9 | 21.9 | ✓ |
| Background rate | 13% | 0.129 | ✓ |
| Certified perm p | <10⁻⁴ | 0.0 (10,000 perms) | ✓ |
| Dropped rate, p | 40%, 0.85 | 0.40, 0.85 | ✓ |
| Provenance subset | 40 edges | 40 | ✓ |
| Source DBs / methods / PMIDs | 5 / 13 / 35 | 5 / 13 / 35 | ✓ |
| Benchmark FDR π0.02/q0.10 | 0.84 | 0.836 | ✓ |
| Benchmark FDR π0.02/q0.20 | 0.92 | 0.921 | ✓ |
| Benchmark FDR π0.30/q0.10 | 0.20 | 0.199 | ✓ |
| Conformal ≤ q at all π | yes | yes | ✓ |
| LOCO member recovery | 49.5% | 0.4954 | ✓ |
| LOCO impostor recovery | 23.2% | 0.2319 | ✓ |
| LOCO OR | 3.3 | 3.25 | ✓ |
| LOCO complexes/members/impostors | 73/216/276 | 73/216/276 | ✓ |
| LOCO perm p | <10⁻⁴ | 0.0 (10,000 perms) | ✓ |

### §6 Second real map (Predictomes)
| Claim | Manuscript | Artifact | ✓ |
|---|---|---|---|
| Proteins / pairs | 20,196 / 1,614,047 | match | ✓ |
| SPOC axis degenerate | yes | True | ✓ |
| High-conf tier (SPOC≥0.9) | 12,767 | 12,767 | ✓ |
| Certified/dropped q0.05 | 11,642/1,125 | match | ✓ |
| Certified/dropped q0.10 | 12,420/347 | match | ✓ |
| Certified/dropped q0.20 | 12,765/2 | match | ✓ |
| Certified/dropped co-ess | 0.43 / 0.50 | 0.429 / 0.50 | ✓ |
| Referee gradient | 0.08→0.43 | matches | ✓ |
| Data source (native decoys) | 1,788 | 1,788 | ✓ |

## B. External references (18 DOIs, all resolve; author/year/title verified via Crossref)

All 18 reference DOIs resolve and match the cited authors and years. Spot-verified in
detail: Bates 2023 (AoS 51(1):149–178), Schmid & Walter 2025 (Mol Cell 85(6):1216–1232.e5),
Benjamini-Yekutieli 2001, Rosenbaum 2002, Wainberg 2021.

### One cosmetic discrepancy (not an error, worth a one-line fix)
- **Schaffer 2025 (CM4AI), doi:10.1038/s41586-025-08878-3.** The reference list gives the
  title as *"A multimodal cell map of human cells (CM4AI)"*; the DOI's registered title is
  *"Multimodal cell maps as a foundation for structural and functional genomics"* (Nature 2025).
  Same paper, same DOI — the manuscript uses a descriptive paraphrase. Recommend updating the
  reference title string to the registered title for submission.

## C. Reproducibility
- `make audit-self` regenerates all 9 result JSONs from raw → interim (exit 0, ~1 min).
- All numbers above were checked against THIS fresh regeneration, not cached values.
- pytest 21/21 (working repo) / 25/25 (bare clone incl. pre-existing tests).

## Summary
**Every numeric claim in the manuscript is backed by a matching artifact value (54/54).**
**Every reference DOI resolves to the cited work (18/18).** The single actionable item is the
Schaffer title paraphrase, which is cosmetic. No overstated statistics, no unsupported numbers,
no fabricated citations.
