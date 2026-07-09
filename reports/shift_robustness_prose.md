## Covariate-shift robustness of the conformal-FDR audit

**The shift is real and significant.** The calibration set (191 CORUM positives + 1,787 native
random decoys) and the wild candidate pool (1,666 AF-Multimer "True" edges being certified) differ
on both axes the marginal conformal guarantee assumes away:

- **Covariate (score) shift.** Comparing the calibration decoys against the wild candidates on the
  AF-Multimer confidence `Score`, the Kolmogorov–Smirnov distance is **D = 0.134 (p = 1.4×10⁻⁹)**.
  The wild pool carries a heavier high-score tail (q90 = 0.437 vs 0.329 for decoys). As a sanity
  check, calibration vs. test decoys give D = 0.026 (p = 0.92) — the decoy population is internally
  exchangeable; the shift is specifically between decoys and wild candidates.
- **Prevalence shift.** Calibration is deliberately enriched (~9.7% positive; the benchmark cutoff
  treats it as roughly balanced), while true interactions are rare in the wild — the classic
  benchmark-FDR failure mode the audit already documented (benchmark realized FDR climbs to 0.45 at
  10% prevalence vs conformal 0.15).

**Weighted Conformal Selection (Jin & Candès 2023) does NOT restore FDR control on this data.** A
logistic density-ratio classifier (calibration-decoy vs wild, on `Score`) gives moderate weights
(max ≈5×, effective-sample-size fraction 0.92). Feeding these into weighted conformal p-values + BH,
with **node-disjoint (protein-disjoint) cal/test splits** so no protein leaks across the split, the
held-out CORUM upper-bound FDR is:

| q | plain conformal | WCS |
|---|---|---|
| 0.05 | 0.29 (14 cert/split) | undefined (0 cert) |
| 0.10 | 0.32 (22 cert/split) | 0.32 (2 cert/split) |
| 0.20 | 0.37 (27 cert/split) | 0.32 (16 cert/split) |

Where WCS actually certifies edges (q=0.10, q=0.20) its held-out FDR is ≈0.32 — essentially
indistinguishable from plain conformal's violation, and well above nominal q. The only q at which WCS
avoids exceeding nominal q is q=0.05, and it does so **vacuously**: it certifies zero edges, so its
FDR is undefined rather than controlled. On the full 1,666-edge pool WCS certifies **zero** edges at
every q ≤ 0.2 (upweighting the high-score decoy tail plus WCS's self-penalty raises the p-value floor
to ≈0.006, and with m=1,666 candidates BH rejects nothing). So WCS trades the FDR violation for a
near-total collapse of power without buying genuine FDR control where it does certify.

**Does the guarantee hold under the real shift? Honest answer: no.** Under protein-disjoint splits the
plain-conformal held-out FDR (0.29–0.37) sits well above nominal q, confirming the exchangeability
assumption is violated by the covariate shift. WCS is a valid *procedure* — becoming conservative
rather than certifying edges it cannot stand behind — but on this data it does not deliver a certified
set at nominal q: it either matches plain conformal's ≈0.32 FDR (where it certifies) or certifies
nothing at all.

**Two important caveats keep this from being fatal to the audit.** (1) The CORUM-based FDR here is an
**upper bound** — CORUM is incomplete, so certified non-CORUM edges are not necessarily false; they
are simply *not supported at q under distribution-free calibration*. (2) The **independent DepMap
co-essentiality referee tells the more favorable story**: the plain-conformal certified set at q=0.10
is co-essentiality-enriched at **frac_coess = 0.41 vs a 0.15 all-wild baseline** — a 2.8× enrichment
that WCS's empty/degenerate set cannot provide. So the plain-conformal certified edges carry real
independent signal even though their CORUM upper-bound FDR exceeds q. The practical recommendation is
to report plain-conformal certifications as *shift-uncorrected* and flank each edge with the DepMap
referee, rather than adopt WCS's zero-power certified set.
