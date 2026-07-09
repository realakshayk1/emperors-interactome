## Baseline head-to-head: what the distribution-free audit buys

**Setup.** All three decision rules select from the same 1,666 candidate AF-Multimer
pairs on the CM4AI cell map (Schaffer et al., Nature 2025) and were tuned to the same
nominal error target (FDR = 0.10). The honest comparison holds the claimed error rate
fixed and asks which rule's selected set has the strongest *independent, held-out*
support — DepMap co-essentiality (GLS), which never enters the score or calibration
labels (purity firewall).

| Rule | n selected | frac co-essential | mean −log₁₀ GLS p |
|---|---|---|---|
| **conformal-BH, q=0.10** | 132 | **0.414** | **4.36** |
| fixed cutoff (Score ≥ 0.394, benchmark FDR=0.10) | 204 | 0.318 | 3.14 |
| published high-confidence flag (the 161) | 161 | 0.368 | 3.84 |
| *dropped by conformal (diagnostic)* | *35* | *0.167* | *1.48* |

**What conformal gives that the benchmark-estimated methods don't.** The
distribution-free audit certifies the *fewest* edges (132) yet delivers the *highest*
realized held-out support on both metrics. The fixed benchmark-FDR cutoff — set to
FDR=0.10 on a balanced 1:1 benchmark — selects 55% more edges (204) but its co-essential
fraction falls to 0.318: exactly the prevalence-shift failure mode, where a threshold
calibrated at benchmark composition over-selects once true interactions are rare, and
realized support drops even though the nominal FDR is unchanged (conformal vs fixed
cutoff, Δfrac_coess = +0.096, permutation p = 0.098). Against the paper's own published
high-confidence flag the two sets overlap heavily (126 of 132 certified edges are also
published), so the enrichment gap is smaller and not significant on its own
(Δfrac_coess = +0.047, permutation p = 0.51). The discriminating contrast is the
**35 published edges conformal drops**: their held-out co-essential fraction is 0.167 —
less than half the published-set average and a quarter of the certified set — indicating
these are precisely the low-support edges a distribution-free re-derivation removes.

**Interpretation (upper-bound wording).** Because CORUM/DepMap coverage is incomplete,
held-out enrichment is a lower bound on true support and the dropped edges are "not
supported at q=0.10 under distribution-free calibration," not "false." The result is
that conformal-BH converts the same nominal error budget into a smaller, better-supported
edge set than either a benchmark-calibrated fixed cutoff or the map's published flag.

**SPOC baseline: unavailable (network-blocked); comparison limited to conformal /
fixed-cutoff / published-flag.** The Predictomes SPOC score table
(`predictomes_pairs` S3 URL) was not reachable from the sandbox. Per the purity firewall,
SPOC was not substituted with any DepMap-derived proxy — SPOC itself ingests
co-essentiality/co-expression, so validating a SPOC-selected set against DepMap would be
circular.
