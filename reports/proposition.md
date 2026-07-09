# Motivating proposition (framing, not a contribution)

> The elementary observation below is textbook Bayes plus a known Benjamini–Hochberg property
> (Bates et al. 2023). We state it only to *motivate* why a distribution-free audit is needed. It is
> **not** claimed as a theorem or a contribution of this work.

## Setup

A structural-proteomics pipeline ranks candidate protein pairs by a confidence score `s` (here
AlphaFold-Multimer ipTM) and reports a **benchmark-estimated FDR**: it fixes a cutoff `t` so that on a
*balanced* benchmark (equal true and decoy pairs) the false-discovery rate is some target `q₀`. Let

- `π` = the true prevalence of interacting pairs in the population the cutoff is applied to,
- `π₀` = the prevalence in the balanced benchmark used to calibrate `t` (here `π₀ = 1/2`),
- `TPR(t)`, `FPR(t)` = the score rule's true- and false-positive rates at cutoff `t`.

## Proposition (two halves)

**(a) Positive half — conformal + BH is prevalence-invariant.** Conformal p-values are calibrated
against an exchangeable null, so under exchangeability the p-values of true nulls are (super-)uniform
*regardless of `π`*. Benjamini–Hochberg on those p-values controls FDR at the nominal level under
independence or PRDS (Benjamini–Yekutieli 2001; Bates et al. 2023). The guarantee does not reference
`π`, so it does not move when prevalence does.

**(b) Negative half — a benchmark-tuned cutoff inflates with the prevalence gap.** The realized FDR of
the fixed cutoff `t` at true prevalence `π` is

```
            (1 − π) · FPR(t)
FDR(π) = ───────────────────────────────
          (1 − π) · FPR(t) + π · TPR(t)
```

Calibrated at `π₀` to give `q₀`, the same `t` applied at `π < π₀` yields a realized FDR that rises
monotonically as `π → 0`, with the inflation governed by the odds ratio `[(1−π)/π] · [π₀/(1−π₀)]`. As
`π → 0` the realized FDR → 1 for any fixed `t` with `FPR(t) > 0`. The benchmark FDR is a
prevalence-fragile *lookup*, not a guarantee.

## Empirical confirmation (already in this repo)

Applying the map's own benchmark-tuned cutoff (`Score ≥ 0.394`, calibrated to "FDR = 0.10" on a
balanced set) to test pools of decreasing true prevalence reproduces (b) exactly, while conformal + BH
tracks (a):

| true prevalence π | benchmark-cutoff realized FDR | conformal + BH realized FDR |
|---|---|---|
| 0.50 | 0.082 | 0.051 |
| 0.30 | 0.173 | 0.083 |
| 0.20 | 0.258 | 0.129 |
| 0.10 | 0.447 | 0.154 |
| 0.05 | 0.633 | 0.000 |
| 0.02 | 0.816 | 0.000 |
| 0.01 | 0.900 | 0.000 |

The benchmark cutoff's realized error climbs to **0.90** as interactions become rare; conformal + BH
stays near the nominal `q` with finite-sample slack and, at very low prevalence, **certifies nothing
rather than exceed `q`** (a conservative failure, not an anti-conservative one). This gap is the wedge
the audit is built on — and the reason a benchmark-estimated FDR on a sparse interactome cannot be
taken at face value. (Figure: `results/figures/prevalence_shift.png`.)
