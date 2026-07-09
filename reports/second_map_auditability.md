# Second map: auditability is a property of the data release

## Claim

A distribution-free conformal FDR audit is not a universal post-hoc tool. It requires an
**exchangeable null** — a matched negative/decoy set run through the *same* prediction pipeline as
the candidates — plus a confidence axis with residual ranking information. Whether a published
AI-interactome satisfies these is a property of *how the data was released*, not of the map's
scientific quality. We show this by attempting to audit two further public deposits and documenting
exactly where each fails.

## Map 1 (auditable) — CM4AI multimodal cell map

CM4AI ships **1,788 native random-pair decoys** alongside its 1,666 candidate pairs, all scored by
the same AlphaFold-Multimer pipeline (confidence axis = ipTM/`Score`). The decoys are the exchangeable
null. The audit runs cleanly and finds **35 of 161 (22%) high-confidence edges not supported at
q=0.10** under distribution-free calibration.

## Map 2a (not auditable) — Krogan host–pathogen mimicry map (Zenodo 15588019)

The deposit ships AF-Multimer ranked-0 PDB models + FASTA + PyMOL sessions — **no tabulated
confidence scores, and no decoy/random null**. We recovered a confidence axis on CPU by computing
**pDockQ** (Bryant et al. 2022; Cβ–Cβ interface contacts × interface pLDDT) directly from each of the
6,542 ranked-0 dimer models — ipTM is unrecoverable because the PAE matrix was not deposited (that
needs GPU recompute).

The result explains why the map cannot be audited post-hoc:
- **Positives-only.** No exchangeable null was released; the 6,542 pairs are the retained predictions.
- **Saturated axis.** pDockQ median = 0.739, IQR [0.723, 0.742]; **93.1% ≥ 0.5**, only **3.9% < 0.23**;
  the 75th–99th percentile span is **0.000**. The score was already used as the selection filter, so
  the released set sits in the regime where pDockQ carries essentially no remaining ranking
  information. Even a synthesised null could not be discriminated.

A distribution-free FDR audit here would require generating a matched decoy set and *folding it*
(GPU-scale) to obtain a true null — deliberately out of scope on the available compute budget, and a
fabricated (un-folded) null would violate the method's exchangeability assumption.

## Map 2b (not auditable, different reason) — Predictomes / SPOC (genome-maintenance screen)

The reachable Predictomes release (`20251110_hs_predictome_pair_scores.csv.gz`, ~40k pairs) exposes
only the **SPOC** and **KIRC** composite scores plus a raw `num_unique_contacts` count. SPOC is
**firewall-forbidden** as an audit target: it ingests DepMap co-essentiality and co-expression as
input features, so auditing SPOC and validating on DepMap would be circular. No raw pDockQ/ipTM column
and no random-pair null are exposed (bucket listing is disabled; no controls file at any probed key).

## Conclusion

Two independent public AI-interactome deposits, two distinct failure modes — positives-only with a
saturated axis (2a); only a DepMap-contaminated composite exposed (2b) — but the **same** outcome:
post-hoc distribution-free auditability is foreclosed. The single map that *is* auditable (CM4AI) is
auditable precisely because it shipped its native decoys.

**Recommendation for the field:** AI-interactome releases should ship (i) per-pair raw confidence
scores and (ii) a matched random/decoy null run through the same pipeline. Without them, an
independent, distribution-free error-rate audit of the published set is not possible — the reader is
left with the authors' self-reported, prevalence-fragile benchmark FDR.

*(Upper-bound wording throughout: "not auditable post-hoc" means the released artifacts do not support
a distribution-free FDR guarantee, not that the maps' interactions are wrong.)*
