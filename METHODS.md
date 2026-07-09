# METHODS — The scientific core

This is the load-bearing document. It specifies the nonconformity score, the conformal FDR
procedure, the physical-validity add-on, the held-out validation design, and the nomination.
Cite the machinery; the novelty is the *first structural-PPI-confidence instantiation + held-out audit*.

---

## 1. Problem setup
We have candidate protein pairs (edges) each with structural confidence metrics from AlphaFold-Multimer.
We want a **distribution-free, FDR-controlled** subset of edges we certify as true interactions, plus a
**calibrated probability** per edge, and an honest demonstration that the certified set reflects real
biology using a channel the model never saw.

Two prediction targets, both handled:
- **Calibration (probabilities):** map raw confidence → P(true interaction). Diagnostic; shows overconfidence.
- **Selection (FDR control):** output a set S of certified edges with FDR ≤ q. The headline result.

---

## 2. Nonconformity score
For edge i, define a scalar nonconformity score `s_i` where **lower = more confident it's a true interaction**.
Base it on raw AF-M interface confidence plus a physical-validity term:

```
s_i = w1 * (1 - iptm_i)
    + w2 * (1 - pdockq2_i)
    + w3 * normalized_interface_PAE_i
    + w4 * phys_penalty_i
```

- Start with equal weights or, better, fit a small monotone/logistic model on the calibration split so
  `s_i` is a learned "improbability of true interaction." Keep it simple and seeded.
- If Table 5 has only one metric (e.g., only pDockQ2), use it alone + the physical-validity term. Don't fabricate metrics.
- `phys_penalty_i` comes from §4.

**Design rule:** `s_i` may use ONLY structural/AF quantities and physical validity. It must NOT use DepMap,
co-expression, GO, or any channel reserved for validation (§5). This is the held-out-purity firewall.

---

## 3. Calibration labels (positives / negatives)
- **Positives:** all within-complex protein pairs from CORUM human core complexes (VERIFIED release 5.3,
  2026-04-14; see DATA.md §2 for the version/host correction) — explode each complex's
  semicolon-separated subunit list into all pairwise combinations). Optionally union with hu.MAP 3.0 high-confidence complexes.
- **Negatives (decoys):** protein pairs not co-occurring in any known complex. Construction matters:
  - Start with **random non-co-complex pairs** matched on degree/abundance where possible, ratio 1:1.
  - Report sensitivity to the positive:negative ratio (1:1, 1:5) — the PPI negative-set problem is real; be explicit.
  - A stricter negative set: pairs whose proteins are each in *some* complex but never the *same* one.
- **Split:** partition labels into a **calibration** set (used to compute conformal p-values) and a held-out
  **CORUM test** set (used only to check empirical FDR ≤ q). Keep complexes disjoint across splits (no leakage of a complex's pairs across cal/test).

---

## 4. Physical-validity add-on (your signature contribution)
Reuse the physical-validity reflex from your `flow-match`/`calibrated-dock` work, transposed to interfaces.
From each deposited AF-M complex structure (the CM4AI map ships 111 high-confidence models; for others
recompute lightly or use interface metrics already tabulated), compute an interface-quality vector:

- **Steric clashes** at the interface (count of atom pairs below VdW threshold; Biopython/`freesasa`).
- **Buried surface area (BSA)** on complex formation — implausibly small BSA ⇒ likely spurious interface.
- **Interface pLDDT** (mean pLDDT of interface residues) and **interface PAE** (mean predicted aligned error
  across the interface) — low-confidence interfaces should be penalized.
- Optional: a strain/energy proxy on interface side chains.

Combine into `phys_penalty_i ∈ [0,1]` (clashes ↑, tiny BSA ↑, high interface-PAE ↑ → larger penalty). This term
is what lets a "high-ipTM but physically impossible" edge be correctly demoted — the crux of the audit.

Tooling: Biopython (`Bio.PDB`), `freesasa` for BSA, numpy for PAE aggregation. PoseBusters is ligand-oriented;
for protein–protein use the clash+BSA+interface-PAE composite above. Metric computation is CPU.

**BioNeMo upgrade (Modal-backed):** if the CM4AI deposited models or their tabulated interface metrics are
incomplete, *generate* the interface structures with the **Boltz-2 / OpenFold3 / AlphaFold2-Multimer BioNeMo
skills** and compute clash/BSA/interface-PAE from those. Budget-limited ($30 Modal ≈ few hundred predictions) —
apply only to the high-leverage subset (the ~111 high-confidence complex interfaces + decision-boundary edges),
NOT the whole map. See DATA.md §5 for the budget. This turns the physical-validity term from a fallback into a
first-class, real-structure signal.

---

## 5. Conformal FDR control (the headline method)
Use **conformal p-values + Benjamini–Hochberg**, following conformal selection / conformal link-prediction FDR
(Marandon 2023, arXiv:2306.14693; Blanchard et al. 2024, arXiv:2404.02542; conformal selection, Jin & Candès 2023).

For each candidate edge j (from the interactome, the "test" pool), and using the calibration **negatives**
{s_k : k in cal, label=0} (size n), compute the conformal p-value:

```
p_j = ( 1 + #{ k in cal_negatives : s_k <= s_j } ) / ( n + 1 )
```

Intuition: p_j is small when edge j looks *more like a true interaction* than the calibration negatives do
(its nonconformity is lower than most negatives). Under exchangeability, p_j is a valid p-value for the null
"j is a non-interaction."

Then apply **Benjamini–Hochberg** to {p_j} at level q → the rejected set S is the **certified interactome**
with FDR ≤ q (distribution-free, marginal). Report:
- the certified set size at q=0.10, plus a sweep q∈{0.05,0.10,0.20};
- how many previously "high-confidence" edges are **dropped** (Branch A headline);
- a comparison to the interactome's own published cutoff / SPOC benchmark-FDR (show ours holds under prevalence shift, theirs doesn't).

**Prevalence-shift centerpiece (the wedge vs SPOC):** benchmark-estimated FDR assumes the test prevalence
equals the benchmark's positive:negative ratio. Re-weight the test pool to a realistic (low) interaction
prevalence and show the benchmark-FDR undershoots true errors while the conformal FDR guarantee still holds
on the held-out CORUM test split. This is the strongest "why conformal, not SPOC" evidence.

**Validity check (unit test):** on synthetic data with known nulls, verify empirical FDR ≤ q across seeds.

---

## 6. Held-out validation design (the referee)
The certified set is only convincing if an *independent* signal agrees. Use **DepMap co-essentiality**
(Wainberg 2021 GLS p-value/sign matrix). Purity firewall: DepMap is used ONLY here, never in §2–§5.

For each edge set (certified@q, raw-high-confidence, dropped-artifacts), compute co-essentiality strength
for its protein pairs (look up GLS value; higher |GLS| / significant = stronger co-essentiality). Then:
- **Enrichment metric:** mean co-essentiality, and AUROC/odds-ratio of "is a strong co-essential pair" for
  certified vs raw vs dropped.
- **Permutation p-value:** shuffle edge-set membership to get a null; report p.
- **Expected result:** certified > raw high-confidence > dropped (certified 2–3× enriched). If the map is
  well-calibrated (Branch B), certified ≈ raw but both > dropped, and the *certified core* is the cleanest —
  still a positive, quotable result.

Purity checklist (put in the notebook): DepMap not in `nonconformity.py`; primary interactome is AP-MS+AF
(DepMap-independent); on Predictomes we audit raw pDockQ2, not SPOC (SPOC ate DepMap). See LEARNINGS.md.

Optional second held-out channel for extra depth: evolutionary coupling / co-conservation (independent of DepMap).

---

## 7. Calibration pre-check (Day-1, decides the headline branch)
Before building the full pipeline, quickly: fit isotonic regression mapping raw ipTM (or pDockQ2) → CORUM label
on the calibration split; plot the reliability diagram on the test split; compute ECE.
- If badly miscalibrated / overconfident → Branch A likely (artifacts exist). 
- If already near-diagonal → Branch B (lead with certified-core enrichment).
Either way keep the reliability diagram — it's demo frame 1.

---

## 8. Missing-member nomination
Pick a **target cancer complex** from the CM4AI map 275 assemblies — ideally one flagged as recurrently mutated
in cancer, with a plausible missing subunit. For each non-member candidate protein c:

```
score(c) = certified_confidence(c ↔ complex members)  ×  held_out_coessentiality(c, complex)
```

- `certified_confidence`: does c have edges to complex members that pass conformal certification (low conf p-value)?
- `held_out_coessentiality`: mean DepMap co-essentiality of c with complex members (the independent vote).
Nominate the top c with certified risk < q AND positive held-out co-essentiality. Output `nomination.json`
and have Claude Science write a literature-grounded rationale (PubMed/bioRxiv connector): what the complex does,
why c is plausible, what's known/unknown. This is the "discrete finding" the research track wants.

---

## 9. Nominee structural validation (now CORE — BioNeMo skill, was Stretch 1)
Because Boltz-2/OpenFold3 are callable BioNeMo skills in Claude Science (Modal-backed), predicting the
{nominated member + top complex partner} structure is a single reliable skill call — promote it from stretch
to **core**. Report ipTM/ligand_iptm + run the §4 physical-validity checks on the predicted interface. A
physically valid, confident predicted interface for a nominee *also* corroborated by held-out co-essentiality
is the wow moment. ~minutes for one complex; ~5–15 predictions well within the $30 Modal budget (DATA.md §5).

## 10. Stretch 2 — second-interactome robustness (Predictomes)
Repeat §2–§6 on Predictomes (~1.6M human pairs; predictomes.org/downloads). **Audit the raw pDockQ2 metric,
NOT SPOC** (SPOC ingests DepMap/co-expression → would break held-out purity). Show the audit conclusion
(overconfidence + certified-core enrichment) generalizes to a second, independently-built map. If the easy
Predictomes CSV lacks pDockQ2, use the small genome-maintenance subset