# HANDOFF — The Emperor's Interactome

**Date:** 2026-07-08 (created, pre-build) · **Updated:** 2026-07-09 (build + extension complete) · **Deadline:** 2026-07-13 21:00 ET

> **STATUS: BUILT AND EXTENDED.** The core audit and all five NEXT_DIRECTIONS extension phases are
> complete, committed, and pushed to `github.com/realakshayk1/emperors-interactome` (main). `make
> reproduce` runs the full CPU pipeline clean end-to-end (~2.5 min); 12 tests pass. The Day-1 gate,
> branch lock, and every "open question" below were resolved during the build — see the resolution
> notes inline and the full record in DECISIONS.md / LEARNINGS.md. The sections below are preserved as
> the original pre-build handoff for provenance; read README.md first for the current state.

**What was actually built (headline results, all verified against committed artifacts):**
- Raw AF-M confidence miscalibrated (AUROC 0.70, ECE 0.18→0.02) → **Branch A locked**.
- **35 of 161 (22%) high-confidence edges fail conformal FDR at q=0.10**; certified {0.05:78, 0.1:132, 0.2:177}.
- Prevalence wedge: benchmark cutoff realized FDR → 0.90 as interactions rare; conformal stays near q.
- Held-out DepMap: certified 41% vs dropped 17% co-essential (perm p=0.016).
- Nominee **KANSL3 → MLL1-WDR5 (CORUM 5386)**, certified risk 0.007, Boltz-2 interface corroboration.
- Extension: Mondrian (structured overconfidence), shift robustness (honest negative), baselines,
  auditability boundary (2nd/3rd maps), cross-arch pilot (conditional GO), FDR-controlled nomination
  sets (393 complexes / 810 nominations), boxed proposition, MLSB/bioRxiv draft, corrected-interactome release.

This doc originally transferred context to *start* building. It is now a build record. Read README.md, then DECISIONS.md.

## Start here (your first move)
0. **Read VERIFY.md and adopt its mandate: question every input in this folder and confirm it firsthand before building on it.** Nothing here is ground truth — it's a plan to verify, not to trust.
1. `@`-reference this folder in Claude Science (or open the repo in Claude Code with `@AGENTS.md` in CLAUDE.md).
2. Create the env (environment.yml) and run `make data` for the small downloads (CM4AI Table 5, CORUM, hu.MAP). These are tiny.
3. **Day-1 gate (do this before anything else):** open CM4AI Suppl. Table 5 and confirm it has per-pair confidence columns (ipTM and/or pDockQ / pDockQ2 / PAE). This one schema fact decides the confidence axis. If present → proceed. If only a pass/fail flag is present → fall back per DATA.md §"Confidence-axis fallback".
4. Run the **calibration pre-check** (METHODS.md §7): is raw ipTM/pDockQ2 miscalibrated vs CORUM? This determines the headline branch (see below).

## Current state
Nothing built yet — this is the initial handoff. All decisions are made and recorded in DECISIONS.md. Datasets are verified available and downloadable (DATA.md). The method is fully specified (METHODS.md).

## The two headline branches (decide on Day 1, keep both figures either way)
- **Branch A ("artifacts"):** if a meaningful fraction of "high-confidence" complexes fail conformal FDR control → headline = "N% of celebrated complexes don't survive honest error control."
- **Branch B ("certified core"):** if the map is well-calibrated → headline = "the certified core is 2–3× more enriched for held-out DepMap co-essentiality than the raw high-confidence set."
Both are strong. You win on either branch — do not let Branch B feel like failure. Report whichever is true; keep both plots.

## Key decisions & rationale (full list in DECISIONS.md)
- **Primary interactome = CM4AI cell map (Schaffer et al., Nature 2025)** (small, clean, cancer-relevant, DepMap-independent). Not Predictomes-as-primary (SPOC uses DepMap → not held-out).
- **Confidence axis = raw AF-M metrics (ipTM/pDockQ2) + physical-validity term**, never SPOC on the primary.
- **Held-out referee = DepMap co-essentiality** (cancer cell lines → matches the U2OS map biology).
- **FDR machinery = conformal p-values + Benjamini-Hochberg** (Marandon 2023 conformal link-prediction FDR). Cited as machinery, not claimed.
- **Scope = core + 2 stretches**: (1) Boltz-2/AF3 structural validation of the nominated member; (2) Predictomes second-interactome robustness (audit raw pDockQ2).

## How to run & verify
- `make reproduce` runs the whole chain. "Green" = every figure in `results/figures/` regenerates from `data/raw/` with no manual step, and acceptance checks in SPEC.md pass.
- Compute division on your 16 GB machine: **everything local/in Claude Science except Boltz-2/AF3 (GPU → Colab/Modal).** DepMap matrix is sliced to interactome genes, so it fits in RAM.

## Open questions / to resolve during build
- [ ] Exact confidence columns in CM4AI Table 5 (Day-1 gate above).
- [ ] Which cancer complex to target for the nomination — pick one of the 275 assemblies that is (a) recurrently mutated in cancer and (b) has a borderline/dropped candidate member with high held-out co-essentiality (METHODS.md §6).
- [ ] Decoy-negative construction ratio for calibration (start 1:1, report sensitivity) — see METHODS.md §3.

## Gotchas already known (full list in LEARNINGS.md)
- Don't materialize the full DepMap matrix (2.5 GB) — slice first.
- Join through UniProt accessions, not gene symbols.
- Predictomes uses UniProt *entry names*; CORUM accessions; DepMap gene symbols.
- Boltz-2's training cutoff means most affinity/structure benchmarks leak — irrelevant here (we don't calibrate on Boltz), but relevant if you extend Stretch 1.

## Next 3 tasks
1. `make data` + verify Table 5 sc
---
## BRANCH LOCK (2026-07-09, build session)
**Headline = Branch A (artifacts exist).** Day-1 calibration pre-check on the held-out
CORUM test split: raw AF-Multimer `Score` is miscalibrated — ECE_raw=0.176 vs
ECE_isotonic=0.016 (AUROC 0.70). Verified schema: CM4AI Table 5 confidence axis = ipTM
(`Score` ≈ mean AF 0.8·ipTM+0.2·pTM); NO pDockQ2/interface-PAE in the table. Native negatives
= 1,788 "random" pairs shipped in Table 5 (same AF-M pipeline, DepMap-independent). Final
Branch-A magnitude = # high-confidence edges dropped by conformal FDR (computed in audit).

## Held-out FDR — report BOTH estimators (no selective reading)
Monte-Carlo over 200 complex-disjoint splits, held-out CORUM label. This is an UPPER
bound on the true FDR (CORUM incompleteness). Report both estimators side by side:
- q=0.10: pooled_FDR=0.157, mean_split_FDR=0.134  (both slightly above nominal)
- q=0.20: pooled_FDR=0.215, mean_split_FDR=0.175  (pooled ~ q, split below q)
So: at q=0.20 the split estimator is controlled and the pooled estimator sits right at
the target; at q=0.10 both run modestly above nominal. The rigorous guarantee is the
synthetic-null unit test (exact control); held-out CORUM corroborates within the
exchangeability/incompleteness caveats, and DepMap is the independent referee.
