# VERIFY — question every input before you build on it

**Standing mandate (read before Day 1, applies every session).** This whole folder was written from a
planning/research session. Treat its factual claims as *hypotheses to confirm*, not ground truth. Datasets
move, schemas differ, papers appear, and the planner may be wrong. Your job is to **independently verify
each load-bearing assumption before code depends on it**, and to surface — not silently absorb — anything
that doesn't check out.

## How to behave
- **Verify, then build.** For every dataset, URL, column name, size, and license used by a stage, confirm it
  firsthand (open the file, print the header, check the shape) before writing logic that assumes it.
- **Cite primary sources.** Don't restate a claim from these docs as fact — trace it to the actual paper/file.
  Claude Science's reviewer agent flags untraceable numbers; lean on it, and never fabricate a figure.
- **Mark uncertainty loudly.** Use `[UNVERIFIED]` / `[NEEDS CLARIFICATION]` inline; keep a running list here.
- **Adversarial self-check.** After each stage, ask: "What would make this result wrong or leaked? What's the
  most likely bug? Would a skeptical reviewer believe the number?" Write the answer in LEARNINGS.md.
- **Escalate contradictions.** If a verified fact contradicts a DECISION, stop and flag it (update DECISIONS.md);
  don't quietly route around it.
- **Prefer falsification.** Actively try to break your own certified set and nomination (permutation nulls,
  swapped labels, ablations) before trusting them. A result that survives your own attack is demo-ready.

## Assumptions to confirm (checklist — each was asserted by the planner, none is guaranteed)

### Data facts
- [ ] CM4AI cell map (Schaffer et al., *Nature* 2025) Suppl. Table 5 actually contains per-pair numeric confidence (ipTM / pDockQ / pDockQ2 / PAE). **This is the Day-1 gate.** Record real column names in DATA.md.
- [ ] The 111 deposited AF-M complex models are downloadable (for the physical-validity term). Find the exact repository.
- [ ] Predictomes flat CSV (`..._pair_scores.csv.gz`) column set — does it carry raw `pdockq2`, or is it SPOC-only? Print the header before relying on it.
- [ ] DepMap GLS files load and are 17,634² with `genes.txt` in matching order; identifiers are HGNC symbols.
- [ ] CORUM 5.0 current download URL on the moved host (mips.helmholtz-munich.de) still serves `coreComplexes.txt` and it's CC BY 4.0.
- [ ] hu.MAP 3.0 files + dual UniProt/gene-name columns are as described.
- [ ] Every dataset's license still permits open-source redistribution of derived results (hackathon requires open source).

### Method / statistics assumptions
- [ ] Conformal p-values assume **exchangeability**; confirm the calibration↔test setup satisfies it (or use the network-adapted conformal-link-prediction FDR variant). State the limitation.
- [ ] BH controls FDR **marginally** under the assumed dependence; sanity-check with a synthetic-null unit test (empirical FDR ≤ q across seeds).
- [ ] The **negative set** is only *presumed* non-interacting — verify results are stable across decoy construction + pos:neg ratio.
- [ ] **Held-out purity**: independently confirm the CM4AI map did NOT use DepMap as an input, and that SPOC DOES use DepMap/co-expression (this is why we audit raw pDockQ2 on Predictomes). If either is wrong, the validation design must change — flag immediately.

### Novelty / prior-art (re-check; the field moves weekly)
- [ ] Re-run a quick prior-art search at build time for "conformal FDR AlphaFold-Multimer interactome / held-out
      co-essentiality audit." The planner found it unoccupied (~July 2026) but with close neighbors (CalPro,
      SPOC, AHMPC, conformal-link-prediction FDR). Confirm nothing new scooped it; adjust framing, don't overclaim.

### Environment / Claude Science capabilities (verify in-app, some were UNVERIFIED)
- [ ] Whether Claude Science auto-reads AGENTS.md/CLAUDE.md (planner could not confirm) — if not, `@`-reference this folder or save a skill.
- [ ] Which connectors are actually available to you (PubMed/bioRxiv/PDB/UniProt) for the literature-grounded rationale + ID mapping.
- [ ] The **BioNeMo skills** (Boltz-2 / OpenFold3 / AlphaFold2-Multimer) are actually callable in your Claude Science session, and produce ipTM/PAE outputs you can parse — test one small prediction before building on it.
- [ ] Your **Modal balance is ~$30** and billing is A100 (~$2.50/hr), not H100. Estimate cost per prediction empirically on the first run; extrapolate before batching. Watch cold-start/MSA overhead.
- [ ] Local RAM headroom on the 16 GB machine when the DepMap slice + working copies are live.

## Definition of "verified enough to ship"
`make reproduce` runs clean from `data/raw/`; 