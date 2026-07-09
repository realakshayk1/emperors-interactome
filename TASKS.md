# TASKS — day-by-day to July 13, 9 PM ET

Today = Tue Jul 8. Submission = Mon Jul 13, 9 PM ET. ~5.5 working days. Check items off with date/commit.
`[P]` = parallel-safe. Each task names a verify step. Core is D1–D3; stretches D4; polish/ship D5–D6.

## Milestone: core audit certified + validated + one nomination (done-criteria = SPEC A1–A4)

### Day 1 (Tue 7/8) — setup, data, the gate, pre-check
- [ ] T1. Create repo + env (`environment.yml`), scaffold dirs, add MIT LICENSE. — verify: `conda activate emperor && python -c "import emperor"`
- [ ] T2. `download.py`: fetch Krogan Table 5, CORUM 5.0, hu.MAP 3.0, DepMap GLS + genes.txt; write provenance rows. — verify: files in `data/raw/` + sha256 logged
- [ ] T3. **DAY-1 GATE:** open Krogan Table 5; record actual confidence column names in DATA.md. If no numeric confidence → confidence-axis fallback (DATA.md). — verify: column names written down; axis chosen
- [ ] T4. `idmap.py`: build UniProt-accession crosswalk (entryname/symbol ↔ acc). — verify: >95% of interactome + CORUM proteins map; unmapped logged
- [ ] T5. `labels.py`: CORUM positives (explode complexes) + decoy negatives (1:1), complex-disjoint cal/test split. — verify: no positive appears as a negative; splits disjoint by complex
- [ ] T6. Calibration pre-check (METHODS §7): isotonic raw-metric→label, reliability diagram + ECE on test. — verify: `results/figures/reliability_raw.png` exists; **lock headline branch** in HANDOFF.md

### Day 2 (Wed 7/9) — nonconformity + conformal FDR
- [ ] T7. `nonconformity.py`: assemble `s_i` from confidence metric(s) (physical-validity term stubbed = 0 for now). — verify: scores finite; monotone vs label on cal set
- [ ] T8. `conformal.py`: conformal p-values from cal negatives + Benjamini-Hochberg at q. — verify (UNIT): synthetic-null FDR ≤ q across 20 seeds (`tests/test_conformal.py`)
- [ ] T9. Run audit on primary interactome; produce `certified.parquet` + FDR sweep {0.05,0.10,0.20}. — verify: certified-set empirical FDR on held-out CORUM test ≤ q
- [ ] T10. Figure: FDR curve / certified-vs-dropped counts; how many "high-confidence" edges drop. — verify: `results/figures/fdr_curve.png`
- [ ] T11. [P] Prevalence-shift experiment (METHODS §5): show benchmark-FDR breaks, conformal holds. — verify: `results/figures/prevalence_shift.png`

### Day 3 (Thu 7/10) — physical-validity + held-out validation
- [ ] T12. Physical-validity term (METHODS §4): interface clashes + BSA + interface-PAE from deposited models → `phys_penalty`. Re-run audit with it in `s_i`. — verify: a high-ipTM-but-clashing edge gets demoted (spot-check one)
- [ ] T13. `validate.py`: DepMap co-essentiality enrichment for certified vs raw-high-conf vs dropped + permutation p. **Purity check**: DepMap absent from T5/T7/T8/T12. — verify: `validation.json`; enrichment(certified) > enrichment(raw), p<0.05
- [ ] T14. Figure: held-out enrichment bar/box + AUROC. — verify: `results/figures/heldout_enrichment.png`
- [ ] T15. Pick the **target cancer complex** from the 275 assemblies (mutated-in-cancer + a borderline candidate). — verify: complex_id chosen + noted in DECISIONS.md

### Day 4 (Fri 7/11) — nomination + structure validation + robustness
- [ ] T16. `nominate.py`: rank non-members by certified-confidence × held-out co-essentiality; output `nomination.json`. — verify: one nominee, certified risk < q, positive co-ess
- [ ] T17. Claude Science: literature-grounded rationale for the nominee (PubMed/bioRxiv). — verify: paragraph in `nomination.json`, citations resolve
- [ ] T18b. [P, optional] Physical-validity upgrade: generate interfaces via BioNeMo skill for the ~111 high-conf complexes + boundary edges; recompute `phys_penalty` from real structures. Budget-gated. — verify: phys term structure-derived on the subset
- [ ] T19b. Track Modal spend (target < $22 used; A100 not H100; keep buffer). — verify: balance checked
- [ ] T18. **CORE (was Stretch 1):** Boltz-2/OpenFold3 **BioNeMo skill** → structure of {nominee + partner}; ipTM + physical validity. Modal A100, ~5–15 preds. — verify: predicted model + `results/figures/nominee_structure.png`
- [ ] T19. **Stretch 2**: Predictomes robustness — audit raw pDockQ2 (NOT SPOC), repeat validate. — verify: conclusion holds; `results/figures/predictomes_robustness.png` (time-boxed; drop if slipping)

### Day 5 (Sat 7/12) — reproducibility, figures, writeup, demo build
- [ ] T20. `make reproduce` from clean `data/raw/` end-to-end; fix any manual steps. — verify: all figures regenerate, smoke test green
- [ ] T21. Assemble `notebooks/1.0-demo.ipynb` = the narrative (frames per DEMO.md). — verify: runs top-to-bottom < 5 min
- [ ] T22. Write `reports/summary.md` (100–200 words). — verify: word count in range
- [ ] T23. Polish figures (titles, captions, colorblind-safe). — verify: figures legible at video resolution
- [ ] T24. [P] Finalize README + push public GitHub repo (MIT). — verify: repo public, license present, `make reproduce` documented

### Day 6 (Sun/Mon 7/13) — record + submit (buffer)
- [ ] T25. Record 3-min demo video (DEMO.md script) → YouTube/Loom unlisted. — verify: ≤ 3:00, link works
- [ ] T26. Final acceptance pass (SPEC A1–A5). — verify: all checks pass
- [ ] T27. Submit on CV platform before 9 PM ET (v