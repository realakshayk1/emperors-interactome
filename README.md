# The Emperor's Interactome

**A distribution-free falsification audit of a celebrated AI protein interactome, refereed by held-out biology.**

Built for the *Built with Claude: Life Sciences* hackathon (Research Track, Gladstone Institutes partner). Submission due **July 13, 2026, 9:00 PM ET**.

---

## The one-paragraph thesis

Structural protein–protein interaction (PPI) maps built from AlphaFold-Multimer are now published as "high-confidence" complexes and used to nominate drug targets. But raw AF-Multimer confidence (ipTM, pDockQ2) is overconfident, and the field's false-discovery rates are *benchmark-estimated lookups* that silently break under prevalence shift. This project applies **distribution-free conformal FDR control** to a published interactome, certifies which complexes survive honest error control, and validates the certified set against **DepMap co-essentiality — evidence the structural model never saw**. On the Krogan/Ideker cancer cell map, we then **nominate a missing member** of a cancer-relevant complex with a certified error rate and a predicted structure.

## Why it wins (judging: Demo 30 / Impact 25 / Claude Use 25 / Depth 20)

- **Demo**: a contrarian, visual story — "the emperor's high-confidence complexes, honestly re-scored" — endable on a live nomination.
- **Impact**: trust layer for AI interactomics + a cancer target nomination.
- **Claude Use**: Claude Science orchestrates literature (PubMed/bioRxiv), pulls structures (PDB/UniProt), runs the pipeline, and writes the interpretation with provenance.
- **Depth**: conformal FDR + a held-out-purity validation design is real statistical rigor.
- **Gladstone Award**: uses Gladstone/Krogan data and executes the "predict a missing component" prompt verbatim.

## Start here (read order)

0. **VERIFY.md** — the verification mandate: question every input in this folder and confirm it firsthand before building on it. Nothing here is ground truth.
1. **HANDOFF.md** — current state + your first move in a fresh Claude Science session.
2. **SPEC.md** — what we're building and the acceptance criteria (the WHAT/WHY).
3. **METHODS.md** — the scientific core: nonconformity score, conformal FDR, held-out validation, nomination.
4. **DATA.md** — exact datasets, download commands, schemas, licenses, ID-mapping plumbing.
5. **TASKS.md** — the day-by-day plan to July 13 (check items off as you go).
6. **PLAN.md** — the HOW: architecture, files, milestones, verification.
7. **DEMO.md** — the frame-by-frame 3-minute demo + narration script.
8. **DECISIONS.md / LEARNINGS.md** — why the project is shaped this way; dead-ends *not* to retry.
9. **AGENTS.md** — operational rules for any coding agent working in this repo.

## Repo layout (target)

```
emperors-interactome/
├── README.md              ← this file
├── AGENTS.md              ← agent operating rules (import into CLAUDE.md)
├── HANDOFF.md             ← context transfer / start-here
├── SPEC.md                ← WHAT/WHY + acceptance criteria
├── PLAN.md                ← HOW: architecture + milestones
├── TASKS.md               ← day-by-day checklist to 7/13
├── DATA.md                ← datasets, downloads, schemas, provenance
├── METHODS.md             ← conformal FDR + held-out validation + physical-validity
├── DEMO.md                ← 3-min demo script
├── DECISIONS.md           ← append-only decision log
├── LEARNINGS.md           ← gotchas + do-not-retry
├── environment.yml        ← pinned conda env
├── requirements.txt       ← pip fallback
├── Makefile               ← `make reproduce` runs the whole pipeline
├── LICENSE                ← MIT (open-source requirement)
├── data/
│   ├── raw/               ← immutable downloads (interactome, CORUM, DepMap, huMAP)
│   ├── interim/           ← ID-mapped, joined tables
│   └── processed/         ← calibration set, certified interactome, results
├── src/emperor/
│   ├── config.py          ← paths, params, seeds, FDR level q
│   ├── download.py        ← fetch + checksum raw data
│   ├── idmap.py           ← UniProt/gene-symbol crosswalk
│   ├── labels.py          ← CORUM/huMAP positives + decoy negatives
│   ├── nonconformity.py   ← ipTM/pDockQ2 + physical-validity score
│   ├── conformal.py       ← conformal p-values + BH FDR control
│   ├── validate.py        ← held-out DepMap co-essentiality enrichment
│   ├── nominate.py        ← missing-member nomination
│   └── plots.py           ← reliability diagram, FDR curve, enrichment, nomination
├── notebooks/
│   └── 1.0-demo.ipynb     ← the demo notebook (also the reproducible narrative)
├── results/figures/       ← generated figures for video + writeup
└── reports/summary.md     ← 100–200 word submission summary
```

## One-command reproduce

```bash
make reproduce   # download → idmap → labels → conformal audit → validate → nominate → figures
```

## Deliverables (hackathon)

- 3-minute demo video (notebook + recorded narration) — see DEMO.md
- Open-source repo (this repo, MIT) with a runnable notebook
- 100–200 word written summary — reports/summary.md

## Config at a glance

- **Primary interactome**: Krogan/Ideker *Nature* 2025 "Multimodal cell maps" (U2OS / cancer), Suppl. Table 5 (1,666 AF-M scored pairs → 111 novel complexes, 275 assemblies).
- **Calibration labels**: CORUM 5.0 (+ hu.MAP 3.0 secondary).
- **Held-out referee**: DepMap co-essentiality (Wainberg 2021 GLS matrix).
- **FDR level**: q = 0.10 (report a sweep).
- **Structure prediction**: Boltz-2 / OpenFold3 **BioNeMo skills** 