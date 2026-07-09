# PLAN — the HOW

Architecture, files, milestones, verification. Generated after SPEC is stable. Keep simple:
fewest components, use libraries directly, no speculative abstraction.

## Technical approach
A linear, cached pipeline of small Python modules, each reading/writing a parquet/json data contract
(SPEC "Data contracts"). Orchestrated by a Makefile so any stage re-runs from its inputs. Claude Science
drives it (writes code, runs cells, renders figures, pulls literature); GPU step offloads to Colab/Modal.

- **Language/stack:** Python 3.11; numpy, pandas, pyarrow, scikit-learn (isotonic/logistic), scipy,
  statsmodels (BH), biopython + freesasa (interface metrics), matplotlib. `ruff` + `pytest`.
- **Why conformal p-values + BH (not a bespoke model):** distribution-free, tiny code, provable FDR,
  and it's the exact wedge vs SPOC's benchmark-FDR. Traced to SPEC R4/G2.
- **Why parquet contracts between stages:** each stage independently re-runnable + inspectable; supports the
  `make <stage>` targets and caching on a 16 GB machine.

## Files to create (responsibilities)
- `src/emperor/config.py` — paths, `SEED=42`, `Q=0.10`, `Q_SWEEP=[0.05,0.10,0.20]`, ratios.
- `src/emperor/download.py` — fetch raw data + checksum + provenance rows (SPEC R1, T2).
- `src/emperor/idmap.py` — UniProt-accession crosswalk (T4).
- `src/emperor/labels.py` — CORUM/huMAP positives + decoys + cal/test split (R2, T5).
- `src/emperor/nonconformity.py` — `s_i` from AF metrics + physical-validity (R3, T7, T12).
- `src/emperor/conformal.py` — conformal p-values + BH; `certify(interactome, cal_neg, q)` (R4, T8–T9).
- `src/emperor/validate.py` — held-out DepMap enrichment + permutation (R6, T13).
- `src/emperor/nominate.py` — missing-member scoring (R7, T16).
- `src/emperor/plots.py` — all figures (reliability, FDR curve, prevalence shift, enrichment, structure).
- `tests/` — conformal validity, BH FDR control, idmap round-trip, decoy purity.
- `notebooks/1.0-demo.ipynb` — orchestrates + displays (no logic).

## Milestones
- **M1 (end D1):** data in, IDs mapped, labels built, reliability pre-check → branch locked. Gate: A1 partial.
- **M2 (end D2):** conformal FDR audit runs, certified set produced, unit tests green. Gate: A2.
- **M3 (end D3):** physical-validity in the score + held-out validation positive. Gate: A3.
- **M4 (end D4):** nomination + rationale + Stretch 1 (Boltz) + Stretch 2 (Predictomes). Gate: A4.
- **M5 (end D5):** `make reproduce` clean, demo notebook, summary, public repo. Gate: A5.
- **M6 (D6):** video recorded, submitted. Gate: all A pass.

## Simplicity gates (check before adding anything)
- Could a plain function replace this class? Use it.
- Are we adding a dependency the stdlib/sklearn already covers? Don't.
- Is this feature in SPEC? If not, it's out of scope.

## Verification (the loop)
`make reproduce` = download → idmap → labels → audit → validate → nominate → figures. "Green" = all
figures regenerate from `data/raw/` + `pytest -q` passes + acceptance checks A1–A5 in SPEC pass.
