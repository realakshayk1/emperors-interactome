# AGENTS.md

Operating rules for any coding agent (Claude Science, Claude Code, etc.) working in this repo.
Keep this file short. If it grows past ~150 lines, move detail into METHODS.md / DATA.md.

> Claude Code note: Claude Code reads `CLAUDE.md`, not `AGENTS.md`. Create a one-line
> `CLAUDE.md` containing `@AGENTS.md` to import this. Claude Science does not confirm
> auto-reading either file — `@`-reference this repo folder at the start of a session,
> or save the pipeline as a Claude Science **skill**.

## Project overview
Distribution-free conformal FDR audit of a published AlphaFold-Multimer protein interactome
(CM4AI cell map (Schaffer et al., *Nature* 2025)), validated against held-out DepMap co-essentiality, ending in a
certified missing-member nomination for a cancer complex. Python, scikit-learn, numpy, pandas,
biopython. All analysis is CPU; structure prediction (nominee validation + physical-validity subset)
uses the **Boltz-2 / OpenFold3 BioNeMo skills in Claude Science, Modal-backed ($30 credit)** — call them
in-app, don't hand off to Colab. Budget is finite: a few hundred predictions; never re-fold the whole map (DATA.md §5).

## Setup commands
- Create env: `conda env create -f environment.yml && conda activate emperor` (or `pip install -r requirements.txt`)
- Configure paths/params: edit `src/emperor/config.py` (data dirs, seed=42, FDR `q=0.10`)

## Build / run / test commands
- Full pipeline: `make reproduce`
- Individual stages: `make data`, `make idmap`, `make labels`, `make audit`, `make validate`, `make nominate`, `make figures`
- Tests: `pytest -q` ; single test: `pytest tests/test_conformal.py -q`
- Lint/format: `ruff check src tests && ruff format src tests`

## Verification mandate (READ VERIFY.md FIRST)
Treat every factual claim in this folder as a hypothesis to confirm, not ground truth — the datasets,
schemas, column names, sizes, licenses, method assumptions, and novelty claims were written by a planning
session and may be stale or wrong. Before any code depends on an input: verify it firsthand (open the file,
print the header/shape), trace numbers to primary sources (never fabricate — Claude Science's reviewer agent
flags untraceable figures), mark unknowns `[UNVERIFIED]`, and try to falsify your own results (permutation
nulls, ablations, swapped labels) before trusting them. Full checklist + standing behavior in **VERIFY.md**.

## The one rule that makes or breaks the science (READ THIS)
**The held-out referee must stay held-out.** DepMap co-essentiality (and any other validation
channel) must NEVER enter the nonconformity score or the calibration labels. On the primary
CM4AI interactome this is automatic (it's an AP-MS + AF map, never trained on DepMap). On the
Predictomes stretch, audit the **raw pDockQ2 metric, NOT the SPOC score** — SPOC ingests DepMap
and co-expression as features, which would make DepMap validation circular. See LEARNINGS.md.

## Code style & conventions
- Python 3.11, `ruff` defaults, type hints on public functions, `snake_case`.
- All randomness seeded via `config.SEED`. No hidden global state.
- Every figure is produced by a function in `src/emperor/plots.py` and written to `results/figures/`.
- Pure functions in `src/`; notebooks only orchestrate + display. Do not bury logic in notebook cells.

## Data rules
- `data/raw/` is **immutable**. Never edit or overwrite a raw download; all cleaning is scripted into `data/interim/` → `data/processed/`.
- Record provenance for every download in DATA.md (URL, date, license, sha256).
- Join proteins through **UniProt accessions**, never gene symbols (alias collisions). See `src/emperor/idmap.py`.

## Guardrails / do NOT
- Do not commit large binaries (`data/raw/*`, `*.npy`, `models/*`) — gitignore them; track by sha256 in DATA.md.
- Do not use any labeled affinity/omics channel that overlaps the confidence score (leakage).
- Do not claim novelty over conformal FDR theory — cite Marandon/Blanchard as machinery; our contribution