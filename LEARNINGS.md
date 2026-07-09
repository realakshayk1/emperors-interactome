# LEARNINGS — gotchas, purity rules, do-not-retry

Highest-value content for a fresh agent. Update as you build. Label dead-ends "DO NOT RETRY."

## Purity rules (violating these silently invalidates the whole result)
- **DepMap is the referee, never a feature.** It must not enter `nonconformity.py`, `labels.py`, or calibration.
  On Predictomes, audit **raw pDockQ2, NOT SPOC** — SPOC ingests DepMap + co-expression, so validating SPOC on
  DepMap is circular. This is the single most important rule in the project.
- **Primary interactome is DepMap-independent** (Krogan AP-MS + AF). Keep it that way; don't "improve" the score with omics.
- **Complex-disjoint cal/test split.** Don't let pairs from the same CORUM complex land in both calibration and
  test — that leaks and inflates coverage.

## Compute / RAM (16 GB machine)
- **DO NOT** `np.load('GLS_p.npy')` fully — it's ~2.49 GB float64 plus working copies → OOM risk. Use
  `mmap_mode='r'` and slice to interactome genes first (DATA.md §3). You only need a few-hundred-gene submatrix.
- Predictomes 116 MB CSV loads fine in RAM; the multi-TB structure archive does NOT — never download it.
- Everything except Boltz-2/AF3 runs locally/in Claude Science. Boltz needs GPU → Colab/Modal; don't attempt on the sandbox.

## ID-mapping traps
- Join through **UniProt accessions**, never gene symbols (alias collisions cause silent misses).
- Predictomes = UniProt **entry names** (`AMBP_HUMAN`); CORUM = accessions; DepMap = gene symbols; Krogan = symbols/acc.
- hu.MAP 3.0 ships dual UniProt+genename columns — use it as a local crosswalk to avoid API calls.
- Log every unmapped ID; a silent 20% dropout will quietly bias the audit.

## Data schema surprises
- CORUM subunits are semicolon-separated multi-accession rows → explode to all within-complex pairs.
- Krogan Table 5 confidence columns are UNVERIFIED until you open it (Day-1 gate). Have the fallback ready (DATA.md).
- Predictomes easy CSV may be SPOC-centric; confirm `pdockq`/`pdockq2` columns before relying on them.

## Statistical gotchas
- Conformal p-values assume **exchangeability**; a PPI network isn't i.i.d. Cite conformal-link-prediction FDR
  (Marandon/Blanchard) which adapts to this; note it as a limitation, don't overclaim.
- The **negative-set problem**: "non-interacting" pairs are only presumed negative. Report sensitivity to the
  positive:negative ratio and to the decoy construction — reviewers will probe this.
- BH controls FDR **marginally**; report it as such.

## Framing / novelty (do not overclaim)
- Prior art that will be cited at you: CalPro (conformal on AF *structure*/pLDDT), SPOC/Predictomes (classifier +
  benchmark-FDR), AHMPC (uses omics as *prior*, not held-out). Your novelty = falsification re-audit + held-out
  corroboration loop + raw-metric conformal FDR on structural PPI. Lead with those; cite the rest as machinery.

## DO NOT RETRY
- Do not make Predictomes the primary and validate on DepMap — purity broken (SPOC ate DepMap). Settled in DECISIONS.
- Do not run Boltz-2 on the local sandbox — no GPU; wastes hours. Use Colab/Modal.
- Do not treat a well-calibrated map (Branch B) as failure — the certified-core enrichment is a real, quotable win.

<!-- Append findings during the build. -->

## Build-session verified findings (2026-07-09)
- **CORUM host + version drift.** Old `mips.helmholtz-munich.de/corum/.../coreComplexes.txt.zip` is now a
  JS SPA returning a 423-byte HTML shell. Real data = FastAPI on `mips.helmholtz-muenchen.de`
  (`-muenchen`): `/fastapi-corum/public/releases/current` → **5.3 / 2026-04-14**; download via
  `/public/file/download_current_file?file_id=human&file_format=txt`. This is CORUM **5.3**, not the 5.0
  the plan assumed — schema is compatible (dual `subunits_uniprot_id` / `subunits_gene_name` columns).
- **CORUM TLS chain is broken from the sandbox** (`CERTIFICATE_VERIFY_FAILED`). `download.py` disables
  verification for that one public CC-BY file only and size/schema-checks the payload after. Do NOT
  generalize the bypass to other hosts.
- **Krogan Table 5 confidence axis = ipTM.** Columns: `ptm_0..4`, `iptm_0..4`, combined `Score`
  (~mean AF 0.8·ipTM+0.2·pTM), the paper's own `FDR`, and boolean flags. NO pDockQ2, NO tabulated
  interface-PAE → single-metric path (METHODS §2). `phys_penalty` stays 0 until BioNeMo structures wired.
- **Native negatives.** Table 5 ships 1,788 "random" pairs alongside 1,666 True candidates — same AF-M
  pipeline, DepMap-independent. Use these as calibration negatives (cleaner exchangeability than ad-hoc decoys).
- **nature.com is NOT reachable from the sandbox** (proxy 403 even after allowlist grant). Fetch the Krogan
  supplement from Europe PMC instead: `ebi.ac.uk/europepmc/webservices/rest/PMC12137143/supplementaryFiles`.

## Phase 2-3 verified findings (2026-07-09)
- **DepMap GLS matrix sliced via HTTP range, never downloaded.** The `.npy` is Fortran-ordered
  (column-contiguous); parallel range requests fetch only interactome-gene columns → 1,049×1,049
  slice in ~100s, ~9 MB. Header verified: `<f8`, shape (17634,17634), data offset 80. TP53-MDM2
  p=4e-18. See `src/emperor/depmap.py`.
- **Held-out validation is the discriminating test.** Certified edges are 41% co-essential vs 17%
  for edges dropped by conformal FDR (permutation p=0.016). Certified-vs-raw-high-conf is NOT
  significant (p=0.51) because certified ⊂ raw — always compare certified vs DROPPED, not vs raw.
- **Nomination self-confirms.** KANSL3 → MLL1-WDR5 (CORUM 5386) recovers a protein that is an
  annotated NSL-complex member (CORUM 7221) but absent from the target entry — a database-internal
  positive control that never touched the pipeline. WDR5's NSL vs MLL/COMPASS interactions are
  mutually exclusive (Dias 2014), which is exactly why an AF-M screen places the NSL scaffold KANSL3
  adjacent to the MLL1-WDR5 entry. Honest scope: essentiality-relevant, NOT a Cancer Gene Census driver.
- **Jupyter kernels cannot bind sockets in this sandbox** (`zmq ... Operation not permitted`, TCP and
  IPC both fail). To produce an executed .ipynb, run cells in the live kernel and inject outputs into
  the notebook JSON via nbformat — `jupyter nbconvert --execute` / nbclient will not start a kernel here.
- **`make` is not preinstalled** in the conda env; `pip install`/conda `make`. Also: the editable
  install can be shadowed by a stray `.venv/emperor`; `pip install -e .` into the active env fixes
  `ModuleNotFoundError: emperor`. The demo notebook adds `src/` to `sys.path` so it runs regardless.
