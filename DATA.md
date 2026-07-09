# DATA — datasets, downloads, schemas, provenance

All datasets verified available and downloadable as of 2026-07. `data/raw/` is immutable; record
URL + date + license + sha256 for each download in the provenance table at the bottom.

Compute note (16 GB RAM): every file below fits comfortably except the DepMap matrix, which you
**slice to interactome genes before loading** (see §3). No cloud needed except Stretch-1 GPU.

---

## 1. Primary interactome — CM4AI multimodal cell map (Schaffer et al., *Nature* 2025)
- **Attribution (VERIFIED firsthand):** Schaffer LV, Hu M, Qian G, … Lundberg E, Ideker T.
  "Multimodal cell maps as a foundation for structural and functional genomics." *Nature* 642,
  222–231 (2025). Senior author **Trey Ideker** (UCSD) — **Krogan is NOT an author**. The map is a
  product of **CM4AI** (Bridge2AI consortium, co-led by Krogan/UCSF-Gladstone, Lundberg, Ideker), the
  indirect Gladstone/Krogan link. Cite as the "CM4AI map (Krogan-co-led consortium)," never "Krogan's map."
- Paper: https://www.nature.com/articles/s41586-025-08878-3 (open access; PMID 40205054)
- Data: **Supplementary Table 5** (AF-Multimer predictions for 1,666 protein pairs → 111 high-confidence
  novel complexes; 275 assemblies). Download the supplementary `.xlsx` from the Nature article page.
  Also mirrored on CM4AI/NDEx (cm4ai.org).
- **Deposited complex structures** (for the physical-validity term): the 111 AF-M models — locate via the
  paper's data-availability section (Zenodo/CM4AI). If per-pair PAE isn't tabulated, recompute interface
  PAE from the model files.
- License: open-access article supplementary. Format: Excel → convert to parquet in `download.py`.
- **DAY-1 GATE:** open Table 5 and confirm it has per-pair confidence columns. Expected: some of
  {ipTM, pDockQ / pDockQ2, PAE, interface residues}. Record the actual column names in this file.

### Confidence-axis fallback (if Table 5 lacks numeric confidence)
1. Recompute pDockQ2 / interface-PAE from the deposited AF-M model files (small; 1,666 pairs).
2. Or use the **Predictomes genome-maintenance 40,176-pair set** as the primary (has pDockQ) and make
   the CM4AI map the robustness check instead. Update DECISIONS.md if you switch.

## 2. Calibration labels — CORUM (VERIFIED: release **5.3**, 2026-04-14)
- **Schema/version drift found at build (2026-07).** The planner assumed CORUM 5.0 at
  `mips.helmholtz-munich.de/corum → coreComplexes.txt.zip`. That host is now a JS single-page
  app (the `.txt.zip` URLs return a 423-byte HTML shell, confirmed firsthand). The real data is
  served by a FastAPI backend on a **different domain**: `https://mips.helmholtz-muenchen.de`
  (`-muenchen`, not `-munich`). `GET /fastapi-corum/public/releases/current` returned
  `{"version":"5.3","date":"2026-04-14","is_current":true}` — so the downloaded release is **5.3**,
  not 5.0. Download: `/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt`.
- **TLS note:** that host presented a broken certificate chain from the sandbox
  (`CERTIFICATE_VERIFY_FAILED`); `download.py` disables verification for this one public CC-BY file
  and size/schema-checks the payload afterward (6.27 MB TSV, 7,867 human complexes).
- Fields (verified): `complex_id, complex_name, organism, cell_line, pmid, subunits_uniprot_id
  (';'-separated), subunits_gene_name (';'-separated), ...` — ships BOTH accessions and gene names,
  so CORUM doubles as the local ID crosswalk (no hu.MAP needed).
- 7,867 human complexes; 6.27 MB; License **CC BY 4.0**.
- Use: explode subunit lists → all within-complex UniProt-accession pairs = positive labels.
- Secondary: **hu.MAP 3.0** — https://humap3.proteincomplexes.org/download →
  `hu.MAP3.0_complexes_wConfidenceScores_total15326_wGenenames_20240922.csv` (cols: `HuMAP3_ID, ComplexConfidence(1–6), Uniprot_ACCs, genenames`). CC0. Ships BOTH UniProt + gene names → use as the ID crosswalk bridge.

## 3. Held-out referee — DepMap co-essentiality (Wainberg 2021)
- Precomputed GLS matrix (fastest path), no login, Stanford host:
  - `http://mitra.stanford.edu/bassik/coessentiality/GLS_p.npy`  (17,634 × 17,634 GLS p-values)
  - `http://mitra.stanford.edu/bassik/coessentiality/GLS_sign.npy` (+1/−1 sign)
  - `http://mitra.stanford.edu/bassik/coessentiality/genes.txt`   (17,634 gene symbols, row/col order)
  - Repo: https://github.com/kundajelab/coessentiality
- License: academic-open. Identifiers: **HGNC gene symbols** (17,634).
- **VERIFIED .npy header (2026-07, HTTP range read):** `{'descr':'<f8','fortran_order':True,
  'shape':(17634,17634)}`, data offset 80 bytes, total 2,487,663,728 bytes; the host supports
  `Accept-Ranges: bytes` (206 partial content).
- **RAM rule (critical on this 8 GB sandbox) — implemented in `src/emperor/depmap.py`:** we NEVER
  download the 2.49 GB array. The `.npy` is **Fortran-ordered**, so each column is contiguous on
  disk. We HTTP-range-fetch only the columns for interactome genes (parallel threads), then subset
  the same rows → a **1,049 × 1,049 submatrix (~9 MB)** written to `data/interim/depmap_slice.npz`.
  Of 1,159 interactome genes, 1,049 (90%) are in DepMap; the rest are symbol-version mismatches
  (CARS1←CARS, DCAF1←VPRBP, …), left unscored. Sanity-checked: matrix symmetric, zero diagonal,
  TP53–MDM2 GLS p=4e-18 (sign −1, anti-correlated — correctly NOT counted as positive co-essential).
  (The `mmap_mode='r'` + `np.ix_` slice below is the equivalent recipe IF the file is downloaded
  locally first; we avoid the download entirely.)
  ```python
  # local-file equivalent (NOT used — we range-fetch instead, see depmap.py):
  genes = [l.strip() for l in open('data/raw/depmap_genes.txt')]
  idx = {g:i for i,g in enumerate(genes)}
  keep = [idx[g] for g in interactome_genes if g in idx]
  M = np.load('GLS_p.npy', mmap_mode='r'); sub = np.asarray(M[np.ix_(keep, keep)])
  ```
- Alt (build-your-own): DepMap `CRISPRGeneEffect.csv` from figshare (24Q2:
  https://plus.figshare.com/articles/dataset/DepMap_24Q2_Public/25880521), correlate columns. CC BY 4.0. Heavier; only if you want a fresher matrix.

## 4. Stretch-2 interactome — Predictomes / SPOC
- Downloads: https://predictomes.org/downloads
- Small file: `20251110_hs_predictome_pair_scores.csv.gz` (28 MB gz / 116 MB; all 1.6M human pairs; open S3, CC BY 4.0).
  Direct: `https://predictomes-hsbps-dataset.s3.us-east-1.amazonaws.com/20251110_hs_predictome_pair_scores.csv.gz`
- **Check the header for `pdockq`/`pdockq2`/`avg_models` columns.** Audit raw pDockQ2, NOT SPOC (LEARNINGS.md).
- IDs are UniProt **entry names** (`AMBP_HUMAN`) — map to accessions before joining.
- If the flat CSV is SPOC-only, use the 40,176-pair genome-maintenance set (SBGrid dataset 1155, CC0).

## 5. Structure prediction — BioNeMo skills (Boltz-2 / OpenFold3 / AlphaFold2-Multimer) via Modal
Claude Science natively exposes **Evo 2, Boltz-2, and OpenFold3** as callable BioNeMo Agent Toolkit skills;
the broader BioNeMo NIM catalog adds **AlphaFold2-Multimer, ESMFold, DiffDock, RFdiffusion**. Call them as
skills in-app — no Colab handoff. Backend compute = **your Modal account** (the hackathon/Claude Science plan
includes **$30 free Modal credit**).
- Boltz-2 outputs: `iptm`, `ligand_iptm`, `plddt`, PAE, `affinity_probability_binary` (affinity not needed here).
- OpenFold3 / AlphaFold2-Multimer: protein-protein complex structures + confidence — use for the
  physical-validity term and the nominee validation.
- Fallback if a skill isn't wired in your session: Boltz weights (MIT) github.com/jwohlwend/boltz +
  Colab github.com/JKourelis/Colab_Boltz-2.

### Modal budget ($30) — spend it surgically
- **Rate:** A100-80GB = $2.50/hr ($0.000694/sec); H100 = $3.95/hr. **Use A100, not H100.** $30 ≈ **12 A100-hrs**.
- **Throughput:** ~1–2 min/complex incl. MSA/overhead → **a few hundred predictions total.** Cold starts + MSA
  dominate; batch, and keep a ~$5–8 buffer.
- **Allocation (suggested):**
  - Nominee validation (nominee + top partners): ~5–15 structures.
  - Physical-validity term: targeted subset only — the ~111 high-confidence complex interfaces + decision-boundary edges (~100–300 structures). **Do NOT re-fold the full interactome** (over budget, and the thesis i
## 6. Structure prediction — nominee interface (VERIFIED, executed 2026-07-09)
- **Tool:** Boltz-2 v2.2.1, architecturally independent of the audited AlphaFold-Multimer.
- **Where:** the user's own Modal account (`byoc:modal`), GPU **A100-80GB**. Image
  `im-HUmTG4YGPRXCeHPD6lJUXV` (bundled env `proteomics_boltz_gpu@fc3f36390514cc77`); Boltz-2
  weights hydrated into the persistent Volume `claude-science-boltz-cache` (`boltz2_conf.ckpt`,
  `boltz2_aff.ckpt`, CCD mols).
- **Input:** full-length KANSL3 (UniProt Q9P2N6, 904 aa) + KANSL1 (Q7Z3B3, 1105 aa) = 2,009-residue
  dimer. MSA generated once via ColabFold (`api.colabfold.com`), cached, then reused so the scored
  run is inference-only (~5 min).
- **Command:** `boltz predict kansl3_kansl1_precomputed.yaml --out_dir ./out --output_format pdb
  --cache /root/.boltz --num_workers 2` (YAML supplies both chains + per-chain precomputed MSA CSVs).
- **Result:** global ipTM 0.469, pTM 0.426; a localized, confidently-placed interface
  KANSL3[279–448]–KANSL1[624–690] (90 inter-chain CB–CB<8Å & PAE<10Å contacts, interface pLDDT
  0.70/0.59, best-contact PAE 4.8 Å). Global confidence is moderate — depressed by the large
  disordered regions of both proteins — so this is honest corroboration of the certified edge, not a
  solved complex. Outputs: `data/structures/kansl3_kansl1_boltz2_model.pdb`,
  `..._confidence.json`, `physical_validity.json`; figure `results/figures/nominee_structure.png`;
  summary folded into `data/processed/nomination.json → physical_validity`.
- **GPU-tier note:** the Modal plan gates GPUs ≥40 GB behind a payment method; 24 GB tiers
  (A10G/L4/T4) OOM on this dimer. A payment method was added to unlock A100-80GB.
