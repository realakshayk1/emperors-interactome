"""Central config: paths, seeds, parameters. Import everywhere; no hidden globals."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "results" / "figures"
for _d in (RAW, INTERIM, PROCESSED, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

SEED = 42

# Conformal / FDR
Q = 0.10                      # primary FDR level
Q_SWEEP = [0.05, 0.10, 0.20]  # report a sweep

# Calibration labels
NEG_POS_RATIO = 1            # decoy negatives per positive; report sensitivity at 5 too
NEG_POS_RATIO_SENSITIVITY = [1, 5]

# Nonconformity weights (start equal; optionally fit on calibration split)
NONCONF_WEIGHTS = dict(iptm=1.0, pdockq2=1.0, pae_int=1.0, phys=1.0)

# Held-out validation
COESS_SIG_THRESHOLD = 0.05   # GLS p-value cutoff for "strong co-essential pair"
N_PERMUTATIONS = 10000

# Nomination target — the cancer complex to nominate a missing member for.
# CHOSEN (not auto-picked): CORUM 5386 "MLL1-WDR5 complex", a leukemia-defining
# H3K4-methyltransferase / chromatin assembly anchored by KMT2A(MLL1) and WDR5.
# It has many conformally-certified edges to non-members AND its members are well
# covered in DepMap. See DECISIONS.md for the rationale and the alternatives ranked.
TARGET_COMPLEX_ID = "5386"

# Data source URLs (see DATA.md for full provenance) — VERIFIED 2026-07 firsthand.
URLS = dict(
    # CM4AI cell map (Schaffer et al., Nature 2025) Suppl. Tables 1-10 workbook, via Europe PMC (PMC12137143).
    cm4ai_supp="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12137143/supplementaryFiles",
    # CORUM 5.3 human complexes via the FastAPI backend (the old .txt.zip host is now a JS SPA).
    corum_human="https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt",
    depmap_gls_p="http://mitra.stanford.edu/bassik/coessentiality/GLS_p.npy",
    depmap_gls_sign="http://mitra.stanford.edu/bassik/coessentiality/GLS_sign.npy",
    depmap_genes="http://mitra.stanford.edu/bassik/coessentiality/genes.txt",
    predictomes_pairs="https://predictomes-hsbps-dataset.s3.us-east-1.amazonaws.com/20251110_hs_predictome_pair_scores.csv.gz",
)

# --- VERIFIED CM4AI (Schaffer et al. 2025) Suppl. Table 5 schema (opened firsthand 2026-07; see DATA.md) ---
CM4AI_SHEET = "Table5"
CM4AI_XLSX_NAME = "NIHMS2086025-supplement-Supplementary_Tables_1_10.xlsx"
# Per-pair AF-Multimer confidence: 5 models of ptm_k / iptm_k, a combined `Score`
# (~mean of AF 0.8*iptm+0.2*ptm), the paper's own benchmark `FDR`, and boolean flags.
# NO pDockQ2 and NO tabulated interface-PAE -> confidence axis is ipTM (`Score`).
CM4AI_CONF_COL = "Score"          # combined AF-M confidence, the nonconformity input
CM4AI_TRUEFLAG_COL = "True vs. random pair"   # True (1666 candidate) vs random (1788 decoy)
CM4AI_HICONF_COL = "High confidence"          # paper's own high-confidence flag (161 True)
CM4AI_NOVEL_COL = "High confidence and novel" # 111 novel complexes
CM4AI_BENCHFDR_COL = "FDR"                     # paper's benchmark-estimated FDR (the wedge target)
