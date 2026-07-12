.PHONY: reproduce data idmap interactome labels calibrate audit prevalence depmap validate nominate figures test lint clean audit-self

# Run modules with src on the path — no editable install required (robust across envs).
PY=PYTHONPATH=src python -m emperor

reproduce: data idmap interactome labels calibrate audit prevalence depmap validate nominate figures ## full pipeline, raw -> figures + certified/validation/nomination

data:      ## download raw datasets + provenance (small; DepMap sliced on use)
	$(PY).download

idmap:     ## build UniProt-accession crosswalk
	$(PY).idmap

interactome: ## parse Krogan Table 5 -> interactome.parquet
	$(PY).interactome

labels:    ## CORUM positives + native decoys, complex-disjoint cal/test split
	$(PY).labels

calibrate: ## Day-1 pre-check: isotonic calibration -> lock Branch A/B
	$(PY).calibrate

audit:     ## nonconformity score + conformal p-values + BH FDR -> certified.parquet
	$(PY).conformal

prevalence: ## prevalence-shift centerpiece (benchmark vs conformal realized FDR)
	$(PY).prevalence

depmap:    ## RAM-safe HTTP-range slice of the DepMap GLS matrix (held-out referee)
	$(PY).depmap

validate:  ## held-out DepMap co-essentiality enrichment (+ permutation p)
	$(PY).validate

nominate:  ## missing-member nomination for the target cancer complex
	$(PY).nominate

figures:   ## regenerate all figures in results/figures/
	$(PY).plots

structure: ## report the nominee's Boltz-2 interface (remote GPU step; see src/emperor/structure.py)
	$(PY).structure

data/interim/labels.parquet:  ## pipeline prerequisites for audit-self (raw -> interim)
	$(MAKE) data idmap interactome labels depmap

audit-self: data/interim/labels.parquet ## PLAN_V3 self-audit (needs interim files: run `make reproduce` first, or this builds them). Dependence-robust FDR, shift attribution, hard-negatives, shift-control gate, sensitivity bound, semi-synthetic benchmark, PPI referee, 2nd real map, LOCO
	$(PY).dependence
	$(PY).shift_attribution
	$(PY).hardnegatives
	$(PY).shift_control
	$(PY).gamma_seed
	$(PY).sensitivity
	$(PY).benchmark_synth
	$(PY).experimental_ppi
	$(PY).secondmap
	$(PY).loco_validation
	@echo "audit-self complete -> data/processed/{dependence_robustness,shift_attribution,hard_negatives,shift_control,sensitivity,benchmark_synth,experimental_ppi_referee,secondmap_audit,loco_validation}.json"

test:      ## unit + integration tests
	pytest -q

lint:      ## format + lint
	ruff format src tests && ruff check src tests

clean:     ## remove interim/processed (keeps immutable data/raw)
	rm -rf data/interim/* data/processed/* results/figures/*

# Stretch 1 (GPU, run in Colab/Modal, not here):
#   see notebooks or DATA.md §5 — Boltz-2 on {nominee + partner}
