.PHONY: reproduce data idmap interactome labels calibrate audit prevalence depmap validate nominate figures test lint clean

# Run modules with src on the path — no editable install required (robust across envs).
PY=PYTHONPATH=src python -m emperor

reproduce: data idmap interactome labels calibrate audit prevalence depmap validate nominate figures ## full pipeline, raw -> figures + certified/validation/nomination

data:      ## download raw datasets + provenance (small; DepMap sliced on use)
	$(PY).download

idmap:     ## build UniProt-accession crosswalk
	$(PY).idmap

interactome: ## parse CM4AI Table 5 -> interactome.parquet
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

test:      ## unit + integration tests
	pytest -q

lint:      ## format + lint
	ruff format src tests && ruff check src tests

clean:     ## remove interim/processed (keeps immutable data/raw)
	rm -rf data/interim/* data/processed/* results/figures/*

# Stretch 1 (GPU, run in Colab/Modal, not here):
#   see notebooks or DATA.md §5 — Boltz-2 on {nominee + partner}
