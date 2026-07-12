# data/external — external inputs (not all tracked)

## predictomes_pair_scores.parquet  (NOT in git — 42 MB, re-downloadable)
Second real map for the G-GENERAL audit: the genome-scale human **Predictomes**
(Schmid & Walter, *Molecular Cell* 2025, 85(6):1216–1232.e5,
doi:10.1016/j.molcel.2025.01.034). 20,196 proteins, 1,614,047 pairs.

Regenerate:
```
curl -L https://predictomes-hsbps-dataset.s3.us-east-1.amazonaws.com/20251110_hs_predictome_pair_scores.csv.gz -o pair_scores.csv.gz
```
then parse to `predictomes_pair_scores.parquet` (columns: complex_name, uniprot_ids,
acc_a, acc_b, sym_a, sym_b, spoc_score, kirc_score, num_unique_contacts,
string_db_score, biogrid_detect_count, intact_db_evidence_count, in_pdb).
`src/emperor/secondmap.py` consumes this file. Source CSV Content-Length 28,427,577 bytes.

## identifying_experiment.json, wcs_results.json  (tracked)
Inputs to `src/emperor/gamma_seed.py` (Γ* seed): the §19 identifying δ→FDR curve and
the weighted-conformal-shift diagnostics from the extended program.
