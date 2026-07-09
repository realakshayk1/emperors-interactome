#!/bin/bash
# Batch Boltz-2 fold of the 12 pilot dimers. Each FASTA staged flat into workdir.
# Uses the free ColabFold MSA server (api.colabfold.com) — only the fold is billed.
set -e
mkdir -p out
for fa in *.fasta; do
  name="${fa%.fasta}"
  echo "=== folding $name ==="
  # boltz writes to a run dir; point it at a per-pair subdir, then copy confidence JSON to out/
  boltz predict "$fa" --use_msa_server --out_dir "run_$name" --output_format pdb --override 2>&1 | tail -3 || echo "FOLD_FAILED $name"
  # harvest the confidence json + top model
  cj=$(find "run_$name" -name 'confidence_*.json' | head -1)
  pdb=$(find "run_$name" -name '*_model_0.pdb' -o -name '*.pdb' | head -1)
  if [ -n "$cj" ]; then cp "$cj" "out/${name}_confidence.json"; echo "  saved confidence for $name"; fi
  if [ -n "$pdb" ]; then cp "$pdb" "out/${name}_model.pdb"; fi
done
echo "=== all folds done ==="
ls -la out/