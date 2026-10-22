#!/usr/bin/env bash
# PBMC14K under MICA's graph mode, the only dataset above the 5,000 cell switch and the
# only one whose labels do not come from clustering (U-32). About 3 hours per seed.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
for s in 0 1 2; do
  echo "######## pbmc14k graph seed $s"
  python -u scripts/run_rule_stability_graph.py "$s" 25 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return "
done
echo "PBMC GRAPH DONE"
