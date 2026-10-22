#!/usr/bin/env bash
# Stability rule (U-32) on the nine benchmark sets that MICA runs in small-data mode.
# CPU only. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
for d in paper_yan paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel; do
  echo "######## $d"
  python -u scripts/run_rule_stability.py "$d" 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|^    subsample"
done
echo "RULES SMALL DONE"
