#!/usr/bin/env bash
# The shared-space silhouette rule on the nine small-data datasets (U-35). CPU only.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
for d in paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel; do
  echo "######## $d"
  python -u scripts/run_rule_silhouette.py "$d" 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch"
done
echo "SILHOUETTE DONE"
