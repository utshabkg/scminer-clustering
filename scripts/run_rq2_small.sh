#!/usr/bin/env bash
# RQ2 prototype (U-19) on the small datasets, CPU only. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
for d in usoskin kolodziejczyk klein zeisel; do
  python scripts/run_rq2_graph.py "$d" 2>&1 | grep -vE "warn|Warn|dispatch"
done
echo "RQ2 SMALL DONE"
