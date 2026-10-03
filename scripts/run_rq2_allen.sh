#!/usr/bin/env bash
# RQ2 at 100,000 cells (U-19), CPU only, 8 workers, about 2 to 4 hours. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
scripts/peak_mem.sh logs/rq2_allen.log python -u scripts/run_rq2_allen.py 8 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return fn" >> logs/rq2_allen.log
echo "RQ2 ALLEN DONE" >> logs/rq2_allen.log
