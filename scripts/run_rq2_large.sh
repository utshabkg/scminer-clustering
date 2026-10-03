#!/usr/bin/env bash
# RQ2 prototype (U-19) on PBMC14K with sampled brute force. CPU only, 8 workers. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
python -u scripts/run_rq2_large.py paper_pbmc14k 8 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return fn"
echo "RQ2 LARGE DONE"
