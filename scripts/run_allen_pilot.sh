#!/usr/bin/env bash
# U-16: a 100,000-cell pilot of the Allen atlas, to measure time and memory before
# the full 2.3 million cell run. CPU only. Usage: scripts/run_allen_pilot.sh
# Writes timing and peak memory to logs/allen_pilot_*.log, results to results/.
set -uo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
conda activate scminer
scripts/peak_mem.sh logs/allen_pilot_prepare.log python -u scripts/prepare_dataset.py allen_100k >> logs/allen_pilot_prepare.log 2>&1
scripts/peak_mem.sh logs/allen_pilot_baselines_all.log python -u scripts/run_baselines.py allen_100k 1 >> logs/allen_pilot_baselines_all.log 2>&1
HVG_FLAVOR=seurat_v3 scripts/peak_mem.sh logs/allen_pilot_baselines_hvg.log python -u scripts/run_baselines.py allen_100k 1 2000 >> logs/allen_pilot_baselines_hvg.log 2>&1
# MICA in graph mode, which is its large-data mode, with the paper's resolution sweep
scripts/peak_mem.sh logs/allen_pilot_mica.log scripts/run_mica_job.sh allen_100k ge >> logs/allen_pilot_mica.log 2>&1
echo "pilot done"
