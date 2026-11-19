#!/usr/bin/env bash
# Run one stage of the Allen pipeline with its time and peak memory recorded. CPU only.
#   scripts/large_scale/run.sh <stage> <scope> [workers] [threads]
#   stages: manifest scan hvg normalize pca kmeans leiden assign collect (or "all" for every stage in order)
#   scopes: full, sample100k, ctxsp
# Start it inside a byobu session (scm-large-<scope>), and check uptime, free -g and nvidia-smi first:
# other people's jobs come first on this machine.
set -euo pipefail
cd "$(dirname "$0")/../.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
stage=$1; scope=${2:-full}; workers=${3:-4}; threads=${4:-16}
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=$threads OPENBLAS_NUM_THREADS=$threads MKL_NUM_THREADS=$threads NUMBA_NUM_THREADS=$threads
mkdir -p logs
run() {
  local s=$1 log=logs/large_scale_${scope}_$1.log prefix=""
  case $s in scan|hvg|normalize) prefix="nice -n 5 ionice -c2 -n7" ;; esac
  echo "== $s ($scope) started $(date)" >> "$log"
  scripts/tools/peak_mem.sh "$log" $prefix python -u scripts/large_scale/allen.py "$s" --scope "$scope" \
    --workers "$workers" --threads "$threads" 2>&1 | grep --line-buffered -vE "Warning|warn\(" | tee -a "$log"
  local st=${PIPESTATUS[0]}
  [ "$st" = 0 ] || { echo "FAILED $s ($scope) exit $st $(date)" | tee -a "$log"; exit "$st"; }
}
if [ "$stage" = all ]; then
  for s in manifest scan hvg normalize pca kmeans leiden assign collect; do
    [ -f "data/runs/large_scale$([ "$scope" = full ] || echo "_$scope")/$s.done" ] && { echo "skip $s, done"; continue; }
    run "$s"
  done
else
  run "$stage"
fi
