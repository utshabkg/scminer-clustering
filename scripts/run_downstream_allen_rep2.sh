#!/usr/bin/env bash
# Second repeat of MICA's node2vec and Louvain on both Allen 100K graphs, to see run-to-run spread (U-19).
# Waits for run_downstream_allen.sh to finish first. CPU only, 8 workers, about 70 minutes each.
set -u
cd "$(dirname "$0")/.."
until grep -q "DOWNSTREAM ALLEN DONE" logs/downstream_allen.log 2>/dev/null; do sleep 60; done
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
name=allen_100k
for g in corr3200 mica_hnsw; do
  conda activate mica-legacy
  python -u scripts/mica_downstream.py data/runs/rq2/${name}_$g data/runs/rq2/downstream_${name}_${g}_rep2 8
  conda activate scminer
  python -u scripts/score_downstream.py $name data/runs/rq2/downstream_${name}_${g}_rep2 $g rq2_downstream_${name}_${g}_rep2
done
echo "DOWNSTREAM ALLEN REP2 DONE"
