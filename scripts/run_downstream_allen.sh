#!/usr/bin/env bash
# MICA's own node2vec and Louvain on our correlation-3200 graph, then on MICA's own graph, Allen 100K (U-19).
# CPU only, 8 workers, about 80 minutes each. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
name=allen_100k; own=data/runs/rq2/${name}_mica_hnsw
mkdir -p $own && cp results/reproduction/mica_${name}/NNDescent_knn_graph_edgelist.txt $own/edgelist.txt && cp data/runs/rq2/${name}_corr3200/cells.txt $own/cells.txt
for g in corr3200 mica_hnsw; do
  conda activate mica-legacy
  python -u scripts/mica_downstream.py data/runs/rq2/${name}_$g data/runs/rq2/downstream_${name}_${g}_rep1 8
  conda activate scminer
  python -u scripts/score_downstream.py $name data/runs/rq2/downstream_${name}_${g}_rep1 $g rq2_downstream_${name}_${g}_rep1
done
echo "DOWNSTREAM ALLEN DONE"
