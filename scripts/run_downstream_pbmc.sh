#!/usr/bin/env bash
# MICA's own node2vec and Louvain on our fast graph and on MICA's own graph, PBMC14K (U-19).
# CPU only, 8 workers, about 40 minutes. Run inside byobu.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
name=paper_pbmc14k; fast=data/runs/rq2/${name}_corr800; own=data/runs/rq2/${name}_mica_hnsw
conda activate scminer
python -u scripts/export_fast_graph.py $name 800 8
mkdir -p $own && cp results/reproduction/mica_${name}/NNDescent_knn_graph_edgelist.txt $own/edgelist.txt && cp $fast/cells.txt $own/cells.txt
for rep in 1 2; do
  for g in corr800 mica_hnsw; do
    conda activate mica-legacy
    python -u scripts/mica_downstream.py data/runs/rq2/${name}_$g data/runs/rq2/downstream_${name}_${g}_rep$rep 8
    conda activate scminer
    python -u scripts/score_downstream.py $name data/runs/rq2/downstream_${name}_${g}_rep$rep $g rq2_downstream_${name}_${g}_rep$rep
  done
done
echo "DOWNSTREAM DONE"
