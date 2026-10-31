#!/usr/bin/env bash
# The stability comparison on matched candidate grids (U-40), one part and one seed per call:
#   scripts/run_rules_matched.sh small <seed>      the nine small datasets, about an hour
#   scripts/run_rules_matched.sh pbmc <seed>       PBMC14K in graph mode, about 2 to 3 hours
#   scripts/run_rules_matched.sh modesize <seed>   the mode test on one 2,500 cell pool
#   scripts/run_rules_matched.sh silhouette 0      the shared-space silhouette rule, nine small datasets
# The grids and the rule for choosing candidates are in src/scminer_clust/grids.py. CPU only.
# Run each call in its own byobu session.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
part=$1 s=$2
case $part in
  small)
    for d in paper_yan paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel; do
      echo "######## matched stability seed $s $d $(date)"
      python -u scripts/run_rule_stability.py "$d" "" seurat "$s" --matched 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|^    subsample"
    done ;;
  pbmc)
    echo "######## matched pbmc14k graph seed $s $(date)"
    python -u scripts/run_rule_stability_graph.py "$s" 25 --matched 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return " ;;
  silhouette)
    for d in paper_yan paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel; do
      echo "######## matched silhouette $d $(date)"
      python -u scripts/run_rule_silhouette.py "$d" "" seurat --matched 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch"
    done ;;
  modesize)
    echo "######## matched mode test pool seed $s $(date)"
    python -u scripts/run_rule_mode_size.py 2500 "$s" 25 --matched 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return " ;;
esac
echo "MATCHED $part SEED $s DONE $(date)"
