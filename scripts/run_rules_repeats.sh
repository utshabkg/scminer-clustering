#!/usr/bin/env bash
# Seeds 1 and 2 of the stability rule, plus MICA's shipped default and the shared-space
# silhouette rule. Waits for the seed 0 batch to finish first. CPU only. Run in byobu.
set -u
cd "$(dirname "$0")/.."
until grep -q "RULES SMALL DONE" logs/rules_small.log 2>/dev/null; do sleep 30; done
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
SETS="paper_yan paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel"
echo "######## MICA shipped default"
python -u scripts/run_rule_default.py 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch"
for d in $SETS; do
  echo "######## silhouette $d"
  python -u scripts/run_rule_silhouette.py "$d" 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch"
done
for s in 1 2; do
  for d in $SETS; do
    echo "######## stability seed $s $d"
    python -u scripts/run_rule_stability.py "$d" "" seurat "$s" 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|^    subsample"
  done
done
echo "RULES REPEATS DONE"
