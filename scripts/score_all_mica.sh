#!/usr/bin/env bash
# Score every finished MICA run. Small-data (mds) runs first get three seeded
# k-means draws on MICA's own embedding, so they compare fairly with the baselines.
# Usage: scripts/score_all_mica.sh
set -uo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
cat logs/mica_jobs_0.txt logs/mica_jobs.txt logs/mica_jobs_2.txt logs/mica_jobs_3.txt 2>/dev/null | sort -u | while read -r name mode; do
  out=results/reproduction/mica_$name
  [[ -d $out ]] && ls $out/clustering_*.txt >/dev/null 2>&1 || { echo "skip $name (no output)"; continue; }
  if [[ $mode != ge ]]; then conda activate mica-legacy; python scripts/mica_seeded_kmeans.py "$out" "$mode"; fi
  conda activate scminer
  if [[ $name == *_hvg2000_* ]]; then ds=${name%%_hvg2000_*}; tag=hvg2000_${name#*_hvg2000_}
    python scripts/score_mica.py "$ds" "$out" "$tag"; python scripts/mica_label_free.py "$ds" "$out" "$tag"
  else python scripts/score_mica.py "$name" "$out"; python scripts/mica_label_free.py "$name" "$out"; fi
done 2>&1 | grep -vi warn
conda activate scminer
# Zeisel in graph mode, which is how the authors' own script ran it
for n in zeisel paper_zeisel; do
  if [[ -d results/reproduction/mica_${n}_ge ]]; then
    python scripts/score_mica.py $n results/reproduction/mica_${n}_ge --exp reproduction_${n}_ge 2>&1 | grep -vi warn
    python scripts/mica_label_free.py $n results/reproduction/mica_${n}_ge --exp labelfree_mica_${n}_ge 2>&1 | grep -vi warn
  fi
done
python scripts/update_tracker.py && python scripts/summarize_gene_sets.py
