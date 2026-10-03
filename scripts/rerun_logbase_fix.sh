#!/usr/bin/env bash
# Redo every run touched by the log-base bug found on Oct 3: on the authors' inputs
# for Pollen, Buettner (log10), Chung and Usoskin (natural log), the seurat and
# seurat_v3 gene sets were picked as if the data were log2. All-gene and
# cell_ranger runs are unaffected. Usage: scripts/rerun_logbase_fix.sh
set -uo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
conda activate scminer
python scripts/prepare_paper_inputs.py
declare -A K=([pollen]=11 [buettner]=3 [chung]=4 [usoskin]=4)
: > logs/mica_jobs_5.txt
for d in pollen buettner chung usoskin; do
  python scripts/export_hvg_subsets.py paper_$d 2000 seurat,seurat_v3
  for f in seurat seurat_v3; do echo "paper_${d}_hvg2000_$f ${K[$d]}" >> logs/mica_jobs_5.txt; done
done
xargs -P 8 -L 1 scripts/run_mica_job.sh < logs/mica_jobs_5.txt
for d in pollen buettner chung usoskin; do
  for f in seurat seurat_v3; do
    HVG_FLAVOR=$f python -u scripts/run_baselines.py paper_$d 3 2000
    HVG_FLAVOR=$f python -u scripts/run_label_free.py paper_$d 3 2000
  done
done
scripts/score_all_mica.sh
python -u scripts/bootstrap_gaps.py 1000
python scripts/update_tracker.py && python scripts/summarize_gene_sets.py && python scripts/summarize_label_free.py
echo "RERUN DONE"
