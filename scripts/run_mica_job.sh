#!/usr/bin/env bash
# Run the authors' MICA on one prepared input. Usage: scripts/run_mica_job.sh <name> <k|ge>
# <name> is an exchange folder name without _exchange, e.g. zeisel or zeisel_hvg2000_seurat.
# An optional third argument suffixes the run name, e.g. "ge" to run a small dataset
# in the large-data mode, which is how the authors' own script ran Zeisel.
# Small datasets use MICA's small-data mode with k set to the true number of
# classes. "ge" uses its large-data mode with the paper's resolution sweep.
set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate mica-legacy
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
name=$1; mode=$2; run=${3:+${name}_$3}; run=${run:-$name}
out=results/reproduction/mica_$run
log=logs/mica_$run.log
if [[ -d "$out" ]]; then   # old output would be pooled with the new run when scoring
  mv "$out" "${out}.stale_$(date +%Y%m%d_%H%M%S)"
fi
input=data/processed/${run}_mica.h5ad   # one input file per run, so parallel repeats never collide
python scripts/export_for_mica.py "$name" "$input" > "$log" 2>&1
start=$(date +%s)
if [[ $mode == ge ]]; then
  mica ge -i $input -o "$out" -minr 0.1 -maxr 4.0 -ss 0.1 -nw 8 >> "$log" 2>&1
else
  mica mds -i $input -o "$out" -nck "$mode" -nw 4 >> "$log" 2>&1
fi
echo "WALL_SECONDS $(( $(date +%s) - start ))" >> "$log"
rm -f $input   # the dense copy is only needed for the run
echo "done $run"
