#!/usr/bin/env bash
# Build the project conda envs. Usage: scripts/setup_envs.sh main|mica
set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

case "${1:-}" in
  main)
    conda env create -y -f "$ROOT/envs/scminer.yml"
    conda activate scminer
    python -c "import scanpy, torch; print('scanpy', scanpy.__version__, 'torch', torch.__version__, 'cuda', torch.cuda.is_available())"
    conda env export -n scminer --no-builds > "$ROOT/envs/scminer.lock.yml"
    ;;
  mica)
    # MICA 1.0.0 from source, following the authors' README (Python 3.9.2)
    [[ -d "$ROOT/external/MICA" ]] || git clone https://github.com/jyyulab/MICA.git "$ROOT/external/MICA"
    conda create -y -n mica-legacy python=3.9.2
    conda activate mica-legacy
    cd "$ROOT/external/MICA"
    pip install .
    cd "$ROOT"
    mica -h | head -5
    python -c "import mihnsw; print('mihnsw ok')"
    conda env export -n mica-legacy --no-builds > "$ROOT/envs/mica-legacy.lock.yml"
    ;;
  *)
    echo "usage: $0 main|mica" >&2; exit 1 ;;
esac
echo "done: $1"
