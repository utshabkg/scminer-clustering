#!/usr/bin/env bash
# Build the project conda envs. Usage: scripts/setup_envs.sh main|mica
set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

case "${1:-}" in
  main)
    conda env list | grep -q "^scminer " || conda env create -y -f "$ROOT/envs/scminer.yml"
    conda activate scminer
    # The default PyPI torch wheel is built for CUDA 13.0, which this driver
    # (CUDA 12.9) cannot load. Install the cu128 build instead.
    pip install --index-url https://download.pytorch.org/whl/cu128 --force-reinstall torch
    # The project itself, so scripts/ notebooks/ and tests/ can import scminer_clust.
    pip install -e "$ROOT"
    python -c "import scanpy, torch; print('scanpy', scanpy.__version__, 'torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.device_count())"
    conda env export -n scminer --no-builds > "$ROOT/envs/scminer.lock.yml"
    # conda records the editable install as `scminer-clust==0.1.0`, a name that exists on no
    # index, which makes the lock file impossible to build from. Put the editable install back.
    sed -i 's/^      - scminer-clust==.*$/      - -e ./' "$ROOT/envs/scminer.lock.yml"
    ;;
  mica)
    # MICA 1.0.0 from source, following the authors' README (Python 3.9.2).
    # pip must be inside the env: MICA's setup.py shells out to `python -m pip`
    # during the build, which fails otherwise.
    [[ -d "$ROOT/external/MICA" ]] || git clone https://github.com/jyyulab/MICA.git "$ROOT/external/MICA"
    conda env list | grep -q "^mica-legacy " || conda create -y -n mica-legacy python=3.9.2 pip
    conda activate mica-legacy
    cd "$ROOT/external/MICA"
    # MICA's setup.py runs `pip install ./mihnsw` while building. The bundled
    # mihnsw imports numpy and pybind11 in its setup.py but declares no build
    # requirements, so that nested isolated build fails. Give it a build config,
    # pinned to the numpy in this env so the extension stays ABI compatible.
    if [[ ! -f mihnsw/pyproject.toml ]]; then
      cat > mihnsw/pyproject.toml <<'TOML'
[build-system]
requires = ["setuptools>=65.6.3", "wheel", "pybind11", "numpy==1.22.4"]
build-backend = "setuptools.build_meta"
TOML
    fi
    python -m pip install -r requirements.txt
    python -m pip install --no-build-isolation .
    cd "$ROOT"
    mica -h | head -5
    python -c "import hnswlib; print('hnswlib (mihnsw) ok')"
    conda env export -n mica-legacy --no-builds > "$ROOT/envs/mica-legacy.lock.yml"
    ;;
  *)
    echo "usage: $0 main|mica" >&2; exit 1 ;;
esac
echo "done: $1"
