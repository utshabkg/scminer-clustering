"""Reference distances from the authors' own MICA functions. Run in mica-legacy.

Usage: conda run -n mica-legacy python scripts/make_mi_reference.py zeisel [n_pairs]
Writes data/processed/mi_reference_<name>.npz, which tests/test_mi_matches_mica.py
compares against our reimplementation in src/scminer_clust/mi.py.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from MICA.lib import distance, mutual_info

PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"


def main(name, n_pairs=10000):
    X = np.asarray(sparse.load_npz(PROCESSED / f"{name}_exchange" / "X.npz").todense(), dtype=np.float32)
    n, g = X.shape
    bins = int(g ** (1 / 3.0))  # MICA's default, as in mica_mds.py
    rng = np.random.default_rng(0)
    pairs = rng.integers(0, n, size=(n_pairs, 2))
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    d_mds = np.empty(len(pairs)); d_ge = np.empty(len(pairs))
    for k, (i, j) in enumerate(pairs):
        d_mds[k] = mutual_info.mi_norm(mutual_info.mi_calc(X[[i, j]], bins))[0, 1]   # small-data mode
        d_ge[k] = distance.numba_calc_mi_dis(X[i], X[j], bins, g)                     # graph mode, numba
    np.savez(PROCESSED / f"mi_reference_{name}.npz", pairs=pairs, d_mds=d_mds, d_ge=d_ge, bins=bins)
    print(f"{name}: {len(pairs)} pairs, {g} genes, {bins} bins")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10000))
