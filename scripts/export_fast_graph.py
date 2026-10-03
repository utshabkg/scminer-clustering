"""Write a fast MI neighbor graph in MICA's edge-list format, so MICA's own later steps can use it (U-19).

Usage: python scripts/export_fast_graph.py paper_pbmc14k 800 [workers]
Builds each cell's 79 nearest neighbors from its top-c correlation candidates,
rescored with the exact graph-mode distance, and writes the undirected graph the
way MICA does: one line per pair, "cell<TAB>neighbor<TAB>distance". Output goes to
data/runs/rq2/<name>_corr<c>/ with edgelist.txt and cells.txt (cell IDs in order).
"""
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse

from scminer_clust import mi_knn
from scminer_clust.paths import DATA, PROCESSED

K = 79  # MICA's 80 neighbors include the cell itself


def main(name, c, workers=8):
    a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    t0 = time.time()
    codes = mi_knn.MICodes(X)
    idx, dst = mi_knn.rerank(codes, mi_knn.candidates_correlation(X, c), K, workers=workers)
    secs = time.time() - t0
    out = DATA / "runs" / "rq2" / f"{name}_corr{c}"
    n_edges = mi_knn.write_mica_edgelist(idx, dst, out, a.obs_names.astype(str))
    (out / "seconds.txt").write_text(f"{secs:.1f}\n")
    print(f"{name}: {n_edges} edges for {len(idx)} cells in {secs:.0f} s, wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 8))
