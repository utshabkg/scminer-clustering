"""Hub cells in the exact mutual information neighbor graph (U-19 side result).

Usage: python scripts/rq2_hubs.py [workers]
PBMC14K: builds (or loads) each cell's exact 79 nearest neighbors under MICA's
graph-mode distance, cached in data/runs/rq2/paper_pbmc14k_exact_k79.npy. A
cached graph is used only after 50 random rows match a fresh brute force. Then it
reports how much of the graph the top 1% of cells by in-degree take, and how
often their links join cells of the same type. Klein uses the cached full
distance matrix from run_rq2_graph.py. Writes results/rq2_hubs/.
"""
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse

from scminer_clust import metrics, mi_knn
from scminer_clust.paths import DATA, PROCESSED

K = 79


def hub_stats(name, idx, labels):
    n = len(idx)
    indeg = np.bincount(idx.ravel(), minlength=n)
    hubs = np.argsort(-indeg)[: n // 100]
    is_hub = np.zeros(n, bool)
    is_hub[hubs] = True
    same = labels[idx] == labels[:, None]
    slot_hub = is_hub[idx]
    return {"dataset": name, "n_cells": n, "k": idx.shape[1], "n_hubs": len(hubs),
            "hub_share_of_slots": round(float(slot_hub.mean()), 4),
            "same_type_share_hub_slots": round(float(same[slot_hub].mean()), 4),
            "same_type_share_other_slots": round(float(same[~slot_hub].mean()), 4),
            "indegree_median": float(np.median(indeg)), "indegree_p99": float(np.percentile(indeg, 99)),
            "indegree_max": int(indeg.max())}


def main(workers=8):
    rows = []
    a = ad.read_h5ad(PROCESSED / "paper_pbmc14k.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    labels = a.obs["label"].to_numpy()
    n = X.shape[0]
    codes = mi_knn.MICodes(X)
    cache = DATA / "runs" / "rq2" / "paper_pbmc14k_exact_k79.npy"
    if cache.exists():
        idx = np.load(cache)
        check = np.sort(np.random.default_rng(1).choice(n, 50, replace=False))
        fresh, _ = mi_knn.exact_knn(codes, check, K, workers)
        assert mi_knn.recall(idx[check], fresh) == 1.0, "cached exact graph does not match brute force"
        print(f"cached exact graph checked on 50 random cells")
    else:
        t0 = time.time()
        idx, _ = mi_knn.exact_knn(codes, np.arange(n), K, workers)
        np.save(cache, idx)
        print(f"exact graph for all {n} cells in {time.time() - t0:.0f} s")
    rows.append(hub_stats("paper_pbmc14k", idx, labels))

    k = ad.read_h5ad(PROCESSED / "klein.h5ad", backed="r")
    D = np.load(DATA / "runs" / "rq2" / "klein_full_ge.npy")
    np.fill_diagonal(D, np.inf)
    rows.append(hub_stats("klein", np.argsort(D, axis=1)[:, :80], k.obs["label"].to_numpy()))
    for r in rows:
        print(r)
    d = metrics.write("rq2_hubs", rows, {"top": "1% of cells by in-degree", "pbmc14k_k": K, "klein_k": 80,
                                          "distance": "graph mode 1 - I/H(x,y)"})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
