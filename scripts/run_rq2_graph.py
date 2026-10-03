"""RQ2 prototype: how many exact comparisons does MICA's neighbor graph really need? (U-19)

Usage: python scripts/run_rq2_graph.py usoskin [k]
Builds the exact graph-mode distance for every pair once (small datasets only,
cached in data/runs/rq2/), then compares retrieve-then-rerank neighbor lists
against it:
  - recall of the k nearest neighbors and the share of all pairs scored
  - the share of cells whose list is certified exact by landmark bounds
  - whether Leiden on the approximate graph clusters like Leiden on the exact one
k defaults to 80, MICA's default neighbor count. Writes results/rq2_graph_<name>/.
"""
import json
import sys
import time

import anndata as ad
import numpy as np
from sklearn.metrics import adjusted_rand_score

from scminer_clust import cluster, metrics, mi, mi_knn
from scminer_clust.paths import DATA, PROCESSED

C_VALUES = (100, 200, 400)


def full_matrix(name, X):
    """Exact graph-mode distances, cached with the seconds they took to compute."""
    cache = DATA / "runs" / "rq2"; cache.mkdir(parents=True, exist_ok=True)
    f, meta = cache / f"{name}_full_ge.npy", cache / f"{name}_full_ge.json"
    if f.exists() and meta.exists():
        return np.load(f), json.loads(meta.read_text())["seconds"]
    t0 = time.time()
    D = mi.distance_matrix(X, "ge")
    secs = round(time.time() - t0, 1)
    np.save(f, D); meta.write_text(json.dumps({"seconds": secs, "shape": list(D.shape)}))
    return D, secs


def main(name, k=80):
    a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    X = np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X, dtype=np.float32)
    truth_labels = a.obs["label"].to_numpy()
    n = X.shape[0]
    D, t_full = full_matrix(name, X)
    np.fill_diagonal(D, np.inf)
    true_idx = np.argsort(D, axis=1)[:, :k]
    true_dst = np.take_along_axis(D, true_idx, axis=1)
    exact_best, exact_parts, exact_info = mi_knn.leiden_on_graph(mi_knn.knn_graph(true_idx, true_dst), truth_labels)
    exact_ari = adjusted_rand_score(truth_labels, exact_best)
    print(f"{name}: {n} cells, full matrix {t_full:.0f} s, Leiden on the exact {k}-neighbor MI graph ARI {exact_ari:.3f}")
    codes = mi_knn.MICodes(X)
    emb = cluster.pca(a, seed=0)
    rows = []
    for src in ("pca", "detected", "union"):
        for c in C_VALUES:
            if c >= n - 1:
                continue
            t1 = time.time()
            if src == "pca":
                cand = mi_knn.candidates_pca(emb, c)
            elif src == "detected":
                cand = mi_knn.candidates_detected(X, c)
            else:  # top c from each source, merged, so between c and 2c per cell
                cand = mi_knn.union(mi_knn.candidates_pca(emb, c), mi_knn.candidates_detected(X, c))
            idx, dst = mi_knn.rerank(codes, cand, k)
            t_search = time.time() - t1
            cert, n_l = mi_knn.landmark_certificate(codes, idx, dst, cand)
            best, parts, info = mi_knn.leiden_on_graph(mi_knn.knn_graph(idx, dst), truth_labels)
            # Same resolution on both graphs: how far apart are the partitions themselves?
            same_res = float(np.median([adjusted_rand_score(exact_parts[r], parts[r]) for r in (0.5, 1.0, 1.5)]))
            m = mi_knn.n_scored(cand)
            r = {"dataset": name, "candidates": src, "c": c, "mean_scored": round(m, 1), "k": k, "n_cells": n,
                 "recall_at_k": round(mi_knn.recall(idx, true_idx), 4), "share_of_pairs_scored": round(m / (n - 1), 4),
                 "certified_exact_share": round(cert, 4), "landmarks": n_l,
                 "leiden_ari_exact_graph": round(exact_ari, 4), "leiden_ari_this_graph": round(adjusted_rand_score(truth_labels, best), 4),
                 "partition_agreement_same_resolution": round(same_res, 4), "hit_true_k": info["hit_true_k"],
                 "seconds_search": round(t_search, 1), "seconds_full_matrix": t_full}
            rows.append(r)
            print(f"  {src:9s} c={c:4d} scored {m:6.0f}  recall@{k} {r['recall_at_k']:.3f}  pairs {r['share_of_pairs_scored']:.2f}"
                  f"  certified {cert:.3f}  Leiden ARI {r['leiden_ari_this_graph']:.3f}  agreement {same_res:.3f}  {t_search:.0f} s")
    d = metrics.write(f"rq2_graph_{name}", rows, {"dataset": name, "k": k, "distance": "graph mode 1 - I/H(x,y)",
                                                   "exact_graph_leiden": exact_info, "candidate_sizes": list(C_VALUES)})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 80))
