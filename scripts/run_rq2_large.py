"""RQ2 prototype on a dataset too big for the full distance matrix (U-19).

Usage: python scripts/run_rq2_large.py paper_pbmc14k [workers]
The exact answer comes from brute force on a random sample of query cells, each
compared with every other cell. Against it we measure:
  - MICA's own approximate graph (HNSW, saved by its graph-mode run)
  - cheap candidates rescored exactly: principal components, shared detected
    genes, correlation over all genes, and merged lists, at several sizes
  - full pipelines for every cell, timed: candidates, exact rescoring, and for
    some, rounds of neighbors of neighbors (as in NN-descent)
Leiden on the correlation graph with 800 candidates is compared with Leiden on
MICA's graph.
Writes results/rq2_large_<name>/.
"""
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from scminer_clust import cluster, metrics, mi_knn
from scminer_clust.paths import PROCESSED, RESULTS

K = 79  # MICA's default of 80 neighbors counts the cell itself
N_QUERIES = 500
C_VALUES = (100, 200, 400)
CORR_C_VALUES = (100, 200, 400, 800, 1600)
# (candidates, per-cell count, refinement rounds) for the timed runs on every cell
PIPELINES = [("pca + detected", 200, 2), ("correlation", 200, 2), ("correlation", 400, 0),
             ("correlation", 800, 0), ("correlation", 1600, 0)]
CLUSTER_ON = ("correlation", 800)
MICA_SECONDS = {"paper_pbmc14k": 903.2}  # HNSW step, 8 workers, from logs/mica_<name>.log


def mica_graph(name, n):
    """MICA's saved graph, read back as each cell's K closest linked cells.

    MICA writes an undirected networkx graph, each pair once, so a cell's own list
    is mixed with the cells that picked it. Taking the K closest links per cell can
    only flatter MICA. Also returns the symmetric weighted graph MICA clusters on.
    """
    e = pd.read_csv(RESULTS / "reproduction" / f"mica_{name}" / "NNDescent_knn_graph_edgelist.txt",
                    sep="\t", header=None, names=["i", "j", "d"])
    both = pd.concat([e, e.rename(columns={"i": "j", "j": "i"})]).drop_duplicates(["i", "j"])
    both = both[both.i != both.j].sort_values(["i", "d"])
    deg = both.groupby("i").size()
    assert len(deg) == n and deg.min() >= K, f"MICA graph has {len(deg)} nodes, min degree {deg.min()}"
    top = both.groupby("i").head(K)
    W = sparse.csr_matrix((1 - both.d.to_numpy(), (both.i.to_numpy(), both.j.to_numpy())), shape=(n, n))
    return top.j.to_numpy().reshape(n, K), top.d.to_numpy().reshape(n, K), e, W


def main(name, workers=8):
    a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    truth_labels = a.obs["label"].to_numpy()
    n = X.shape[0]
    codes = mi_knn.MICodes(X)
    q = np.sort(np.random.default_rng(0).choice(n, N_QUERIES, replace=False))

    # MICA's own graph, and a check that our distance matches the weights it saved
    m_idx, m_dst, edges, W_mica = mica_graph(name, n)
    s = edges.sample(2000, random_state=0)
    ours = np.array([codes.dist(int(i), np.array([int(j)]))[0] for i, j in zip(s.i, s.j)])
    max_diff = float(np.abs(ours - s.d.to_numpy()).max())
    print(f"{name}: {n} cells, {X.shape[1]} genes. Our distance against MICA's saved weights: max difference {max_diff:.2e}")

    t0 = time.time()
    true_idx, _ = mi_knn.exact_knn(codes, q, K, workers)
    t_exact = time.time() - t0
    est_all = t_exact * n / N_QUERIES
    print(f"  brute force for {N_QUERIES} cells: {t_exact:.0f} s, so about {est_all / 60:.0f} min for all cells")
    rows = [{"dataset": name, "stage": "query sweep", "candidates": "brute force (estimated)", "k": K,
             "recall_at_k": 1.0, "share_of_pairs_scored": 1.0, "seconds_all_cells": round(est_all, 1), "workers": workers},
            {"dataset": name, "stage": "query sweep", "candidates": "MICA HNSW, K closest links", "k": K,
             "recall_at_k": round(mi_knn.recall(m_idx[q], true_idx), 4),
             "seconds_all_cells": MICA_SECONDS.get(name), "workers": 8}]
    print(f"  MICA HNSW                recall@{K} {rows[-1]['recall_at_k']:.3f}  {MICA_SECONDS.get(name)} s")

    emb = cluster.pca(a, seed=0)
    sources = {
        "pca": lambda c, r: mi_knn.candidates_pca(emb, c, rows=r),
        "detected": lambda c, r: mi_knn.candidates_detected(X, c, rows=r),
        "correlation": lambda c, r: mi_knn.candidates_correlation(X, c, rows=r),
        "pca + detected": lambda c, r: mi_knn.union(mi_knn.candidates_pca(emb, c, rows=r), mi_knn.candidates_detected(X, c, rows=r)),
        "correlation + detected": lambda c, r: mi_knn.union(mi_knn.candidates_correlation(X, c, rows=r), mi_knn.candidates_detected(X, c, rows=r)),
    }
    for src, make in sources.items():
        for c in (CORR_C_VALUES if src == "correlation" else C_VALUES):
            cand = make(c, q)
            idx, _ = mi_knn.rerank(codes, cand, K, rows=q, workers=workers)
            m = mi_knn.n_scored(cand)
            rows.append({"dataset": name, "stage": "query sweep", "candidates": src, "c": c, "mean_scored": round(m, 1), "k": K,
                         "recall_at_k": round(mi_knn.recall(idx, true_idx), 4), "share_of_pairs_scored": round(m / (n - 1), 4)})
            print(f"  {src:24s} c={c:4d} scored {m:6.0f}  recall@{K} {rows[-1]['recall_at_k']:.3f}")

    # Full pipelines for every cell: candidates, exact rescoring, refinement rounds
    graphs = {}
    for src, full_c, rounds in PIPELINES:
        t0 = time.time()
        cand = sources[src](full_c, None)
        t_cand = time.time() - t0
        idx, dst = mi_knn.rerank(codes, cand, K, workers=workers)
        scored = mi_knn.n_scored(cand)
        for rnd in range(rounds + 1):
            if rnd:
                idx, dst, extra = mi_knn.refine(codes, idx, K, workers=workers)
                scored += extra
            secs = time.time() - t0
            rows.append({"dataset": name, "stage": f"all cells, round {rnd}", "candidates": src, "c": full_c,
                         "mean_scored": round(scored, 1), "k": K, "recall_at_k": round(mi_knn.recall(idx[q], true_idx), 4),
                         "share_of_pairs_scored": round(scored / (n - 1), 4), "seconds_all_cells": round(secs, 1),
                         "seconds_candidates": round(t_cand, 1), "workers": workers})
            print(f"  all cells, {src:24s} c={full_c:4d} round {rnd}: recall@{K} {rows[-1]['recall_at_k']:.3f}, "
                  f"{scored:.0f} scored per cell, {secs:.0f} s")
        graphs[(src, full_c)] = (idx, dst, rows[-1]["recall_at_k"])

    best = f"{CLUSTER_ON[0]} {CLUSTER_ON[1]}"
    idx, dst, _ = graphs[CLUSTER_ON]
    fast_best, fast_parts, fast_info = mi_knn.leiden_on_graph(mi_knn.knn_graph(idx, dst), truth_labels)
    mica_best, mica_parts, mica_info = mi_knn.leiden_on_graph(W_mica, truth_labels)
    agree = float(np.median([adjusted_rand_score(mica_parts[r], fast_parts[r]) for r in (0.5, 1.0, 1.5)]))
    rows.append({"dataset": name, "stage": "clustering", "candidates": best, "k": K,
                 "leiden_ari_fast_graph": round(adjusted_rand_score(truth_labels, fast_best), 4),
                 "leiden_ari_mica_graph": round(adjusted_rand_score(truth_labels, mica_best), 4),
                 "partition_agreement_same_resolution": round(agree, 4),
                 "note": "Leiden stands in for MICA's node2vec and Louvain"})
    print(f"  Leiden ARI: fast graph ({best}) {rows[-1]['leiden_ari_fast_graph']:.3f}, "
          f"MICA graph {rows[-1]['leiden_ari_mica_graph']:.3f}, agreement {agree:.3f}")
    d = metrics.write(f"rq2_large_{name}", rows, {"dataset": name, "k": K, "queries": N_QUERIES, "query_seed": 0,
                                                    "distance": "graph mode 1 - I/H(x,y)", "check_max_diff_vs_mica": max_diff,
                                                    "refine": "neighbors of neighbors, first 20 of the first 20",
                                                    "pipelines": [list(p) for p in PIPELINES],
                                                    "leiden_fast": fast_info, "leiden_mica": mica_info})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 8))
