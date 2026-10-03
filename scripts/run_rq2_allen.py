"""RQ2 at 100,000 cells: correlation candidates, exact rescoring, against MICA's HNSW (U-19).

Usage: python scripts/run_rq2_allen.py [workers]
On the Allen pilot, MICA's HNSW step took 5.2 hours. This measures, on the same
cells and the same distance:
  1. the exact neighbors of 500 random cells by brute force (cached after the first run)
  2. how many of them MICA's saved graph found
  3. the same for correlation candidates of several sizes, rescored exactly
  4. full graphs for every cell at FULL_C candidates, timed
  5. Leiden on the largest one and on MICA's graph, on a coarse grid of 20 resolutions
The largest graph is also written in MICA's edge-list format, so MICA's own
node2vec and Louvain can run on it (scripts/mica_downstream.py).
Exact scoring uses the sparse path in mi_knn.MICodes. The first, dense run is in
logs/rq2_allen_dense_partial.log. Writes results/rq2_allen_100k/.
"""
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from scminer_clust import metrics, mi_knn
from scminer_clust.paths import DATA, PROCESSED

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from run_rq2_large import K, mica_graph  # noqa: E402

NAME = "allen_100k"
N_QUERIES = 500
CORR_C = (200, 400, 800, 1600, 3200)
FULL_C = (1600, 3200)
MICA_SECONDS = 18887.9  # "kNN/ANN costs" in logs/mica_allen_100k.log, 8 workers
DENSE_BRUTE_SECONDS = 1832  # the same 500 cells with the dense exact path, logs/rq2_allen_dense_partial.log
RESOLUTIONS = np.round(np.arange(0.1, 2.0 + 1e-9, 0.1), 1)


def main(workers=8):
    t_start = time.time()
    a = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    truth_labels = a.obs["label"].to_numpy()
    cell_names = a.obs_names.astype(str).tolist()
    n = X.shape[0]
    del a
    codes = mi_knn.MICodes(X)
    q = np.sort(np.random.default_rng(0).choice(n, N_QUERIES, replace=False))
    print(f"{NAME}: {n} cells, {X.shape[1]} genes, {np.mean(codes.codes > 0):.3f} of codes nonzero, "
          f"loaded and binned in {time.time() - t_start:.0f} s")

    m_idx, _, edges, W_mica = mica_graph(NAME, n)
    s = edges.sample(2000, random_state=0)
    ours = np.array([codes.dist(int(i), np.array([int(j)]))[0] for i, j in zip(s.i, s.j)])
    max_diff = float(np.abs(ours - s.d.to_numpy()).max())
    print(f"  our distance against MICA's saved weights: max difference {max_diff:.2e}")
    del edges

    cache = DATA / "runs" / "rq2" / f"{NAME}_exact_q{N_QUERIES}.npz"
    t0 = time.time()
    true_idx, _ = mi_knn.exact_knn(codes, q, K, workers)
    t_exact = time.time() - t0
    if cache.exists():
        prev = np.load(cache)
        assert np.array_equal(prev["q"], q) and mi_knn.recall(true_idx, prev["idx"]) > 0.999, "brute force changed"
    np.savez(cache, q=q, idx=true_idx)
    est_all = t_exact * n / N_QUERIES
    mica_recall = mi_knn.recall(m_idx[q], true_idx)
    print(f"  brute force for {N_QUERIES} cells: {t_exact:.0f} s (dense path took {DENSE_BRUTE_SECONDS} s), "
          f"so about {est_all / 3600:.1f} h for all cells")
    print(f"  MICA HNSW: recall@{K} {mica_recall:.3f} in {MICA_SECONDS / 3600:.1f} h")
    rows = [{"dataset": NAME, "stage": "query sweep", "candidates": "brute force (estimated)", "k": K, "recall_at_k": 1.0,
             "share_of_pairs_scored": 1.0, "seconds_all_cells": round(est_all, 1), "workers": workers,
             "seconds_all_cells_dense_path": round(DENSE_BRUTE_SECONDS * n / N_QUERIES, 1)},
            {"dataset": NAME, "stage": "query sweep", "candidates": "MICA HNSW, K closest links", "k": K,
             "recall_at_k": round(mica_recall, 4), "seconds_all_cells": MICA_SECONDS, "workers": 8}]

    cq = mi_knn.candidates_correlation(X, max(CORR_C), rows=q)
    for c in CORR_C:
        idx, _ = mi_knn.rerank(codes, cq[:, :c], K, rows=q, workers=workers)
        r = mi_knn.recall(idx, true_idx)
        rows.append({"dataset": NAME, "stage": "query sweep", "candidates": "correlation", "c": c, "mean_scored": c, "k": K,
                     "recall_at_k": round(r, 4), "share_of_pairs_scored": round(c / (n - 1), 4)})
        print(f"  correlation c={c:5d}: recall@{K} {r:.3f}, {c / (n - 1):.2%} of pairs")

    t0 = time.time()
    cand = mi_knn.candidates_correlation(X, max(FULL_C))
    t_cand = time.time() - t0
    del X
    print(f"  correlation candidates for all cells: {t_cand / 60:.0f} min")
    graphs = {}
    for c in FULL_C:
        t0 = time.time()
        idx, dst = mi_knn.rerank(codes, cand[:, :c], K, workers=workers)
        t_rr = time.time() - t0
        r = mi_knn.recall(idx[q], true_idx)
        graphs[c] = (idx, dst)
        rows.append({"dataset": NAME, "stage": "all cells, round 0", "candidates": "correlation", "c": c, "mean_scored": c,
                     "k": K, "recall_at_k": round(r, 4), "share_of_pairs_scored": round(c / (n - 1), 4),
                     "seconds_all_cells": round(t_cand + t_rr, 1), "seconds_candidates": round(t_cand, 1),
                     "seconds_rescoring": round(t_rr, 1), "workers": workers})
        print(f"  all cells, correlation c={c}: recall@{K} {r:.3f} in {(t_cand + t_rr) / 3600:.2f} h "
              f"({t_cand / 60:.0f} min candidates, {t_rr / 60:.0f} min exact rescoring), against MICA's {MICA_SECONDS / 3600:.1f} h")
    del cand

    big = max(FULL_C)
    idx, dst = graphs[big]
    n_edges = mi_knn.write_mica_edgelist(idx, dst, DATA / "runs" / "rq2" / f"{NAME}_corr{big}", cell_names)
    print(f"  wrote the c={big} graph in MICA's format, {n_edges} edges")
    fast_best, fast_parts, fast_info = mi_knn.leiden_on_graph(mi_knn.knn_graph(idx, dst), truth_labels, resolutions=RESOLUTIONS)
    mica_best, mica_parts, mica_info = mi_knn.leiden_on_graph(W_mica, truth_labels, resolutions=RESOLUTIONS)
    agree = float(np.median([adjusted_rand_score(mica_parts[x], fast_parts[x]) for x in (0.5, 1.0, 1.5)]))
    rows.append({"dataset": NAME, "stage": "clustering", "candidates": f"correlation {big}", "k": K,
                 "leiden_ari_fast_graph": round(adjusted_rand_score(truth_labels, fast_best), 4),
                 "leiden_ari_mica_graph": round(adjusted_rand_score(truth_labels, mica_best), 4),
                 "partition_agreement_same_resolution": round(agree, 4),
                 "note": "class level, Leiden on 20 resolutions stands in for MICA's node2vec and Louvain"})
    print(f"  Leiden ARI (class level): fast graph {rows[-1]['leiden_ari_fast_graph']:.3f}, "
          f"MICA graph {rows[-1]['leiden_ari_mica_graph']:.3f}, agreement {agree:.3f}")
    d = metrics.write(f"rq2_{NAME}", rows, {"dataset": NAME, "k": K, "queries": N_QUERIES, "query_seed": 0,
                                             "distance": "graph mode 1 - I/H(x,y)", "check_max_diff_vs_mica": max_diff,
                                             "exact_path": "sparse" if codes.sparse else "dense",
                                             "resolutions": RESOLUTIONS.tolist(), "leiden_fast": fast_info, "leiden_mica": mica_info,
                                             "total_seconds": round(time.time() - t_start, 1)})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
