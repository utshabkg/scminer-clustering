"""Peak memory and time to build the correlation-shortlist graph, for the cost comparison.

Usage: python scripts/measure_graph_cost.py <dataset> <c> [workers] [dense]

The Assignment 2 runs measured recall and time for each graph but not memory, which the
feedback asked for. This repeats the all-cells pipeline of run_rq2_large.py (PBMC14K) and
run_rq2_allen.py (Allen) with memory recorded after each stage:

  1. load the matrix and bin every cell once (mi_knn.MICodes)
  2. shortlist c candidates per cell by correlation over all genes
  3. rescore each shortlist with the exact MICA distance and keep the K closest

Nothing guarantees the shortlist contains every true neighbour, so recall is measured, not
assumed.

Peak resident memory is read from the kernel and is a high-water mark, so the value after
loading is what loading alone cost, and the final value covers the whole build. The matrix
is dropped once the shortlist is built, as run_rq2_allen.py does, since rescoring needs only
the binned codes.

Recall is checked here too, against brute force cached in data/runs/rq2/: every cell where
the full exact graph exists (PBMC14K), otherwise the same 500 query cells as Assignment 2.

Passing "dense" forces the dense exact distance. The default sparse path skips genes where a
cell sits in its lowest bin and is exact by construction. Running both and comparing the
neighbour sets shows that directly, and accounts for the gap between today's timing and the
dense-path timing in Assignment 2.

The MICA side is scripts/measure_mica_hnsw.py, which has to run in the authors' environment.
Writes results/graph_cost_<dataset>_corr<c>[_dense]/.
"""
import resource
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse

from scminer_clust import metrics, mi_knn
from scminer_clust.paths import DATA, PROCESSED

K = 79  # MICA's 80 neighbours count the cell itself


def peak_gb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2, 2)


def recall_vs_brute_force(name, idx):
    """Recall on every cell if the full exact graph is cached, else on the Assignment 2 queries."""
    full = DATA / "runs" / "rq2" / f"{name}_exact_k{K}.npy"
    q500 = DATA / "runs" / "rq2" / f"{name}_exact_q500.npz"
    if full.exists():
        truth = np.load(full)
        q = np.sort(np.random.default_rng(0).choice(len(truth), 500, replace=False))
        return mi_knn.recall(idx, truth), mi_knn.recall(idx[q], truth[q]), "all cells"
    if q500.exists():
        z = np.load(q500)
        return None, mi_knn.recall(idx[z["q"]], z["idx"]), "500 query cells"
    return None, None, "no cached brute force"


def main(name, c, workers=8, dense=""):
    c, workers, dense = int(c), int(workers), dense == "dense"
    t0 = time.time()
    a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    n, g = X.shape
    del a
    codes = mi_knn.MICodes(X, sparse_codes=False if dense else None)
    t_load, m_load = time.time() - t0, peak_gb()
    print(f"{name}: {n} cells, {g} genes, loaded and binned in {t_load:.0f} s, peak {m_load} GB", flush=True)

    t0 = time.time()
    cand = mi_knn.candidates_correlation(X, c)
    t_cand, m_cand = time.time() - t0, peak_gb()
    del X
    print(f"  shortlist of {c} by correlation: {t_cand:.0f} s, peak {m_cand} GB", flush=True)

    t0 = time.time()
    idx, dst = mi_knn.rerank(codes, cand, K, workers=workers)
    t_rr, m_all = time.time() - t0, peak_gb()
    assert idx.shape == (n, K), idx.shape
    print(f"  exact rescoring ({'sparse' if codes.sparse else 'dense'} path): {t_rr:.0f} s, peak {m_all} GB",
          flush=True)
    r_all, r_q, r_on = recall_vs_brute_force(name, idx)
    print(f"  recall@{K}: {r_all if r_all is None else round(r_all, 4)} over all cells, "
          f"{r_q if r_q is None else round(r_q, 4)} on the 500 query cells ({r_on})", flush=True)
    tag = f"corr{c}" + ("_dense" if dense else "")
    out_idx = DATA / "runs" / "rq2" / f"graph_cost_{name}_{tag}_idx.npy"
    np.save(out_idx, idx)  # kept outside the repo, for comparing the two exact paths

    row = {"dataset": name, "graph": tag, "exact_path": "sparse" if codes.sparse else "dense",
           "c": c, "k": K, "workers": workers, "n_cells": n,
           "n_genes": g, "seconds_load_and_bin": round(t_load, 1), "seconds_candidates": round(t_cand, 1),
           "seconds_rescoring": round(t_rr, 1), "seconds_graph": round(t_cand + t_rr, 1),
           "peak_rss_gb_after_load": m_load, "peak_rss_gb_after_candidates": m_cand,
           "peak_rss_gb": m_all, "recall_all_cells": None if r_all is None else round(r_all, 4),
           "recall_500_queries": None if r_q is None else round(r_q, 4)}
    d = metrics.write(f"graph_cost_{name}_{tag}", [row],
                      {"dataset": name, "c": c, "k": K, "workers": workers,
                       "exact_path": "sparse" if codes.sparse else "dense",
                       "pipeline": "MICodes, then candidates_correlation, then rerank. X dropped before rescoring",
                       "memory": "resource.getrusage(RUSAGE_SELF).ru_maxrss, a high-water mark"})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
