"""Retrieve then rerank: MICA's mutual information neighbor graph without all pairs (RQ2, U-19).

MICA's graph mode needs each cell's k nearest neighbors under the distance
D = 1 - I(x, y) / H(x, y). Exact search compares every pair, which is quadratic.
Here, as in a search engine, a cheap first stage proposes candidates and only those
are scored exactly:
  - principal components: nearest neighbors on 50 PCs of the log data
  - detected genes: Jaccard overlap of the sets of genes detected in each cell
    (co-detection carries much of MICA's similarity, see diagnostics.py)
Because D is a metric (Vinh et al. 2010), distances to a few landmark cells give
a lower bound |D(i, l) - D(j, l)| <= D(i, j) for any pair. A cell's list is
certified exact when every non-candidate's lower bound already exceeds the
distance of its k-th neighbor. This CPU version is for small datasets; the GPU
version (U-18) will reuse the same logic.
"""
import multiprocessing as mp

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from . import mi

BLOCK = 512  # cells per joint-histogram block, keeps memory near BLOCK x genes x 8 bytes


class MICodes:
    """Cells binned once, with marginal entropies, ready for exact pair scoring.

    With sparse=True (the default when under half the codes are nonzero), the
    joint histogram of cells i and j is built only from the genes where j's code
    is nonzero. The column for j's code 0 then follows from i's own bin counts,
    so the result is exact, and the cost scales with j's nonzero genes instead of
    all genes. On 10x data about 6 percent of codes are nonzero.
    """

    def __init__(self, X, bins=None, rule="floor", sparse_codes=None):
        X = np.asarray(X)
        self.bins = bins or mi.n_bins(X.shape[1])
        self.codes = mi.bin_codes(X, self.bins, rule)
        self.h = mi.marginal_entropy(self.codes, self.bins)
        self.n_genes = self.codes.shape[1]
        if sparse_codes is None:
            sparse_codes = (self.codes > 0).mean() < 0.5
        self.sparse = bool(sparse_codes)
        if self.sparse:
            self.csr = sparse.csr_matrix(self.codes)
            self.csr.eliminate_zeros()
            self.counts = np.stack([np.bincount(r, minlength=self.bins) for r in self.codes]).astype(np.int64)

    def _joint_entropy_sparse(self, i, js):
        b = self.bins
        ci = self.codes[i].astype(np.int64)
        sub = self.csr[js]
        rows = np.repeat(np.arange(len(js)), np.diff(sub.indptr))
        comb = rows * b * b + ci[sub.indices] * b + sub.data.astype(np.int64)
        A = np.bincount(comb, minlength=len(js) * b * b).reshape(len(js), b, b)
        A[:, :, 0] = self.counts[i][None, :] - A[:, :, 1:].sum(axis=2)  # genes where j's code is 0
        return mi._entropy_from_counts(A.reshape(len(js), b * b), self.n_genes)

    def dist(self, i, js):
        """Exact graph-mode distance from cell i to the cells js."""
        js = np.asarray(js)
        if not len(js):
            return np.empty(0)
        if self.sparse:
            hij = np.concatenate([self._joint_entropy_sparse(i, js[s:s + 4 * BLOCK]) for s in range(0, len(js), 4 * BLOCK)])
        else:
            hij = np.concatenate([mi.joint_entropy(self.codes[i], self.codes[js[s:s + BLOCK]], self.bins)
                                  for s in range(0, len(js), BLOCK)])
        return 1 - (self.h[i] + self.h[js] - hij) / hij


def candidates_pca(emb, c, rows=None):
    """The c nearest cells on principal components, for every cell or only `rows`."""
    q = np.arange(len(emb)) if rows is None else np.asarray(rows)
    nn = NearestNeighbors(n_neighbors=c + 1).fit(emb)
    idx = nn.kneighbors(emb[q], return_distance=False)
    return np.stack([r[r != i][:c] for r, i in zip(idx, q)])


def candidates_detected(X, c, rows=None, block=1000):
    """Top c cells by Jaccard overlap of detected genes, a block of query cells at a time."""
    B = sparse.csr_matrix((np.asarray(X) > 0).astype(np.float32))
    size = np.asarray(B.sum(axis=1)).ravel()
    q = np.arange(B.shape[0]) if rows is None else np.asarray(rows)
    out = np.empty((len(q), c), dtype=np.int64)
    for s in range(0, len(q), block):
        r = q[s:s + block]
        inter = (B[r] @ B.T).toarray()
        jac = inter / (size[r][:, None] + size[None, :] - inter)
        jac[np.arange(len(r)), r] = -1
        top = np.argpartition(-jac, c, axis=1)[:, :c]
        order = np.argsort(-np.take_along_axis(jac, top, axis=1), axis=1)
        out[s:s + block] = np.take_along_axis(top, order, axis=1)
    return out


def candidates_correlation(X, c, rows=None, block=2000):
    """Top c cells by Pearson correlation over all genes, a block of query cells at a time.

    No dimension reduction: on PBMC14K, 50 principal components lose most of the
    neighbor structure MICA sees, while plain correlation keeps much of it.
    """
    Z = np.asarray(X, dtype=np.float32)
    Z = Z - Z.mean(axis=1, keepdims=True)
    Z /= np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12
    q = np.arange(Z.shape[0]) if rows is None else np.asarray(rows)
    out = np.empty((len(q), c), dtype=np.int64)
    for s in range(0, len(q), block):
        r = q[s:s + block]
        sim = Z[r] @ Z.T
        sim[np.arange(len(r)), r] = -np.inf
        top = np.argpartition(-sim, c, axis=1)[:, :c]
        order = np.argsort(-np.take_along_axis(sim, top, axis=1), axis=1)
        out[s:s + block] = np.take_along_axis(top, order, axis=1)
    return out


def refine(codes, idx, k, fan=20, workers=1):
    """One round of neighbors of neighbors (as in NN-descent), rescored exactly.

    Each cell's new candidates are its current list plus the first `fan` neighbors
    of each of its first `fan` neighbors. Returns the new lists and how many were scored.
    """
    n = idx.shape[0]
    cand = []
    for i in range(n):
        two = idx[idx[i, :fan], :fan].ravel()
        cand.append(np.unique(np.concatenate([idx[i], two[two != i]])))
    new_idx, new_dst = rerank(codes, cand, k, workers=workers)
    return new_idx, new_dst, n_scored(cand)


def write_mica_edgelist(idx, dst, out_dir, cell_names):
    """Write neighbor lists as MICA writes its graph: each undirected pair once, cell, neighbor, distance."""
    from pathlib import Path

    n, k = idx.shape
    i = np.repeat(np.arange(n), k)
    j, d = idx.ravel(), dst.ravel()
    lo, hi = np.minimum(i, j), np.maximum(i, j)
    _, keep = np.unique(lo * n + hi, return_index=True)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "edgelist.txt", "w") as f:
        for a, b, w in zip(i[keep], j[keep], d[keep]):
            f.write(f"{a}\t{b}\t{w}\n")
    (out / "cells.txt").write_text("\n".join(map(str, cell_names)) + "\n")
    return len(keep)


def union(*lists):
    """Per-cell union of several candidate lists. Rows can differ in length."""
    return [np.unique(np.concatenate(rows)) for rows in zip(*lists)]


_SHARED = {}  # set before forking, so worker processes share the codes without copying


def _top_k(task):
    i, js, k = task
    d = _SHARED["codes"].dist(i, js)
    o = np.argsort(d)[:k]
    return js[o], d[o]


def _run(codes, tasks, workers):
    if workers <= 1:
        _SHARED["codes"] = codes
        return [_top_k(t) for t in tasks]
    _SHARED["codes"] = codes
    with mp.get_context("fork").Pool(workers) as pool:
        return pool.map(_top_k, tasks, chunksize=max(1, len(tasks) // (workers * 8)))


def rerank(codes, cand, k, rows=None, workers=1):
    """Score every candidate exactly and keep the k closest.

    cand is an array or a list of rows, one per cell in `rows` (all cells if None).
    """
    q = np.arange(len(cand)) if rows is None else np.asarray(rows)
    res = _run(codes, [(int(i), np.asarray(c), k) for i, c in zip(q, cand)], workers)
    return np.stack([r[0] for r in res]), np.stack([r[1] for r in res])


def exact_knn(codes, rows, k, workers=1):
    """Brute force: each query cell against every other cell."""
    n = codes.codes.shape[0]
    allc = np.arange(n)
    res = _run(codes, [(int(i), allc[allc != i], k) for i in rows], workers)
    return np.stack([r[0] for r in res]), np.stack([r[1] for r in res])


def n_scored(cand):
    return float(np.mean([len(c) for c in cand]))


def knn_graph(idx, dst):
    """Symmetric weighted graph from neighbor lists. Weight is the similarity I/H(x,y) = 1 - D."""
    n, k = idx.shape
    W = sparse.csr_matrix((np.clip(1 - dst, 0, None).ravel(), (np.repeat(np.arange(n), k), idx.ravel())), shape=(n, n))
    return W.maximum(W.T)


def leiden_on_graph(W, truth, seed=0, resolutions=None):
    """Leiden on a given graph, swept and picked by the paper's rule (or over `resolutions`).

    A stand-in for MICA's node2vec plus Louvain, used only to ask whether an
    approximate graph clusters the same way as the exact one.
    """
    import anndata as ad
    import scanpy as sc

    from .cluster import PAPER_RESOLUTIONS, pick_paper_protocol

    a = ad.AnnData(obs={"i": np.arange(W.shape[0]).astype(str)})
    a.obs_names = a.obs["i"]
    partitions = {}
    for res in (PAPER_RESOLUTIONS if resolutions is None else resolutions):
        sc.tl.leiden(a, resolution=float(res), adjacency=W, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        partitions[float(res)] = a.obs["leiden"].to_numpy().astype(str)
    res, info = pick_paper_protocol(partitions, truth)
    return partitions[res], partitions, info


def recall(found, truth):
    k = truth.shape[1]
    return float(np.mean([len(set(a) & set(b)) / k for a, b in zip(found, truth)]))


def landmark_certificate(codes, idx, dst, cand, n_landmarks=64, seed=0):
    """Share of cells whose k-neighbor list is provably exact, via landmark lower bounds.

    Landmarks are picked by farthest-first traversal, which spreads them out. A
    non-candidate j can be ruled out for cell i when max over landmarks of
    |D(i, l) - D(j, l)| is larger than the distance of i's k-th neighbor.
    """
    n = codes.codes.shape[0]
    rng = np.random.default_rng(seed)
    L = [int(rng.integers(n))]
    dl = codes.dist(L[0], np.arange(n))[None, :]
    near = dl[0].copy()
    while len(L) < n_landmarks:
        nxt = int(np.argmax(near))
        L.append(nxt)
        row = codes.dist(nxt, np.arange(n))
        dl = np.vstack([dl, row])
        near = np.minimum(near, row)
    certified = 0
    for i in range(n):
        lb = np.abs(dl[:, i][:, None] - dl).max(axis=0)  # lower bound to every cell
        mask = np.ones(n, dtype=bool); mask[cand[i]] = False; mask[i] = False
        if (lb[mask] > dst[i, -1]).all():
            certified += 1
    return certified / n, len(L)
