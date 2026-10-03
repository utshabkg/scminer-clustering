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
import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from . import mi


class MICodes:
    """Cells binned once, with marginal entropies, ready for exact pair scoring."""

    def __init__(self, X, bins=None, rule="floor"):
        X = np.asarray(X)
        self.bins = bins or mi.n_bins(X.shape[1])
        self.codes = mi.bin_codes(X, self.bins, rule)
        self.h = mi.marginal_entropy(self.codes, self.bins)

    def dist(self, i, js):
        """Exact graph-mode distance from cell i to the cells js."""
        js = np.asarray(js)
        hij = mi.joint_entropy(self.codes[i], self.codes[js], self.bins)
        return 1 - (self.h[i] + self.h[js] - hij) / hij


def candidates_pca(emb, c):
    nn = NearestNeighbors(n_neighbors=c + 1).fit(emb)
    return nn.kneighbors(emb, return_distance=False)[:, 1:]


def candidates_detected(X, c):
    """Top c cells by Jaccard overlap of detected genes. Dense n x n, so small data only."""
    B = sparse.csr_matrix((np.asarray(X) > 0).astype(np.float32))
    inter = (B @ B.T).toarray()
    size = np.asarray(B.sum(axis=1)).ravel()
    jac = inter / (size[:, None] + size[None, :] - inter)
    np.fill_diagonal(jac, -1)
    return np.argsort(-jac, axis=1)[:, :c]


def union(*lists):
    """Per-cell union of several candidate lists. Rows can differ in length."""
    return [np.unique(np.concatenate(rows)) for rows in zip(*lists)]


def rerank(codes, cand, k):
    """Score every candidate exactly and keep the k closest. cand is an array or a list of rows."""
    n = len(cand)
    idx = np.empty((n, k), dtype=np.int64)
    dst = np.empty((n, k))
    for i in range(n):
        ci = np.asarray(cand[i])
        d = codes.dist(i, ci)
        o = np.argsort(d)[:k]
        idx[i], dst[i] = ci[o], d[o]
    return idx, dst


def n_scored(cand):
    return float(np.mean([len(c) for c in cand]))


def knn_graph(idx, dst):
    """Symmetric weighted graph from neighbor lists. Weight is the similarity I/H(x,y) = 1 - D."""
    n, k = idx.shape
    W = sparse.csr_matrix((np.clip(1 - dst, 0, None).ravel(), (np.repeat(np.arange(n), k), idx.ravel())), shape=(n, n))
    return W.maximum(W.T)


def leiden_on_graph(W, truth, seed=0):
    """Leiden on a given graph, swept and picked by the paper's rule.

    A stand-in for MICA's node2vec plus Louvain, used only to ask whether an
    approximate graph clusters the same way as the exact one.
    """
    import anndata as ad
    import scanpy as sc

    from .cluster import PAPER_RESOLUTIONS, pick_paper_protocol

    a = ad.AnnData(obs={"i": np.arange(W.shape[0]).astype(str)})
    a.obs_names = a.obs["i"]
    partitions = {}
    for res in PAPER_RESOLUTIONS:
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
