"""Baseline clustering methods, all on the same preprocessed data."""
import numpy as np
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering, KMeans


def pca(adata, n_comps=50, seed=0):
    """Principal components on the log normalized matrix, cached in .obsm."""
    key = f"X_pca_{n_comps}_{seed}"
    if key not in adata.obsm:
        a = adata.copy()
        sc.pp.scale(a, max_value=10)
        sc.tl.pca(a, n_comps=n_comps, svd_solver="arpack", random_state=seed)
        adata.obsm[key] = a.obsm["X_pca"]
    return adata.obsm[key]


def kmeans_pca(adata, k, n_comps=50, seed=0):
    emb = pca(adata, n_comps, seed)
    labels = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(emb)
    return labels.astype(str), emb


def hierarchical_pca(adata, k, n_comps=50, seed=0):
    emb = pca(adata, n_comps, seed)
    labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(emb)
    return labels.astype(str), emb


def leiden_pca(adata, resolution=1.0, n_comps=50, n_neighbors=15, seed=0):
    """The standard Scanpy pipeline: principal components, neighbor graph, Leiden."""
    a = adata.copy()
    a.obsm["X_pca"] = pca(adata, n_comps, seed)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return a.obs["leiden"].to_numpy().astype(str), a.obsm["X_pca"]


def leiden_at_k(adata, k, n_comps=50, n_neighbors=15, seed=0, lo=0.05, hi=4.0, steps=40):
    """Sweep the resolution until the cluster count matches k.

    This is the setting the scMINER paper used for its benchmark, where the true
    number of cell types decides the parameter. We keep it so we can reproduce
    their protocol, and compare it against choosing the number without labels.
    """
    best = None
    for res in np.linspace(lo, hi, steps):
        labels, emb = leiden_pca(adata, resolution=float(res), n_comps=n_comps, n_neighbors=n_neighbors, seed=seed)
        n = len(np.unique(labels))
        if best is None or abs(n - k) < abs(best[2] - k):
            best = (labels, emb, n, float(res))
        if n == k:
            break
    labels, emb, n, res = best
    return labels, emb, {"resolution": res, "n_clusters": n}
