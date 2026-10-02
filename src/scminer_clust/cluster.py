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


def _graph(adata, n_comps=50, n_neighbors=15, seed=0):
    """Neighbor graph on the principal components, held in a small AnnData.

    Only the embedding is carried, not the expression matrix, so a resolution
    sweep does not copy the data once per step.
    """
    import anndata as ad

    emb = pca(adata, n_comps, seed)
    a = ad.AnnData(obs=adata.obs[[]].copy())
    a.obsm["X_pca"] = emb
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_pca", random_state=seed)
    return a, emb


def leiden_pca(adata, resolution=1.0, n_comps=50, n_neighbors=15, seed=0):
    """The standard Scanpy pipeline: principal components, neighbor graph, Leiden."""
    a, emb = _graph(adata, n_comps, n_neighbors, seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return a.obs["leiden"].to_numpy().astype(str), emb


PAPER_RESOLUTIONS = np.round(np.arange(0.1, 4.0 + 1e-9, 0.1), 1)


def pick_paper_protocol(partitions, truth):
    """The scMINER benchmark rule, which uses the labels twice.

    Among partitions whose cluster count equals the true number of types, keep
    the one with the highest ARI. If none hits it exactly, keep the closest count,
    breaking ties by ARI. `partitions` maps a resolution to its labels.
    """
    from sklearn.metrics import adjusted_rand_score

    k = len(np.unique(truth))
    scored = []
    for res, labels in partitions.items():
        n = len(np.unique(labels))
        scored.append((abs(n - k), -adjusted_rand_score(truth, labels), res, n))
    gap, neg_ari, res, n = min(scored)
    return res, {"resolution": float(res), "n_clusters": int(n), "hit_true_k": gap == 0}


def leiden_paper_protocol(adata, truth, n_comps=50, n_neighbors=15, seed=0):
    """Leiden swept over the paper's resolution grid, picked by the paper's rule."""
    a, emb = _graph(adata, n_comps, n_neighbors, seed)
    partitions = {}
    for res in PAPER_RESOLUTIONS:
        sc.tl.leiden(a, resolution=float(res), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        partitions[float(res)] = a.obs["leiden"].to_numpy().astype(str)
    res, info = pick_paper_protocol(partitions, truth)
    return partitions[res], emb, info
