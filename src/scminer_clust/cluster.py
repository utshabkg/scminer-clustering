"""Baseline clustering methods, all on the same preprocessed data."""
import numpy as np
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering, KMeans


N_HVG = None  # set by run scripts; None means use every gene
HVG_FLAVOR = "seurat"  # "seurat" and "cell_ranger" use log data, "seurat_v3" uses raw counts


def select_hvg(adata, n_hvg, flavor):
    """Names of the n_hvg highly variable genes, picked with the given method.

    Used by both the baselines and the MICA export, so both see exactly the same
    genes. "seurat" and "cell_ranger" work on the log data, "seurat_v3" on the
    raw values in the counts layer.
    """
    a = adata.copy()
    # Our matrix is log2(CPM + 1). Without this, the seurat method treats it as
    # natural log and picks 19 to 47 percent different genes (found in the Oct 2 audit).
    a.uns["log1p"] = {"base": 2}
    if flavor == "seurat_v3":
        sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    else:
        sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor=flavor)
    return a.var_names[a.var["highly_variable"]].tolist()


def pca(adata, n_comps=50, seed=0, n_hvg=None):
    """Principal components on the log normalized matrix, cached in .obsm.

    n_hvg keeps only that many highly variable genes first. The scMINER paper ran
    every baseline on 2,000 of them, so we need this to reproduce their setup.
    """
    n_hvg = n_hvg if n_hvg is not None else N_HVG
    key = f"X_pca_{n_comps}_{seed}_{n_hvg or 'all'}_{HVG_FLAVOR if n_hvg else ''}"
    if key not in adata.obsm:
        a = adata.copy()
        if n_hvg:
            a = a[:, select_hvg(adata, n_hvg, HVG_FLAVOR)].copy()
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


# The paper's text says the baselines' Leiden resolution was swept from 0.1 to 4.0
# in steps of 0.1, but the authors' own Scanpy code (scMINER_analysis,
# Clustering_HPC/Scanpy/Scanpy_clustering.py) sweeps 0.01 to 2.0 in steps of 0.01.
# The coarse grid often cannot hit the true number of clusters (found Oct 2), so
# we follow the code. MICA's own sweep keeps its 0.1 grid, as in the paper.
PAPER_RESOLUTIONS = np.round(np.arange(0.01, 2.0 + 1e-9, 0.01), 2)


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
    hits = [-s[1] for s in scored if s[0] == 0]
    # How many settings hit the true count, and their median ARI, show whether the
    # pick is a best-of-many (an advantage) or the only option.
    return res, {"resolution": float(res), "n_clusters": int(n), "hit_true_k": gap == 0,
                 "n_hits": len(hits), "median_ari_of_hits": float(np.median(hits)) if hits else None}


def leiden_paper_protocol(adata, truth, n_comps=50, n_neighbors=15, seed=0):
    """Leiden swept over the paper's resolution grid, picked by the paper's rule."""
    a, emb = _graph(adata, n_comps, n_neighbors, seed)
    partitions = {}
    for res in PAPER_RESOLUTIONS:
        sc.tl.leiden(a, resolution=float(res), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        partitions[float(res)] = a.obs["leiden"].to_numpy().astype(str)
    res, info = pick_paper_protocol(partitions, truth)
    return partitions[res], emb, info
