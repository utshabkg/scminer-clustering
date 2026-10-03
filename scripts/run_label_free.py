"""Baselines with the number of clusters chosen without labels (U-12, U-13).

Usage: python scripts/run_label_free.py paper_zeisel [n_seeds] [n_hvg]
k-means and hierarchical try k = 2..20 and keep the best silhouette on the
principal components. Leiden tries the full resolution sweep and keeps the best
silhouette, and is also reported at scanpy's default resolution of 1.0. Same
gene settings as run_baselines.py (HVG_FLAVOR picks the selection method).
"""
import os
import sys

import anndata as ad
import numpy as np
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering, KMeans

from scminer_clust import cluster, metrics, selection
from scminer_clust.paths import PROCESSED


def main(name, n_seeds=3, n_hvg=None):
    cluster.N_HVG = n_hvg
    cluster.HVG_FLAVOR = os.environ.get("HVG_FLAVOR", "seurat")
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    truth = adata.obs["label"].to_numpy()
    ks = range(2, selection.K_MAX + 1)
    rows = []
    for seed in range(n_seeds):
        emb = cluster.pca(adata, seed=seed)
        cands = {
            "kmeans_pca": {k: KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(emb).astype(str) for k in ks},
            "hierarchical_pca": {k: AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(emb).astype(str) for k in ks},
        }
        g, _ = cluster._graph(adata, seed=seed)
        leiden = {}
        for res in cluster.PAPER_RESOLUTIONS:
            sc.tl.leiden(g, resolution=float(res), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
            leiden[float(res)] = g.obs["leiden"].to_numpy().astype(str)
        cands["leiden_pca"] = leiden
        sc.tl.leiden(g, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        default = g.obs["leiden"].to_numpy().astype(str)
        for method, parts in cands.items():
            key, info = selection.pick_by_silhouette(emb, parts, seed)
            rec = {"dataset": name, "method": method, "seed": seed, "k_setting": "silhouette",
                   "genes": n_hvg or "all", "hvg_flavor": cluster.HVG_FLAVOR if n_hvg else "", **info}
            rec.update(metrics.score(truth, parts[key]))
            rows.append(rec)
        rec = {"dataset": name, "method": "leiden_pca", "seed": seed, "k_setting": "default_resolution_1",
               "genes": n_hvg or "all", "hvg_flavor": cluster.HVG_FLAVOR if n_hvg else ""}
        rec.update(metrics.score(truth, default))
        rows.append(rec)
        print(f"  seed {seed}: " + ", ".join(f"{r['method']}/{r['k_setting'][:4]} k={r['n_clusters']} ARI {r['ari']:.3f}" for r in rows[-4:]))
    exp = f"labelfree_{name}" + (f"_hvg{n_hvg}" if n_hvg else "") + (f"_{cluster.HVG_FLAVOR}" if n_hvg and cluster.HVG_FLAVOR != "seurat" else "")
    d = metrics.write(exp, rows, {"dataset": name, "rule": "max silhouette over k=2..20 or the resolution sweep",
                                  "n_hvg": n_hvg, "hvg_flavor": cluster.HVG_FLAVOR if n_hvg else None})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    n_hvg = int(sys.argv[3]) if len(sys.argv) > 3 else None
    sys.exit(main(sys.argv[1], n_seeds, n_hvg))
