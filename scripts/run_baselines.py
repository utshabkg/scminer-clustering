"""Run the baseline clustering methods on a prepared dataset.

Usage: python scripts/run_baselines.py zeisel [n_seeds] [n_hvg]

n_hvg=2000 reproduces the paper's baseline setup. Leave it out to use every gene.
The gene selection method comes from the HVG_FLAVOR environment variable (default seurat).

Each method is run with the number of clusters set to the true number of classes,
which is the protocol the scMINER paper used. The label-free setting comes later.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

from scminer_clust import cluster, metrics
from scminer_clust.paths import DATA, PROCESSED


def exp_name(name, n_hvg):
    exp = f"baselines_{name}" + (f"_hvg{n_hvg}" if n_hvg else "")
    if n_hvg and cluster.HVG_FLAVOR != "seurat":
        exp += f"_{cluster.HVG_FLAVOR}"
    return exp


def main(name, n_seeds=3, n_hvg=None):
    cluster.N_HVG = n_hvg
    cluster.HVG_FLAVOR = os.environ.get("HVG_FLAVOR", "seurat")
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    truth = adata.obs["label"].to_numpy()
    k = len(np.unique(truth))
    print(f"{name}: {adata.n_obs} cells, {adata.n_vars} genes, {k} reference classes")

    rows = []
    for seed in range(n_seeds):
        # Ward hierarchical clustering needs a full cell-by-cell matrix, which does
        # not fit beyond a few tens of thousands of cells, so it is skipped there.
        methods = ("kmeans_pca", "hierarchical_pca", "leiden_pca") if adata.n_obs <= 20_000 else ("kmeans_pca", "leiden_pca")
        if seed == 0 and len(methods) == 2:
            print(f"  skipping hierarchical: {adata.n_obs} cells is too many for a full distance matrix")
        for method in methods:
            rec = {"dataset": name, "method": method, "seed": seed, "k_setting": "paper_protocol",
                   "genes": n_hvg or "all", "hvg_flavor": cluster.HVG_FLAVOR if n_hvg else ""}
            with metrics.timed(rec):
                if method == "kmeans_pca":
                    labels, emb = cluster.kmeans_pca(adata, k, seed=seed)
                elif method == "hierarchical_pca":
                    labels, emb = cluster.hierarchical_pca(adata, k, seed=seed)
                else:
                    labels, emb, info = cluster.leiden_paper_protocol(adata, truth, seed=seed)
                    rec["resolution"] = info["resolution"]
                    rec["n_resolutions_hit_k"] = info["n_hits"]
                    rec["median_ari_of_hits"] = info["median_ari_of_hits"]
            rec.update(metrics.score(truth, labels, embedding=emb))
            runs = DATA / "runs" / exp_name(name, n_hvg)
            runs.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"label": labels}, index=adata.obs_names).to_csv(runs / f"{method}_seed{seed}.csv")
            rec["n_cells"], rec["n_genes"] = adata.n_obs, (n_hvg or adata.n_vars)
            rows.append(rec)
            print(f"  seed {seed} {method:18s} ARI {rec['ari']:.3f}  AMI {rec['ami']:.3f}  "
                  f"clusters {rec['n_clusters']}  {rec['seconds']}s")

    exp = exp_name(name, n_hvg)
    d = metrics.write(exp, rows, {"dataset": name, "n_seeds": n_seeds, "true_k": int(k), "n_hvg": n_hvg,
                                  "hvg_flavor": cluster.HVG_FLAVOR if n_hvg else None,
                                  "n_cells": adata.n_obs, "n_genes_total": adata.n_vars})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    n_hvg = int(sys.argv[3]) if len(sys.argv) > 3 else None
    sys.exit(main(sys.argv[1], n_seeds, n_hvg))
