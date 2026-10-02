"""Run the baseline clustering methods on a prepared dataset.

Usage: python scripts/run_baselines.py zeisel [n_seeds] [n_hvg]

n_hvg=2000 reproduces the paper's baseline setup. Leave it out to use every gene.

Each method is run with the number of clusters set to the true number of classes,
which is the protocol the scMINER paper used. The label-free setting comes later.
"""
import sys

import anndata as ad
import numpy as np

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED


def main(name, n_seeds=3, n_hvg=None):
    cluster.N_HVG = n_hvg
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    truth = adata.obs["label"].to_numpy()
    k = len(np.unique(truth))
    print(f"{name}: {adata.n_obs} cells, {adata.n_vars} genes, {k} reference classes")

    rows = []
    for seed in range(n_seeds):
        for method in ("kmeans_pca", "hierarchical_pca", "leiden_pca"):
            rec = {"dataset": name, "method": method, "seed": seed, "k_setting": "paper_protocol",
                   "genes": n_hvg or "all"}
            with metrics.timed(rec):
                if method == "kmeans_pca":
                    labels, emb = cluster.kmeans_pca(adata, k, seed=seed)
                elif method == "hierarchical_pca":
                    labels, emb = cluster.hierarchical_pca(adata, k, seed=seed)
                else:
                    labels, emb, info = cluster.leiden_paper_protocol(adata, truth, seed=seed)
                    rec["resolution"] = info["resolution"]
            rec.update(metrics.score(truth, labels, embedding=emb))
            rows.append(rec)
            print(f"  seed {seed} {method:18s} ARI {rec['ari']:.3f}  AMI {rec['ami']:.3f}  "
                  f"clusters {rec['n_clusters']}  {rec['seconds']}s")

    exp = f"baselines_{name}" + (f"_hvg{n_hvg}" if n_hvg else "")
    d = metrics.write(exp, rows, {"dataset": name, "n_seeds": n_seeds, "true_k": int(k), "n_hvg": n_hvg})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    n_hvg = int(sys.argv[3]) if len(sys.argv) > 3 else None
    sys.exit(main(sys.argv[1], n_seeds, n_hvg))
