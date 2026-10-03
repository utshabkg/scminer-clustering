"""How much does Leiden on tiny datasets depend on the number of principal components?

Usage: python scripts/pc_sensitivity.py
The paper's baselines use 50 PCs. On datasets of a few hundred cells that may
be too many, so this reruns Leiden under the paper's rule with 5 to 50 PCs, for
every gene setting, on the authors' Buettner and Chung inputs. Writes
results/pc_sensitivity/.
"""
import sys

import anndata as ad
import numpy as np

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED

DATASETS = ("buettner", "chung")
N_PCS = (5, 10, 20, 30, 50)
FLAVORS = (None, "seurat", "cell_ranger", "seurat_v3")


def main():
    rows = []
    for ds in DATASETS:
        a = ad.read_h5ad(PROCESSED / f"paper_{ds}.h5ad")
        truth = a.obs["label"].to_numpy()
        for flavor in FLAVORS:
            cluster.N_HVG = 2000 if flavor else None
            cluster.HVG_FLAVOR = flavor or "seurat"
            for nc in N_PCS:
                for seed in range(3):
                    labels, emb, info = cluster.leiden_paper_protocol(a, truth, n_comps=nc, seed=seed)
                    r = {"dataset": f"paper_{ds}", "method": "leiden_pca", "genes": 2000 if flavor else "all",
                         "hvg_flavor": flavor or "", "n_pcs": nc, "seed": seed, "k_setting": "pc_sensitivity"}
                    r.update(metrics.score(truth, labels))
                    rows.append(r)
            line = " | ".join(f"{nc} PCs {np.mean([r['ari'] for r in rows if r['dataset'] == f'paper_{ds}' and r['hvg_flavor'] == (flavor or '') and r['n_pcs'] == nc]):.2f}" for nc in N_PCS)
            print(f"{ds:9s} {flavor or 'all genes':12s} {line}")
    d = metrics.write("pc_sensitivity", rows, {"datasets": list(DATASETS), "n_pcs": list(N_PCS), "seeds": 3,
                                               "rule": "paper protocol, resolutions 0.01 to 2.0"})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
