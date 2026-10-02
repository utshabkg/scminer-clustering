"""Write exchange folders holding only the 2,000 highly variable genes, for MICA.

Usage: python scripts/export_hvg_subsets.py zeisel [n_hvg]
Writes <name>_hvg<n>_<flavor>_exchange for each of the three selection methods,
using the same select_hvg() the baselines use, so MICA and the baselines see the
same genes. Then run export_for_mica.py on each folder name inside mica-legacy.
"""
import shutil
import sys

import anndata as ad
from scipy import sparse

from scminer_clust import cluster
from scminer_clust.paths import PROCESSED

FLAVORS = ("seurat", "cell_ranger", "seurat_v3")


def main(name, n_hvg=2000, flavors=FLAVORS):
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    for flavor in flavors:
        genes = cluster.select_hvg(adata, n_hvg, flavor)
        sub = adata[:, genes]
        exch = PROCESSED / f"{name}_hvg{n_hvg}_{flavor}_exchange"
        exch.mkdir(exist_ok=True)
        sparse.save_npz(exch / "X.npz", sparse.csr_matrix(sub.X))
        (exch / "cells.txt").write_text("\n".join(sub.obs_names))
        (exch / "genes.txt").write_text("\n".join(sub.var_names))
        shutil.copy(PROCESSED / f"{name}_exchange" / "labels.csv", exch / "labels.csv")
        print(f"{name} {flavor}: {sub.n_vars} genes -> {exch.name}")
    return 0


if __name__ == "__main__":
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    sys.exit(main(sys.argv[1], n, tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else FLAVORS))
