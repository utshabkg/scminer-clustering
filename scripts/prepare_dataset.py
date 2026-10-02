"""Preprocess a dataset and save it for clustering.

Usage: python scripts/prepare_dataset.py zeisel
Writes <processed>/<name>.h5ad and a MICA input file in the same folder.
"""
import sys

from scipy import sparse

from scminer_clust import data, preprocess
from scminer_clust.paths import PROCESSED


def main(name):
    adata = data.load(name)
    print(f"loaded {name}: {adata.n_obs} cells x {adata.n_vars} genes")
    adata = preprocess.run(adata)
    print(f"after quality control: {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs["label"].value_counts().to_string())

    PROCESSED.mkdir(parents=True, exist_ok=True)
    out = PROCESSED / f"{name}.h5ad"
    adata.write_h5ad(out)
    print(f"wrote {out}")

    # The authors' MICA runs on anndata 0.8, which cannot read files written by
    # current anndata. Write a plain exchange copy that scripts/export_for_mica.py
    # turns into an h5ad inside the mica-legacy environment.
    exch = PROCESSED / f"{name}_exchange"
    exch.mkdir(exist_ok=True)
    sparse.save_npz(exch / "X.npz", sparse.csr_matrix(adata.X))
    (exch / "cells.txt").write_text("\n".join(adata.obs_names))
    (exch / "genes.txt").write_text("\n".join(adata.var_names))
    adata.obs[["label", "label_fine"]].to_csv(exch / "labels.csv")
    print(f"wrote exchange files to {exch}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
