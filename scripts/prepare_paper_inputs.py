"""Save the authors' own processed MICA inputs in our format.

Usage: python scripts/prepare_paper_inputs.py
Writes data/processed/paper_<name>.h5ad and paper_<name>_exchange/ for each
benchmark dataset, with the authors' values left untouched, so MICA can be run on
exactly what the paper ran it on.
"""
import sys

from scipy import sparse

from scminer_clust import data
from scminer_clust.paths import PROCESSED


def main():
    for name in data.PAPER_INPUTS:
        a = data.load_paper_input(name)
        out = f"paper_{name}"
        a.write_h5ad(PROCESSED / f"{out}.h5ad")
        exch = PROCESSED / f"{out}_exchange"
        exch.mkdir(exist_ok=True)
        sparse.save_npz(exch / "X.npz", sparse.csr_matrix(a.X))
        (exch / "cells.txt").write_text("\n".join(a.obs_names))
        (exch / "genes.txt").write_text("\n".join(a.var_names))
        a.obs[["label", "label_fine"]].to_csv(exch / "labels.csv")
        # the processed file is what the baselines and gene selection read
        assert a.obs["label"].notna().all() and not a.obs["label"].isin(["nan", "None"]).any(), out
        print(f"{out}: {a.n_obs} cells x {a.n_vars} genes, {a.obs.label.nunique()} classes from {a.uns['label_column']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
