"""A clean sheet for the team's pivot table (Rhea's R-1).

Usage: python scripts/export_pivot.py   (after update_tracker.py)
docs/results_tracker.csv holds every row from every experiment, including
label-free runs, repeats and diagnostics. This keeps only the main comparison:
each method told the true number of clusters (the paper's rule), one row per
run, with plain column names. Writes docs/results_for_pivot.csv.
"""
import sys

import pandas as pd

from scminer_clust.paths import REPO

METHODS = {"mica": "MICA", "leiden_pca": "Leiden (Scanpy)", "kmeans_pca": "k-means", "hierarchical_pca": "hierarchical"}
GENES = {"all genes": "all genes", "2000 genes (seurat)": "2000 HVG seurat",
         "2000 genes (cell_ranger)": "2000 HVG cell_ranger", "2000 genes (seurat_v3)": "2000 HVG seurat_v3"}


def main():
    t = pd.read_csv(REPO / "docs" / "results_tracker.csv")
    # Graph-mode MICA (_ge) and repeat runs are reported on their own, not in the pivot.
    keep = (t.k_setting == "paper_protocol") & t.method.isin(METHODS) \
        & ~t.experiment.str.endswith("_ge") & ~t.experiment.str.startswith("repeat_")
    t = t[keep].copy()
    out = pd.DataFrame({
        "dataset": t.dataset.str.removeprefix("paper_"),
        "input": t.dataset.map(lambda d: "authors' data" if d.startswith("paper_") else "our preprocessing"),
        "genes": t.gene_setting.map(GENES),
        "method": t.method.map(METHODS),
        "seed": t.seed,
        "clusters_found": t.n_clusters,
        "true_clusters": t.n_true_classes,
        "ari": t.ari.round(3),
        "ami": t.ami.round(3),
        "experiment": t.experiment,
    }).sort_values(["input", "dataset", "genes", "method", "seed"])
    dup = out.duplicated(["dataset", "input", "genes", "method", "seed"], keep=False)
    if dup.any():
        print("duplicate runs, check before sharing:")
        print(out[dup].to_string(index=False))
        return 1
    f = REPO / "docs" / "results_for_pivot.csv"
    out.to_csv(f, index=False)
    n = out.groupby(["input", "dataset", "genes", "method"]).ngroups
    print(f"wrote {f}: {len(out)} runs, {n} dataset, gene and method combinations")
    return 0


if __name__ == "__main__":
    sys.exit(main())
