"""Which cluster counts a Leiden resolution grid reaches on the full data, with no labels.

Usage: python scripts/probe_grid_coverage.py <grid_name> [seeds]

The stability comparison caps both methods at 2 to 20 clusters, but a cap is not a match:
on the old 0.05 to 2.0 grid the baseline never produced more than 8 clusters on Yan or
fewer than 4 on Pollen, while MICA's k grid reaches every count. This sweeps a candidate
grid on each dataset's full data, exactly as the stability runners build their graphs,
and records the count each resolution gives. It reads no labels and scores nothing, so it
can be used to fix the grid before the stability runs.

Writes results/grid_coverage_<grid_name>/ with one row per dataset, seed and resolution,
and a summary of the counts from 2 to 20 each dataset never reaches.
"""
import sys

import anndata as ad
import numpy as np
import scanpy as sc

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED

DATASETS = ["paper_yan", "paper_goolam", "paper_buettner", "paper_pollen", "paper_chung", "paper_usoskin",
            "paper_kolodziejczyk", "paper_klein", "paper_zeisel", "paper_pbmc14k"]
GRIDS = {
    "old": np.round(np.arange(0.05, 2.0 + 1e-9, 0.05), 2),
    "fine": np.round(np.geomspace(0.001, 50, 300), 5),
}


def main(grid_name, seeds="0,1,2"):
    grid = GRIDS[grid_name]
    rows, summary = [], []
    for name in DATASETS:
        a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
        for seed in [int(s) for s in seeds.split(",")]:
            g, _ = cluster._graph(a.copy(), seed=seed)  # a fresh copy, as the runners make per run
            counts = {}
            for r in grid:
                sc.tl.leiden(g, resolution=float(r), random_state=seed, flavor="igraph", n_iterations=2,
                             directed=False)
                counts[float(r)] = int(g.obs["leiden"].nunique())
                rows.append({"dataset": name, "seed": seed, "resolution": float(r), "n_clusters": counts[float(r)]})
            reached = set(counts.values())
            missing = [k for k in range(2, 21) if k not in reached]
            summary.append({"dataset": name, "seed": seed, "n_cells": a.n_obs, "grid": grid_name,
                            "n_resolutions": len(grid), "min_clusters": min(reached), "max_clusters": max(reached),
                            "missing_2_to_20": " ".join(map(str, missing))})
            print(f"{name} seed {seed}: {min(reached)} to {max(reached)} clusters, missing {missing}", flush=True)
    d = metrics.write(f"grid_coverage_{grid_name}", rows,
                      {"grid": grid_name, "resolutions": grid.tolist(), "seeds": seeds,
                       "graph": "cluster._graph, as in run_rule_stability.py", "labels_used": False})
    import pandas as pd
    pd.DataFrame(summary).to_csv(d / "summary.csv", index=False)
    print(f"wrote {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
