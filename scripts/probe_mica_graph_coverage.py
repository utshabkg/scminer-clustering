"""Which cluster counts MICA's graph-mode Louvain grid reaches, with no labels.

Usage: python scripts/probe_mica_graph_coverage.py [max_resolution] [workers]

MICA's own grid runs Louvain from 0.1 to 4.0, and on PBMC14K it never gave fewer than 5
clusters, so MICA could not be offered the coarse counts the baseline reaches. In MICA's
Louvain a higher resolution gives fewer clusters. This runs MICA's graph pipeline once on
the full PBMC14K data and once on each of the three 2,500 cell pools of the mode test,
sweeping 0.1 to max_resolution, and records the count at each resolution. Nothing is
scored, so it can fix the grid before the matched stability runs.

Writes results/grid_coverage_mica_graph/.
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import run_rule_mode_size as mode_size
import run_rule_stability_graph as graph
from scminer_clust import metrics
from scminer_clust.paths import DATA, PROCESSED


def main(max_res=30.0, workers=8):
    graph.MICA_MAX_RES = float(max_res)
    graph.WORK.mkdir(parents=True, exist_ok=True)
    n = ad.read_h5ad(PROCESSED / f"{graph.NAME}.h5ad", backed="r").n_obs
    D = np.load(DATA / "runs" / "rules" / f"{graph.NAME}_ge_distance.npy", mmap_mode="r")
    inputs = [("full", np.arange(n))] + [(f"pool2500_seed{s}", mode_size.pool_cells(n, 2500, s)) for s in (0, 1, 2)]
    rows, summary = [], []
    for name, cells in inputs:
        # the pools index the full matrix directly, as the mode test does
        parts = graph.mica_graph_partitions(D, cells, f"probe_{name}", workers=int(workers))
        counts = {r: len(np.unique(lab)) for r, lab in parts.items()}
        rows += [{"input": name, "n_cells": len(cells), "resolution": r, "n_clusters": k} for r, k in sorted(counts.items())]
        reached = set(counts.values())
        missing = [k for k in range(2, 21) if k not in reached]
        first2 = min((r for r, k in counts.items() if k <= 2), default=None)
        summary.append({"input": name, "n_cells": len(cells), "max_resolution": float(max_res),
                        "min_clusters": min(reached), "max_clusters": max(reached),
                        "first_resolution_at_2_or_fewer": first2, "missing_2_to_20": " ".join(map(str, missing))})
        print(f"{name}: {min(reached)} to {max(reached)} clusters, first at 2 or fewer: {first2}, "
              f"missing {missing}", flush=True)
    d = metrics.write("grid_coverage_mica_graph", rows,
                      {"max_resolution": float(max_res), "step": 0.1, "labels_used": False,
                       "pipeline": "run_rule_stability_graph.mica_graph_partitions"})
    pd.DataFrame(summary).to_csv(d / "summary.csv", index=False)
    print(f"wrote {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
